function init(model, par, side)

	-- =========================
	-- Muscles & Actuators
	-- =========================
	gas_r = model:find_muscle("gastroc_r")
	gas_l = model:find_muscle("gastroc_l")
	sol_r = model:find_muscle("soleus_r")
	sol_l = model:find_muscle("soleus_l")

	p_motor_Str = model:find_actuator("ankle_angle_r")
	p_motor_Stl = model:find_actuator("ankle_angle_l")

	-- =========================
	-- Timers
	-- =========================
	t0  = 0.0
	t00 = 0.0

	-- =========================
	-- Gains
	-- =========================
	KpGas_ES = par:create_from_mean_std("KpGas_ES", 0.01, 0.5, 0, 0.5)
	KpSol_ES = par:create_from_mean_std("KpSol_ES", 0.01, 0.5, 0, 0.5)

	KpGas_LS = par:create_from_mean_std("KpGas_LS", 0.04, 0.1, 0, 1)
	KpSol_LS = par:create_from_mean_std("KpSol_LS", 0.08, 0.1, 0, 1)

	KpGas_LO = par:create_from_mean_std("KpGas_LO", 0.1, 0.1, 0, 1)
	KpSol_LO = par:create_from_mean_std("KpSol_LO", 0.15, 0.2, 0, 1)

	KpGas_S  = par:create_from_mean_std("KpGas_S", 0.001, 0.01, 0, 0.01)
	KpSol_S  = par:create_from_mean_std("KpSol_S", 0.001, 0.001, 0, 0.01)

	KpGas_L  = par:create_from_mean_std("KpGas_L", 0.01, 0.01, 0, 0.01)
	KpSol_L  = par:create_from_mean_std("KpSol_L", 0.01, 0.001, 0, 0.01)

end


-- =========================================================
-- Fourier series 
-- =========================================================
local function fourier(x, a, b, w)
	return a[1]
	+ a[2]*math.cos(x*w) + b[1]*math.sin(x*w)
	+ a[3]*math.cos(2*x*w) + b[2]*math.sin(2*x*w)
	+ a[4]*math.cos(3*x*w) + b[3]*math.sin(3*x*w)
	+ a[5]*math.cos(4*x*w) + b[4]*math.sin(4*x*w)
	+ a[6]*math.cos(5*x*w) + b[5]*math.sin(5*x*w)
	+ a[7]*math.cos(6*x*w) + b[6]*math.sin(6*x*w)
	+ a[8]*math.cos(7*x*w) + b[7]*math.sin(7*x*w)
	+ a[9]*math.cos(8*x*w) + b[8]*math.sin(8*x*w)
end


-- =========================================================
-- Gain scheduling 
-- =========================================================
local function select_gain(phi, ES, LS, LO, S, L)
	if phi < 0.24 then return ES, ES
	elseif phi < 0.46 then return LS, LS
	elseif phi < 0.56 then return LO, LO
	elseif phi < 0.74 then return S, S
	else return L, L end
end


function update(model)

	-- =========================
	-- Time
	-- =========================
	local t  = model:time()
	local dt = model:delta_time()

	local tr = t
	local tl = t
	local tb = t
	local tbl = t

	-- =========================
	-- Muscle forces
	-- =========================
	local gas_r_F = gas_r:force()
	local gas_l_F = gas_l:force()
	local sol_r_F = sol_r:force()
	local sol_l_F = sol_l:force()

	-- =========================
	-- Contact
	-- =========================
	local calcnr = model:find_body("calcn_r")
	local calcnl = model:find_body("calcn_l")

	local GRFr = calcnr:contact_force().y
	local GRFl = calcnl:contact_force().y

	local calcnr_y = calcnr:com_pos().y
	local calcnl_y = calcnl:com_pos().y

	local toesr_y = model:find_body("toes_r"):com_pos().y
	local toesl_y = model:find_body("toes_l"):com_pos().y

	-- =========================
	-- Phase reset
	-- =========================
	if calcnr_y < toesr_y and GRFl < GRFr and GRFr < 830 and GRFr > 800 and GRFl < 400 and GRFl > 350 then
		tr = 0
	end

	if calcnl_y < toesl_y and GRFr < GRFl and GRFl < 830 and GRFl > 800 and GRFr < 400 and GRFr > 350 then
		tl = 0
	end

	if tb - tr > 0.01 then t0 = tb end
	if tbl - tl > 0.01 then t00 = tbl end

	tb  = tb - t0
	tbl = tbl - t00

	-- =========================
	-- Phase (PHI)
	-- =========================
	local p = {-1.294, 2.525, 0.8035, -4.065, 6.509, 0.08953}

	local PHI1  = (p[1]*tb^5  + p[2]*tb^4  + p[3]*tb^3  + p[4]*tb^2  + p[5]*tb  + p[6]) / (2*math.pi)
	local PHI1l = (p[1]*tbl^5 + p[2]*tbl^4 + p[3]*tbl^3 + p[4]*tbl^2 + p[5]*tbl + p[6]) / (2*math.pi)

	-- =========================
	-- Fourier coefficients
	-- =========================
	local a = {493.7, -214.7, 58.88, -60.25, -20.01, -15.54, 8.291, -15.36, 12.69}
	local b = {285.3, -189.9, 137.4, -96.27, 11.51, -18.02, -6.974, -9.077}
	local w = 6.306

	local x_r = PHI1 - 0.1
	local x_l = PHI1l - 0.1

	local GFNormr = fourier(x_r, a, b, w)
	local GFNorml = fourier(x_l, a, b, w)

	-- =========================
	-- SF
	-- =========================
	local a2 = {649.3, -580.3, 239.2, -99.22, -20.09, 7.659, 15.92, -37.2, 26.16}
	local b2 = {299.7, -243.1, 188.6, -113.1, 48.59, -29.47, -0.2349, -6.094}
	local w2 = 6.473

	local SFNormr = fourier(x_r, a2, b2, w2)
	local SFNorml = fourier(x_l, a2, b2, w2)

	-- =========================
	-- Torque terms
	-- =========================
	local TGasNr = -0.043 * GFNormr
	local TGasNl = -0.043 * GFNorml

	local TGasPFr = -0.043 * gas_r_F
	local TGasPFl = -0.043 * gas_l_F

	local TSolNr = -0.040 * SFNormr
	local TSolNl = -0.040 * SFNorml

	local TSolPFr = -0.040 * sol_r_F
	local TSolPFl = -0.040 * sol_l_F

	-- =========================
	-- Gain scheduling
	-- =========================
	local KpGasr, KpSolr = select_gain(PHI1,  KpGas_ES, KpGas_LS, KpGas_LO, KpGas_S, KpGas_L)
	local KpGasl, KpSoll = select_gain(PHI1l, KpGas_ES, KpGas_LS, KpGas_LO, KpGas_S, KpGas_L)

	-- =========================
	-- Control law
	-- =========================
	local GTorquer = KpGasr * (TGasNr - TGasPFr)
	local STorquer = KpSolr * (TSolNr - TSolPFr)

	local GTorquel = KpGasl * (TGasNl - TGasPFl)
	local STorquel = KpSoll * (TSolNl - TSolPFl)

	-- saturation (same logic)
	if GTorquer > 0 then GTorquer = 0 end
	if STorquer > 0 then STorquer = 0 end
	if GTorquel > 0 then GTorquel = 0 end
	if STorquel > 0 then STorquel = 0 end

	local Torquer = GTorquer + STorquer
	local Torquel = GTorquel + STorquel

	-- =========================
	-- Actuation
	-- =========================
	if t > 2 then
		p_motor_Str:add_input(Torquer)
		p_motor_Stl:add_input(Torquel)
	end

	_last_Torquer = Torquer
	_last_Torquel = Torquel
end


function store_data(frame)
	frame:set_value("Torquer", _last_Torquer or 0)
	frame:set_value("Torquel", _last_Torquel or 0)
end