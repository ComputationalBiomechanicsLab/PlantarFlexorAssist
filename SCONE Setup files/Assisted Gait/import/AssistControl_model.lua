function init(model, par, side)
	
	-- =========================
	-- Muscles & actuators
	-- =========================
	gas_r = model:find_muscle("gastroc_r")
	gas_l = model:find_muscle("gastroc_l")
	
	sol_r = model:find_muscle("soleus_r")
	sol_l = model:find_muscle("soleus_l")
	
	p_motor_Str = model:find_actuator("ankle_angle_r")
	p_motor_Stl = model:find_actuator("ankle_angle_l")
	
	-- =========================
	-- Initial time 
	-- =========================
	t0  = 0.0
	t00 = 0.0
	
	-- =========================
	-- Controller gains
	-- =========================
	KpGas_ES = par:create_from_mean_std("KpGas_ES", 0,   0.5, 0, 0.5)
	KpSol_ES = par:create_from_mean_std("KpSol_ES", 0,   0.5, 0, 0.5)
	
	KpGas_LS = par:create_from_mean_std("KpGas_LS", 0.3, 0.5, 0, 1)
	KpSol_LS = par:create_from_mean_std("KpSol_LS", 0.5, 0.5, 0, 1)
	
	KpGas_LO = par:create_from_mean_std("KpGas_LO", 0.4, 0.5, 0, 1)
	KpSol_LO = par:create_from_mean_std("KpSol_LO", 0.5, 0.5, 0, 1)
	
	KpGas_S  = par:create_from_mean_std("KpGas_S",  0,   0.5, 0, 0.5)
	KpSol_S  = par:create_from_mean_std("KpSol_S",  0,   0.5, 0, 0.5)
	
	KpGas_L  = par:create_from_mean_std("KpGas_L",  0,   0.5, 0, 0.5)
	KpSol_L  = par:create_from_mean_std("KpSol_L",  0,   0.5, 0, 0.5)
	
end


function update(model)
	
	local t = model:time()
	t = model:time()
	tr = model:time()
	tl = model:time()
	tb = model:time()
	tbl = model:time()
	
		-- =========================
	-- Muscle forces
	-- =========================
	gas_r_F = gas_r:force()
	gas_l_F = gas_l:force()
	
	sol_r_F = sol_r:force()
	sol_l_F = sol_l:force()
	
	-- =========================
	-- Contact model
	-- =========================
	calcnr = model:find_body("calcn_r")
	calcnl = model:find_body("calcn_l")
	
	GRFr = calcnr:contact_force().y
	GRFl = calcnl:contact_force().y
	
	calcnr_y = calcnr:com_pos().y
	calcnl_y = calcnl:com_pos().y
	
	toesr = model:find_body("toes_r")
	toesl = model:find_body("toes_l")
	
	toesr_y = toesr:com_pos().y
	toesl_y = toesl:com_pos().y
	
	-- =========================
	-- Phase timers
	-- =========================
	if calcnr_y < toesr_y and GRFl > GRFr and GRFr > 90 and GRFr < 180 and GRFl > 500 then
		tr = 0
	end
	
	if tb - tr > 0.01 then
		t0 = tb
	end
	tb = tb - t0
	
	if calcnl_y < toesl_y and GRFr > GRFl and GRFl > 90 and GRFl < 180 and GRFr > 500 then
		tl = 0
	end
	
	if tbl - tl > 0.01 then
		t00 = tbl
	end
	tbl = tbl - t00
	
	-- =========================
	-- Phase angle (PHI)
	-- =========================
	local p1,p2,p3,p4,p5,p6 = -1.986, 2.101, 4.98, -8.596, 8.72, 0.02334
	
	PHI1  = (p1*tb^5 + p2*tb^4 + p3*tb^3 + p4*tb^2 + p5*tb + p6) / (2*3.14)
	PHI1l = (p1*tbl^5 + p2*tbl^4 + p3*tbl^3 + p4*tbl^2 + p5*tbl + p6) / (2*3.14)
	
	-- =========================
	-- Fourier series (GRF normalization)
	-- =========================
	local a0,a1,b1 = 624.4, -526.1, 320.9
	local a2,b2 = 160.4, 47.45
	local a3,b3 = -112.9, -113.8
	local a4,b4 = 5.281, 31.15
	local a5,b5 = 30.31, -1.966
	local a6,b6 = 16.73, -20.21
	local a7,b7 = -6.167, 20.18
	local a8,b8 = 10.33, 1.37
	local w = 5.821
	
	local function fourier(x)
	return a0 +
	a1*math.cos(x*w) + b1*math.sin(x*w) +
	a2*math.cos(2*x*w) + b2*math.sin(2*x*w) +
	a3*math.cos(3*x*w) + b3*math.sin(3*x*w) +
	a4*math.cos(4*x*w) + b4*math.sin(4*x*w) +
	a5*math.cos(5*x*w) + b5*math.sin(5*x*w) +
	a6*math.cos(6*x*w) + b6*math.sin(6*x*w) +
	a7*math.cos(7*x*w) + b7*math.sin(7*x*w) +
	a8*math.cos(8*x*w) + b8*math.sin(8*x*w)
end

GFNormr = fourier(PHI1)
GFNorml = fourier(PHI1l)

-- =========================
-- Torque terms
-- =========================
TGasNr = -0.0165 * GFNormr
TGasNl = -0.0165 * GFNorml

TSolNr = -0.0412 * GFNormr
TSolNl = -0.0412 * GFNorml

TGasPFr = -0.0392 * gas_r_F
TGasPFl = -0.0392 * gas_l_F

TSolPFr = -0.0364 * sol_r_F
TSolPFl = -0.0364 * sol_l_F

-- =========================
-- Gain scheduling (right/left)
-- =========================
local function select_gain(phi, ES, LS, LO, S, L)
if phi < 0.279 then return ES, ES
elseif phi < 0.540 then return LS, LS
elseif phi < 0.586 then return LO, LO
elseif phi < 0.786 then return S,  S
else return L, L end
end

KpGasr, KpSolr = select_gain(PHI1,
KpGas_ES, KpGas_LS, KpGas_LO, KpGas_S, KpGas_L)

KpGasl, KpSoll = select_gain(PHI1l,
KpGas_ES, KpGas_LS, KpGas_LO, KpGas_S, KpGas_L)

-- =========================
-- Control law
-- =========================
GTorquer = KpGasr * (TGasNr - TGasPFr)
STorquer = KpSolr * (TSolNr - TSolPFr)

GTorquel = KpGasl * (TGasNl - TGasPFl)
STorquel = KpSoll * (TSolNl - TSolPFl)

Torquer = GTorquer + STorquer
Torquel = GTorquel + STorquel

-- =========================
-- Actuation
-- =========================
if t > 2 then
	p_motor_Str:add_input(Torquer)
	p_motor_Stl:add_input(Torquel)
end

end


function store_data(frame)
frame:set_value("Torquer", Tr)
frame:set_value("Torquel", Tl)
end