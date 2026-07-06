% This script generates the figures and analyses presented in the Results
% and Supporting Information sections of the study:
% "Quantifying the Need for Assistance in Patients with Plantar Flexor
% Muscle Weakness to Improve Gait Outcomes" by Amini et al.

% It is provided to facilitate reproducibility,
% further analyses, and adaptation for related research studies.

%% INIT
clear all;
close all;
clc
addpath(genpath(pwd))
resultDir = fullfile(pwd,'SCONE Results');
%% SCONE result file names
SconeFile(1)={'Unimpaired_gait.sto'};
SconeFile(2)={'Unimpaired_gait_model%95.sto'};
SconeFile(3)={'Unimpaired_gait_model%5.sto'};

SconeFile(4)={'PF_gait.sto'};
SconeFile(5)={'PF_gait_model%95.sto'};
SconeFile(6)={'PF_gait_model%5.sto'};

SconeFile(7)={'Assisted_gait.sto'};
SconeFile(8)={'Assisted_gait_model%95.sto'};
SconeFile(9)={'Assisted_gait_model%5.sto'};
%% Find spec data in SCONE results

NK=length(SconeFile);
aaf=strings(NK,1000);

for i=1:NK

    filename = fullfile(resultDir,SconeFile{i});
    fileID = fopen(filename);
    aa = textscan(fileID,'%s');
    fclose(fileID);
    for k=1:1000
        aaf(i,k)=aa{1,1}{k+7,1};
    end
 
    HAMS_i=find(strcmp(aaf(i,:),'hamstrings_r.mtu_force'))+1;
    BFSH_i=find(strcmp(aaf(i,:),'bifemsh_r.mtu_force'))+1;
    GMAX_i=find(strcmp(aaf(i,:),'glut_max_r.mtu_force'))+1;
    ILPSO_i=find(strcmp(aaf(i,:),'iliopsoas_r.mtu_force'))+1;
    RF_i=find(strcmp(aaf(i,:),'rect_fem_r.mtu_force'))+1;
    VAS_i=find(strcmp(aaf(i,:),'vasti_r.mtu_force'))+1; 
    GAS_i=find(strcmp(aaf(i,:),'gastroc_r.mtu_force'))+1;
    SOL_i=find(strcmp(aaf(i,:),'soleus_r.mtu_force'))+1;
    TA_i=find(strcmp(aaf(i,:),'tib_ant_r.mtu_force'))+1;

    hip_r_i=find(strcmp(aaf(i,:),'hip_flexion_r'))+1;
    knee_r_i=find(strcmp(aaf(i,:),'knee_angle_r'))+1;
    ankle_r_i=find(strcmp(aaf(i,:),'ankle_angle_r'))+1;

    hip_r_moment_i=find(strcmp(aaf(i,:),'hip_flexion_r.moment'))+1;
    knee_r_moment_i=find(strcmp(aaf(i,:),'knee_angle_r.moment'))+1;
    ankle_r_moment_i=find(strcmp(aaf(i,:),'ankle_angle_r.moment'))+1;

    grf_r_x_i=find(strcmp(aaf(i,:),'leg1_r.grf_x'))+1;
    grf_r_y_i=find(strcmp(aaf(i,:),'leg1_r.grf_y'))+1;

    com_x_i=find(strcmp(aaf(i,:),'com_x'))+1;

    assist_torque_i = find(strcmp(aaf(i,:),'Torquer'))+1; 

    temp = dlmread(cell2mat(SconeFile(i)),'',7,0);
    res(i,1:size(temp,1),1:size(temp,2)) = temp;  
    time(i,:)=res(i,:,1);

    %%%%%%%%%%%%%%%% filter data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lowpassFreq = 10;
    % Construct a 4th-order lowpass Butterworth filter based on the 'lowpassFreq' argument. 
    timeStep(i) = time(i,2) - time(i,1);
    sampleRate(i) = 1 / timeStep(i);
    halfSampleRate(i) = sampleRate(i) / 2;
    cutoffFreq(i) = lowpassFreq / halfSampleRate(i);
    [B(i,:),A(i,:)] = butter(1, cutoffFreq(i), 'low');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Muscle Forces 
    HAMS(i,:)=res(i,:,HAMS_i);HAMS(i,:)= filtfilt(B(i,:), A(i,:), HAMS(i,:));
    BFSH(i,:)=res(i,:,BFSH_i);BFSH(i,:)= filtfilt(B(i,:), A(i,:), BFSH(i,:));
    GMAX(i,:)=res(i,:,GMAX_i);GMAX(i,:)= filtfilt(B(i,:), A(i,:), GMAX(i,:));
    ILPSO(i,:)=res(i,:,ILPSO_i);ILPSO(i,:)= filtfilt(B(i,:), A(i,:), ILPSO(i,:));
    RF(i,:)=res(i,:,RF_i);RF(i,:)= filtfilt(B(i,:), A(i,:), RF(i,:));
    VAS(i,:)=res(i,:,VAS_i);VAS(i,:)= filtfilt(B(i,:), A(i,:), VAS(i,:));
    GAS(i,:)=res(i,:,GAS_i);GAS(i,:)= filtfilt(B(i,:), A(i,:), GAS(i,:));
    SOL(i,:)=res(i,:,SOL_i);SOL(i,:)= filtfilt(B(i,:), A(i,:), SOL(i,:));
    TA(i,:)=res(i,:,TA_i);TA(i,:)= filtfilt(B(i,:), A(i,:), TA(i,:));
    
    % angles    
    hip_r(i,:)=res(i,:,hip_r_i)/pi*180;hip_r(i,:)= filtfilt(B(i,:), A(i,:), hip_r(i,:));
    knee_r(i,:)=res(i,:,knee_r_i)/pi*180;knee_r(i,:)= filtfilt(B(i,:), A(i,:), knee_r(i,:));
    ankle_r(i,:)=res(i,:,ankle_r_i)/pi*180;ankle_r(i,:)= filtfilt(B(i,:), A(i,:), ankle_r(i,:));

    % Joint moments 
    hip_r_moment(i,:)=res(i,:,hip_r_moment_i);hip_r_moment(i,:)= filtfilt(B(i,:), A(i,:), hip_r_moment(i,:));
    knee_r_moment(i,:)=res(i,:,knee_r_moment_i);knee_r_moment(i,:)= filtfilt(B(i,:), A(i,:), knee_r_moment(i,:));
    ankle_r_moment(i,:)=res(i,:,ankle_r_moment_i);ankle_r_moment(i,:)= filtfilt(B(i,:), A(i,:), ankle_r_moment(i,:));

    % GRF 
    grf_r_x(i,:)=res(i,:,grf_r_x_i);grf_r_x(i,:)= filtfilt(B(i,:), A(i,:), grf_r_x(i,:));
    grf_r_y(i,:)=res(i,:,grf_r_y_i);grf_r_y(i,:)= filtfilt(B(i,:), A(i,:), grf_r_y(i,:));

    % step length 
    com_x(i,:)=res(i,:,com_x_i);com_x(i,:)= filtfilt(B(i,:), A(i,:), com_x(i,:));

    % Asssitance torque 
    if ismember(i, [7 8 9])
        assist_torque(i,:)=res(i,:,assist_torque_i);assist_torque(i,:)= filtfilt(B(i,:), A(i,:), assist_torque(i,:));
    end 
end 

%% Representative Gait Cycle for All Models
% Since the simulation data were generated over the same complete gait
% sequence for all models, a single representative gait cycle was selected
% for analysis. The cycle was cropped from one initial contact (IC) to the
% subsequent initial contact, representing a complete gait cycle.

CT_un = [366:492];CT_un_95 = [566:710];CT_un_05 = [322:436];
CT_pf = [489:616];CT_pf_95 = [417:564];CT_pf_05 = [437:549];
CT_a = [761:889];CT_a_95 = [421:567];CT_a_05 = [97:209];

CT = {CT_un,CT_un_95,CT_un_05,CT_pf,CT_pf_95,CT_pf_05,CT_a,CT_a_95,CT_a_05};


for i=1:NK
    % Crop cycle gait
    time_c{i,:} = time(i,CT{i}); time_c{i,:} = time_c{i,:} - time_c{i}(1);
    range{i,:} = linspace(time_c{i}(1),time_c{i}(end),100);

    % Muscle Forces 
    HAMS_c{i,:} = HAMS(i,CT{i});HAMS_norm{i,:} = interp1(time_c{i,:},HAMS_c{i,:},range{i,:});
    BFSH_c{i,:} = BFSH(i,CT{i});BFSH_norm{i,:} = interp1(time_c{i,:},BFSH_c{i,:},range{i,:});
    GMAX_c{i,:} = GMAX(i,CT{i});GMAX_norm{i,:} = interp1(time_c{i,:},GMAX_c{i,:},range{i,:});
    ILPSO_c{i,:} = ILPSO(i,CT{i});ILPSO_norm{i,:} = interp1(time_c{i,:},ILPSO_c{i,:},range{i,:});
    RF_c{i,:} = RF(i,CT{i});RF_norm{i,:} = interp1(time_c{i,:},RF_c{i,:},range{i,:});
    VAS_c{i,:} = VAS(i,CT{i});VAS_norm{i,:} = interp1(time_c{i,:},VAS_c{i,:},range{i,:});
    GAS_c{i,:} = GAS(i,CT{i});GAS_norm{i,:} = interp1(time_c{i,:},GAS_c{i,:},range{i,:});
    SOL_c{i,:} = SOL(i,CT{i});SOL_norm{i,:} = interp1(time_c{i,:},SOL_c{i,:},range{i,:});
    TA_c{i,:} = TA(i,CT{i});TA_norm{i,:} = interp1(time_c{i,:},TA_c{i,:},range{i,:});    
    
    
    % Angles 
    hip_r_c{i,:} = hip_r(i,CT{i});hip_r_norm{i,:} = interp1(time_c{i,:},hip_r_c{i,:},range{i,:});
    knee_r_c{i,:} = knee_r(i,CT{i});knee_r_norm{i,:} = interp1(time_c{i,:},knee_r_c{i,:},range{i,:});
    ankle_r_c{i,:} = ankle_r(i,CT{i});ankle_r_norm{i,:} = interp1(time_c{i,:},ankle_r_c{i,:},range{i,:});

     % Joint moments 
    hip_r_moment_c{i,:} = hip_r_moment(i,CT{i});hip_r_moment_norm{i,:} = interp1(time_c{i,:},hip_r_moment_c{i,:},range{i,:});
    knee_r_moment_c{i,:} = knee_r_moment(i,CT{i});knee_r_moment_norm{i,:} = interp1(time_c{i,:},knee_r_moment_c{i,:},range{i,:});
    ankle_r_moment_c{i,:} = ankle_r_moment(i,CT{i});ankle_r_moment_norm{i,:} = interp1(time_c{i,:},ankle_r_moment_c{i,:},range{i,:});

    % GRF
    grf_r_x_c{i,:} = grf_r_x(i,CT{i});grf_r_x_norm{i,:} = interp1(time_c{i,:},grf_r_x_c{i,:},range{i,:});
    grf_r_y_c{i,:} = grf_r_y(i,CT{i});grf_r_y_norm{i,:} = interp1(time_c{i,:},grf_r_y_c{i,:},range{i,:});

    % step length 
    com_x_c{i,:} = com_x(i,CT{i});step_length{i} = com_x_c{i}(end)-com_x_c{i}(1);
    com_x_norm{i,:} = interp1(time_c{i,:},com_x_c{i,:},range{i,:});

    % Gait Speed
    Gait_speed{i} = step_length{i}/time_c{i}(end);

    % Assistance torque 
    if ismember(i, [7 8 9])
        assist_torque_c{i,:} = assist_torque(i,CT{i});assist_torque_norm{i,:} = interp1(time_c{i,:},assist_torque_c{i,:},range{i,:});
    end 
end 

%% Figure 3
subplot(2,1,1)
range1 = (0:0.01:0.99)';
colors = {[0.1641 0.6133 0.9531], [0.0938 0.4805 0.8008], [0.0664 0.4023 0.6914]};
styles = {':','-','-.'};


hold on
for i = 7:9
    idx = i - 6;
    plot(range1, -assist_torque_norm{i,:}, ...
        'Color', colors{idx}, ...
        'LineStyle', styles{idx}, ...
        'LineWidth', 2);
end
plot(range1, zeros(size(range1)), 'k--', 'LineWidth', 1);

xticks(0:0.2:1)
xticklabels({'0','20','40','60','80','100'})
xlabel('Percent Gait Cycle (%)')
ylabel('Torque (N.m)')
set(gca,'FontSize',10,'TickLength',[0.02 0.02],'LineWidth',1)
box off
legend({'Original model','95th impaired male model','5th impaired female model'},'FontSize',10)
yl = ylim;   
hPatch = patch([0.43 0.64 0.64 0.43], [yl(1) yl(1) yl(2) yl(2)], [0.8 0.8 0.8], 'EdgeColor','none', 'FaceAlpha',0.3); 
uistack(findobj(gca,'Type','patch'),'bottom');
hPatch.Annotation.LegendInformation.IconDisplayStyle = 'off';

exportgraphics(gcf,'Fig 3.tif','Resolution',600);
%% Figure 4
data2 = {ankle_r_norm,knee_r_norm,hip_r_norm, ankle_r_moment_norm, knee_r_moment_norm,hip_r_moment_norm, grf_r_y_norm, grf_r_x_norm};
M = [75, 96, 49]; % model mass

for k = 1:3
    for i = 1:length(data2)
        Y = data2{i};
        signFactor = 1;
            if ismember(i, [2 4])
                  signFactor = -1;
            end 
            if ismember(i, [4 5 6])
                de = M(k);
            elseif ismember(i, [7 8])
                de = M(k) * 9.81;
            else
                de = 1;
            end
         Y{k,:} = (Y{k,:}/de)*signFactor;
         Y{k+3,:} = (Y{k+3,:}/de)*signFactor;
         Y{k+6,:} = (Y{k+6,:}/de)*signFactor;
         data2{i} = Y;
    end   
end 


for i = 1:length(data2)
        mean_unimpaired_data2(i,:) = mean([cell2mat(data2{i}(1,:)); cell2mat(data2{i}(2,:)); cell2mat(data2{i}(3,:))]);
        std_unimpaired_data2(i,:)  = std([cell2mat(data2{i}(1,:)); cell2mat(data2{i}(2,:)); cell2mat(data2{i}(3,:))]);

        mean_pf_data2(i,:) = mean([cell2mat(data2{i}(4,:)); cell2mat(data2{i}(5,:)); cell2mat(data2{i}(6,:))]);
        std_pf_data2(i,:)  = std([cell2mat(data2{i}(4,:)); cell2mat(data2{i}(5,:)); cell2mat(data2{i}(6,:))]);

        mean_assist_data2(i,:) = mean([cell2mat(data2{i}(7,:)); cell2mat(data2{i}(8,:)); cell2mat(data2{i}(9,:))]);
        std_assist_data2(i,:)  = std([cell2mat(data2{i}(7,:)); cell2mat(data2{i}(8,:)); cell2mat(data2{i}(9,:))]);
end 

figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.7]); tiledlayout(3,3,'Padding','compact','TileSpacing','compact');
range1 = (0:0.01:0.99)';
for i = 1:length(data2)
    tileMap = [1 4 7 2 5 8 3 6 9];
    nexttile(tileMap(i));hold on;
    mu_PF = mean_pf_data2(i,:); sd_PF = std_pf_data2(i,:);
    upper_PF = mu_PF + sd_PF;lower_PF = mu_PF - sd_PF;
    x2 = [range1' fliplr(range1')];inPF = [lower_PF fliplr(upper_PF)];
    h_imp = plot(range1', mu_PF, 'Color',[0 0 1], 'LineWidth', 2);
    fill(x2, inPF, [0 0 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    mu_AS = mean_assist_data2(i,:);sd_AS = std_assist_data2(i,:);
    upper_AS = mu_AS + sd_AS;lower_AS = mu_AS - sd_AS;
    inAS = [lower_AS fliplr(upper_AS)];
    h_unimp = plot(range1', mu_AS, 'Color',[1 0 0], 'LineWidth', 2);
    fill(x2, inAS, [1 0 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    xticks(0:0.2:1);
           yLabels = {'Ankle dorsiflexion','Knee flexion','Hip flexion'};
        if mod(tileMap(i)-1,3) == 0
             ylabel(yLabels{(tileMap(i)-1)/3 + 1}, 'FontSize', 10);
        end
        tTitles = {'Angles(deg)','Moments(N.m/kg)','GRF(BW)'};
        t = tileMap(i);
        

        if ismember(t, [1 2 3])
           title(tTitles{t}, 'FontSize', 15, 'FontWeight', 'bold');
        end
        if ismember(tileMap(i), [7 8])
            xticks(0:0.2:1);
            xticklabels({'0','20','40','60','80','100'});
            xlabel('Percent Gait Cycle (%)');
        else
            xticklabels({});
        end
         if mod(tileMap(i)-1,3) == 2   
            if tileMap(i) == 3 || tileMap(i) == 6
                if tileMap(i) == 3
                    ylabel('Vertical', 'FontSize', 10);
                elseif tileMap(i) == 6
                    ylabel('Horizontal', 'FontSize', 10);
                end
            end
        end
    ax = gca;
    box off
    set(ax,'FontSize',15,'TickLength',[0.02 0.02],'LineWidth',1);
end
lgd = legend([h_imp h_unimp], {'Impaired gait','Unimpaired gait'});
exportgraphics(gcf,'Fig 4.tif','Resolution',600);
%% Figure 5

M = [75, 96, 49]; % model mass 
muscles = {'HAMS','BFSH','GMAX','ILPSO','RF','VAS','GAS','SOL','TA'};
data = {HAMS_norm, BFSH_norm, GMAX_norm, ILPSO_norm, RF_norm, VAS_norm, GAS_norm, SOL_norm, TA_norm};

for i = 1:length(muscles)
        mean_unimpaired_muscles(i,:) = mean([cell2mat(data{i}(1,:))/M(1); cell2mat(data{i}(2,:))/M(2); cell2mat(data{i}(3,:))/M(3)]);
        std_unimpaired_muscles(i,:)  = std([cell2mat(data{i}(1,:))/M(1); cell2mat(data{i}(2,:))/M(2); cell2mat(data{i}(3,:))/M(3)]);

        mean_pf_muscles(i,:) = mean([cell2mat(data{i}(4,:))/M(1); cell2mat(data{i}(5,:))/M(2); cell2mat(data{i}(6,:))/M(3)]);
        std_pf_muscles(i,:)  = std([cell2mat(data{i}(4,:))/M(1); cell2mat(data{i}(5,:))/M(2); cell2mat(data{i}(6,:))/M(3)]);

        mean_assist_muscles(i,:) = mean([cell2mat(data{i}(7,:))/M(1); cell2mat(data{i}(8,:))/M(2); cell2mat(data{i}(9,:))/M(3)]);
        std_assist_muscles(i,:)  = std([cell2mat(data{i}(7,:))/M(1); cell2mat(data{i}(8,:))/M(2); cell2mat(data{i}(9,:))/M(3)]);
end 

figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.7]); tiledlayout(3,3,'Padding','compact','TileSpacing','compact');

for i = 1:9
    nexttile;hold on;
    mu_PF = mean_pf_muscles(i,:); sd_PF = std_pf_muscles(i,:);
    upper_PF = mu_PF + sd_PF;lower_PF = mu_PF - sd_PF;
    x2 = [range1' fliplr(range1')];inPF = [lower_PF fliplr(upper_PF)];
    h_imp = plot(range1', mu_PF, 'Color',[0 0 1], 'LineWidth', 2);
    fill(x2, inPF, [0 0 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    mu_AS = mean_assist_muscles(i,:);sd_AS = std_assist_muscles(i,:);
    upper_AS = mu_AS + sd_AS;lower_AS = mu_AS - sd_AS;
    inAS = [lower_AS fliplr(upper_AS)];
    h_unimp = plot(range1', mu_AS, 'Color',[1 0 0], 'LineWidth', 2);
    fill(x2, inAS, [1 0 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    xticks(0:0.2:1);
    if i <= 6
        xticklabels({});
    else
        xticklabels({'0','20','40','60','80','100'});xlabel('Percent Gait Cycle (%)');
    end
    if mod(i-1,3) == 0
       ylabel('N/kg', 'FontSize', 10);
    end
    
    title(muscles{i}, 'FontSize', 15, 'FontWeight', 'bold');
    ax = gca;
    box off
    set(ax,'FontSize',15,'TickLength',[0.02 0.02],'LineWidth',1);
end
lgd = legend([h_imp h_unimp], {'Impaired gait','Unimpaired gait'});
exportgraphics(gcf,'Fig 5.tif','Resolution',600);

%% Figure 6

figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.7]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

%%% ============= Muscle cost of transport ============
nexttile([1 3]);
M = [75, 96, 49]; % model mass 
muscles = {'HAMS','BFSH','GMAX','ILPSO','RF','VAS','GAS','SOL','TA'};
data = {HAMS_norm, BFSH_norm, GMAX_norm, ILPSO_norm, RF_norm, VAS_norm, GAS_norm, SOL_norm, TA_norm};

for i = 1:length(muscles)
        mean_unimpaired_muscles(i,:) = mean([cell2mat(data{i}(1,:))/M(1); cell2mat(data{i}(2,:))/M(2); cell2mat(data{i}(3,:))/M(3)]);
        std_unimpaired_muscles(i,:)  = std([cell2mat(data{i}(1,:))/M(1); cell2mat(data{i}(2,:))/M(2); cell2mat(data{i}(3,:))/M(3)]);

        mean_pf_muscles(i,:) = mean([cell2mat(data{i}(4,:))/M(1); cell2mat(data{i}(5,:))/M(2); cell2mat(data{i}(6,:))/M(3)]);
        std_pf_muscles(i,:)  = std([cell2mat(data{i}(4,:))/M(1); cell2mat(data{i}(5,:))/M(2); cell2mat(data{i}(6,:))/M(3)]);

        mean_assist_muscles(i,:) = mean([cell2mat(data{i}(7,:))/M(1); cell2mat(data{i}(8,:))/M(2); cell2mat(data{i}(9,:))/M(3)]);
        std_assist_muscles(i,:)  = std([cell2mat(data{i}(7,:))/M(1); cell2mat(data{i}(8,:))/M(2); cell2mat(data{i}(9,:))/M(3)]);
end 


X = categorical(muscles);
X = reordercats(X,muscles);

for i = 1:length(muscles)
    Y4(i,:)   = [mean(mean_pf_muscles(i,:)) mean(mean_assist_muscles(i,:))]/10;
    std4(i,:) = [mean(std_pf_muscles(i,:))  mean(std_assist_muscles(i,:))]/10;
end

b = bar(X,Y4); hold on
b(1).FaceColor = [0 0 1];
b(2).FaceColor = [1 0 0];

ng = size(Y4,1);
nb = size(Y4,2);
gw = min(0.8, nb/(nb+1.5));

for i = 1:nb
    x = (1:ng) - gw/2 + (2*i-1)*gw/(2*nb);
    errorbar(x,Y4(:,i),std4(:,i),'.k','LineWidth',1.5);
end

legend({'Unassisted gait','Assisted gait'})
ylabel('Muscle Cost of Transport (J/kg/m)')
set(gca,'FontSize',15,'TickLength',[0.02 0.02],'LineWidth',1)
box off

%%% ============= Cost of Transport ============
M = [75, 96, 49]; % model mass 
muscles = {'HAMS','BFSH','GMAX','ILPSO','RF','VAS','GAS','SOL','TA'};
data = {HAMS_norm, BFSH_norm, GMAX_norm, ILPSO_norm, RF_norm, VAS_norm, GAS_norm, SOL_norm, TA_norm};

for m=4:NK
    tmp = 0;
    for k=1:length(muscles)
        tmp = tmp + cell2mat(data{k}(m,:));
    end 
    COT(m) = mean(tmp);   
end 

COT(4)=COT(4)/(M(1)*Gait_speed{4});
COT(5)=COT(5)/(M(2)*Gait_speed{5});
COT(6)=COT(6)/(M(3)*Gait_speed{6});
COT(7)=COT(7)/(M(1)*Gait_speed{7});
COT(8)=COT(8)/(M(2)*Gait_speed{8});
COT(9)=COT(9)/(M(3)*Gait_speed{9});

mean_COT_pf = mean(COT(4:6));
std_COT_pf = std(COT(4:6));

mean_COT_assist = mean(COT(7:9));
std_COT_assist  = std(COT(7:9));

X = categorical({'Unassist','Assist'});
X = reordercats(X,{'Unassist','Assist'});

Y = [mean_COT_pf mean_COT_assist]/10;
S = [std_COT_pf std_COT_assist]/10;

nexttile
X_num = double(X);
bar(X_num(1)+0.2, Y(1), 'barWidth',0.4,'FaceColor',[0 0 1],'EdgeColor','flat','LineWidth',1);hold on  
bar(X_num(2)-0.2, Y(2), 'barWidth',0.4,'FaceColor',[1 0 0],'EdgeColor','flat','LineWidth',1);  

errorbar(X_num(1)+0.2, Y(1), S(1), 'k', 'LineWidth', 1.5, 'CapSize', 10);  
errorbar(X_num(2)-0.2, Y(2), S(2), 'k', 'LineWidth', 1.5, 'CapSize', 10); 

xticks([1, 2]); 
xticklabels({'Unassisted', 'Assisted'});  

% Perform two-sample t-test
[~, p_value] = ttest2(COT(4:6), COT(7:9));
max_y = max(Y + S);y_bracket = max_y + 0.5;
plot([1.2, 1.2, 1.8, 1.8], [y_bracket, y_bracket+0.2,y_bracket+0.2, y_bracket], '-k', 'LineWidth', 1.5); 

if p_value < 0.001
    significance = '***';  
elseif p_value < 0.01
    significance = '**';  
elseif p_value < 0.05
    significance = '*';    
else
    significance = 'n.s.'; 
end
text(1.5, y_bracket + 0.3, significance, 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
ylim([0, y_bracket + 1]);
ylabel('Cost of Transport(J/kg/m)')
set(gca,'fontsize',15,'TickLength',[0.02 0.02],'LineWidth', 1);box off

%%% ============= Gait Speed ============
nexttile;hold on
colors  = {[0 0.6 0], [0.7 0.3 0.5], [1 0.5 0.1]};  
markers = {'o','^','s'};
for k = 4:6
    idx = k-3;   
    h(idx) = plot([1 2], [Gait_speed{k}, Gait_speed{k+3}], 'Color', colors{idx},'LineWidth', 2);
    plot(1, Gait_speed{k}, markers{idx}, 'MarkerSize',10,'MarkerFaceColor',colors{idx},'MarkerEdgeColor',colors{idx});
    plot(2, Gait_speed{k+3},markers{idx}, 'MarkerSize',10,'MarkerFaceColor',colors{idx},'MarkerEdgeColor',colors{idx});
end
axis([0.5 2.5 0.7 1])
xticks([1 2])
xticklabels({'Unassisted','Assisted'})
ylabel('Gait Speed (m/s)')
box off;set(gca,'FontSize',15,'TickLength',[0.02 0.02],'LineWidth',1)

%%% ============= Speed Length ============
nexttile
hold on
for k = 4:6
    idx = k-3;
    h(idx) = plot([1 2], [step_length{k}, step_length{k+3}], 'Color', colors{idx}, 'LineWidth', 2);
    plot(1, step_length{k}, markers{idx}, 'MarkerSize',10,'MarkerFaceColor',colors{idx},'MarkerEdgeColor',colors{idx});
    plot(2, step_length{k+3}, markers{idx}, 'MarkerSize',10,'MarkerFaceColor',colors{idx},'MarkerEdgeColor',colors{idx});
end
axis([0.5 2.5 0.95 1.25])
xticks([1 2])
xticklabels({'Unassisted','Assisted'})
ylabel('Step Length (m)')
set(gca,'FontSize',15,'TickLength',[0.02 0.02],'LineWidth',1);box off
legend(h, {'Original model', '95th-tile male model','5th-tile female model'},'fontsize',15);
exportgraphics(gcf,'Fig 6.tif','Resolution',600);

%% Supplementary Figures

%%% ============= S1.Fig ============
colors = {[0 0 0], [0 0 1], [1 0 0]};styles = {':', '-', '-'};
data2 = {ankle_r_norm,knee_r_norm,hip_r_norm, ankle_r_moment_norm, knee_r_moment_norm,hip_r_moment_norm, grf_r_y_norm, grf_r_x_norm};
for k = 1:3
    figure('Units', 'normalized','Position', [0.1 0.1 0.8 0.7]);tiledlayout(3,3, 'Padding','compact', 'TileSpacing','compact');
    for i = 1:length(data2)
        tileMap = [1 4 7 2 5 8 3 6 9];
        nexttile(tileMap(i));
        hold on;
        Y = data2{i};
            if ismember(i, [4 5 6])
                de = M(k);
            elseif ismember(i, [7 8])
                de = M(k) * 9.81;
            else
                de = 1;
            end
        plot(range1, Y{k,:}/de, 'Color', colors{1}, 'LineWidth', 2, 'LineStyle', styles{1});
        plot(range1, Y{k+3,:}/de, 'Color', colors{2}, 'LineWidth', 2, 'LineStyle', styles{2});
        plot(range1, Y{k+6,:}/de,'Color', colors{3}, 'LineWidth', 2, 'LineStyle', styles{3});
        xticks(0:0.2:1);
        yLabels = {'Ankle dorsiflexion','Knee flexion','Hip flexion'};
        if mod(tileMap(i)-1,3) == 0
             ylabel(yLabels{(tileMap(i)-1)/3 + 1}, 'FontSize', 10);
        end
         if mod(tileMap(i)-1,3) == 2   
            if tileMap(i) == 3 || tileMap(i) == 6
                if tileMap(i) == 3
                    ylabel('Vertical', 'FontSize', 10);
                elseif tileMap(i) == 6
                    ylabel('Horizontal', 'FontSize', 10);
                end
            end
        end

        tTitles = {'Angles(deg)','Moments(N.m/kg)','GRF(BW)'};
        t = tileMap(i);

        if ismember(t, [1 2 3])
           title(tTitles{t}, 'FontSize', 15, 'FontWeight', 'bold');
        end
        if ismember(tileMap(i), [7 8])
            xticks(0:0.2:1);
            xticklabels({'0','20','40','60','80','100'});
            xlabel('Percent Gait Cycle (%)');
        else
            xticklabels({});
        end

        ax = gca;
        box off
        set(ax,'FontSize',15,'TickLength',[0.02 0.02], 'LineWidth',1);
    end
    switch k
        case 1
            legendLabels = {'Unimpaired Gait', 'PF Gait', 'Assisted Gait'};
        case 2
            legendLabels = {'Unimpaired Gait-Model 95%','PF Gait-Model 95%', 'Assisted Gait-Model 95%'};
        case 3
            legendLabels = {'Unimpaired Gait-Model 5%', 'PF Gait-Model 5%', 'Assisted Gait-Model 5%'};
    end
    lgd = legend(legendLabels);
    switch k
        case 1
            exportgraphics(gcf,'S1A.Fig.tif','Resolution',600);
        case 2
            exportgraphics(gcf,'S1B.Fig.tif','Resolution',600);
        case 3
            exportgraphics(gcf,'S1C.Fig.tif','Resolution',600);
    end
end



%%% ============= S2.Fig ============
M = [75, 96, 49]; % model mass 
muscles = {'HAMS','BFSH','GMAX','ILPSO','RF','VAS','GAS','SOL','TA'};
data = {HAMS_norm, BFSH_norm, GMAX_norm, ILPSO_norm, RF_norm, VAS_norm, GAS_norm, SOL_norm, TA_norm};
range1 = (0:0.01:0.99)';colors = {[0 0 0], [0 0 1], [1 0 0]};styles = {':', '-', '-'};

for k = 1:3
    figure('Units', 'normalized','Position', [0.1 0.1 0.8 0.7]);tiledlayout(3,3, 'Padding','compact', 'TileSpacing','compact');
    for i = 1:9
        nexttile(i);
        hold on;
        Y = data{i};
        plot(range1, Y{k,:}/M(k), 'Color', colors{1}, 'LineWidth', 2, 'LineStyle', styles{1});
        plot(range1, Y{k+3,:}/M(k), 'Color', colors{2}, 'LineWidth', 2, 'LineStyle', styles{2});
        plot(range1, Y{k+6,:}/M(k),'Color', colors{3}, 'LineWidth', 2, 'LineStyle', styles{3});
        title(muscles{i}, 'FontSize', 15, 'FontWeight', 'bold');
        xticks(0:0.2:1);
        if i <= 6
            xticklabels({});
        else
            xticklabels({'0','20','40','60','80','100'});
            xlabel('Percent Gait Cycle (%)');
        end
        if mod(i-1,3) == 0
            ylabel('N/kg', 'FontSize', 10);
        end
        ax = gca;
        box off
        set(ax,'FontSize',15,'TickLength',[0.02 0.02], 'LineWidth',1);
    end
    switch k
        case 1
            legendLabels = {'Unimpaired Gait', 'PF Gait', 'Assisted Gait'};
        case 2
            legendLabels = {'Unimpaired Gait-Model 95%','PF Gait-Model 95%', 'Assisted Gait-Model 95%'};
        case 3
            legendLabels = {'Unimpaired Gait-Model 5%', 'PF Gait-Model 5%', 'Assisted Gait-Model 5%'};
    end
    lgd = legend(legendLabels);

    switch k
        case 1
            exportgraphics(gcf,'S2A.Fig.tif','Resolution',600);
        case 2
            exportgraphics(gcf,'S2B.Fig.tif','Resolution',600);
        case 3
            exportgraphics(gcf,'S2C.Fig.tif','Resolution',600);
    end
end

