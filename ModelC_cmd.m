%% --------------- Matlab codes of "Oscillation" ------------------------------------------------------------------ %%
close all
clear all
clc
tspan = [0:100:1000000]; 
%

%% --------------- PARAMETERS ------------------------------------------------------------------------------------- %%
param.k0 = 81;
param.k1 = 20;
param.k2a = 22;
param.k2b = 1;
param.k2c = 60;
param.Km2c = 3;
param.k3 = 0.012;
param.k4f = 1.4e-3;
param.k4b = 1e-3;
param.k4c = 6e-3;
param.k4d = 2.91e-4;
param.k5 = 110;
param.k6 = 100;
param.k7 = 3000;
param.k8 = 0.0058;
param.k9 = 3.85e-4;
param.k10 = 0.01185;
param.k11 = 1.25e-6;
param.k12f = 7.43e-5;
param.k12b = 1.25e-6;
param.k13 = 7.43e-5;
param.k14 = 7.43e-5;
param.k15 = 7.43e-5;
param.k16 = 3.85e-4;
param.HSP90 = 1;
param.V_ratio = 0.3;
param.PRXtot = 100;
param.TRXtot = 30;
param.TRtot = 5.52;
default_param = param;
%

%% --------------- INITIAL CONDITION ------------------------------------------------------------------------------ %%
init.H2O2_cyto = 0;
init.H2O2_mito = 0.2;
init.PRXSH = 100;
init.PRXSOH = 0;
init.PRXSO2H = 0;
init.SRX_mito = 0.1;
init.SRXOH_cyto = 0.1;
init.SRX_cyto = 13;
init.SRXOH_HSP90_cyto = 0;
init.PRXSO2H_SRX = 0;
init.TRXSS = 0;
default_init = init;
%

%% --------------- Run to stable oscillation ---------------------------------------------------------------------- %%
tspan0 = [0:100:3600*24*100]; 
y0 = cell2mat(struct2cell(init));
[t,y] = ode15s('ModelC_ode', tspan0, y0, [], param); 
%

%% --------------- Generate final_fig. 10A Time course for oscillation -------------------------------------------- %%
% run to stable oscillation
tspan0 = [0:100:3600*24*5.5]; 
y0 = cell2mat(struct2cell(init));
[t2,y2] = ode15s('ModelC_ode', tspan0, y0, [], param);

% start stable oscillation
tspan3 = [0:10:3600*24*5]; 
y0 = y2(end,:); 
[t3,y3] = ode15s('ModelC_ode', tspan3, y0, [], param);

pb_ratio = 1.8;
font_size = 22;
x_lim = 24*5;
line_width = 3;


% % Non_PRXSO2H = PRX + PRXSOH + PRXSS
% figure(11111)
% plot(t3/3600, param.PRXtot - (y3(:,7) + y3(:,10)),'LineWidth', line_width)
% set(gca, 'xtick', 0:24:240)
% xlim([0,x_lim])
% % ylim([0,8e-3])
% xlabel ('Time (h)')
% ylabel ('Non PRXSO2H (uM)')
% set(gca,'fontsize',font_size);
% set(get(gca,'XLabel'),'FontSize',font_size);
% set(get(gca,'YLabel'),'FontSize',font_size);
% pbaspect([pb_ratio 1 1])
% hold on


% H2O2_mito
figure(1101)
plot(t3/3600, y3(:,2),'LineWidth', line_width)
set(gca, 'xtick', 0:24:240)
xlim([0,x_lim])
% ylim([0,8e-3])
xlabel ('Time (h)')
ylabel ('Mito H2O2 (uM)')
set(gca,'fontsize',font_size);
set(get(gca,'XLabel'),'FontSize',font_size);
set(get(gca,'YLabel'),'FontSize',font_size);
pbaspect([pb_ratio 1 1])
hold on

% H2O2_cyto
figure(1102)
plot(t3/3600, y3(:,1),'LineWidth', line_width)
set(gca, 'xtick', 0:24:240)
xlim([0,x_lim])
% ylim([0,8e-3])
xlabel ('Time (h)')
ylabel ('Cyto H2O2 (uM)')
set(gca,'fontsize',font_size);
set(get(gca,'XLabel'),'FontSize',font_size);
set(get(gca,'YLabel'),'FontSize',font_size);
pbaspect([pb_ratio 1 1])
ax=gca; ax.YAxis.Exponent = -2;
hold on

% PRXSOH
figure(1103)
plot(t3/3600, y3(:,9),'LineWidth', line_width)
set(gca, 'xtick', 0:24:240)
xlim([0,x_lim])
% ylim([0,8e-3])
xlabel ('Time (h)')
ylabel ('PRXSOH (uM)')
set(gca,'fontsize',font_size);
set(get(gca,'XLabel'),'FontSize',font_size);
set(get(gca,'YLabel'),'FontSize',font_size);
pbaspect([pb_ratio 1 1])
hold on

% PRXSS = PRXtot - PRXSH - PRXSOH - PRXSO2H - PRXSO2H_SRX
figure(1104)
plot(t3/3600, (param.PRXtot - y3(:,8) - y3(:,9) - y3(:,10) - y3(:,7)),'LineWidth', line_width)
set(gca, 'xtick', 0:24:240)
xlim([0,x_lim])
% ylim([0,8e-3])
xlabel ('Time (h)')
ylabel ('PRXSS (uM)')
set(gca,'fontsize',font_size);
set(get(gca,'XLabel'),'FontSize',font_size);
set(get(gca,'YLabel'),'FontSize',font_size);
pbaspect([pb_ratio 1 1])
hold on

% PRXSO2H_tot = PRXSO2H + PS;
figure(1105)
plot(t3/3600, y3(:,7) + y3(:,10),'LineWidth', line_width)
set(gca, 'xtick', 0:24:240)
xlim([0,x_lim])
% ylim([0,8e-3])
xlabel ('Time (h)')
ylabel ('PRXSO2H tot (uM)')
set(gca,'fontsize',font_size);
set(get(gca,'XLabel'),'FontSize',font_size);
set(get(gca,'YLabel'),'FontSize',font_size);
pbaspect([pb_ratio 1 1])
hold on

% TRXSH = TRXtot - TRXSS;
figure(1106)
plot(t3/3600, (param.TRXtot - y3(:,11)),'LineWidth', line_width)
set(gca, 'xtick', 0:24:240)
xlim([0,x_lim])
% ylim([0,8e-3])
xlabel ('Time (h)')
ylabel ('TRXSH (uM)')
set(gca,'fontsize',font_size);
set(get(gca,'XLabel'),'FontSize',font_size);
set(get(gca,'YLabel'),'FontSize',font_size);
pbaspect([pb_ratio 1 1])
hold on

% TRXSS;
figure(1107)
plot(t3/3600, y3(:,11),'LineWidth', line_width)
set(gca, 'xtick', 0:24:240)
xlim([0,x_lim])
% ylim([0,8e-3])
xlabel ('Time (h)')
ylabel ('TRXSS (uM)')
set(gca,'fontsize',font_size);
set(get(gca,'XLabel'),'FontSize',font_size);
set(get(gca,'YLabel'),'FontSize',font_size);
pbaspect([pb_ratio 1 1])
hold on

% SRX_mito_tot = SRX_mito + PS
figure(1108)
plot(t3/3600, y3(:,6) + y3(:,7),'LineWidth', line_width)
set(gca, 'xtick', 0:24:240)
xlim([0,x_lim])
% ylim([0,8e-3])
xlabel ('Time (h)')
ylabel ('Mito SRX tot (uM)')
set(gca,'fontsize',font_size);
set(get(gca,'XLabel'),'FontSize',font_size);
set(get(gca,'YLabel'),'FontSize',font_size);
pbaspect([pb_ratio 1 1])
hold on

% SRX_cyto_tot = SRX_cyto + SRXOH_cyto + SRXOH_HSP90_cyto;
figure(1109)
plot(t3/3600, y3(:,3) + y3(:,4) + y3(:,5),'LineWidth', line_width)
set(gca, 'xtick', 0:24:240)
xlim([0,x_lim])
% ylim([0,8e-3])
xlabel ('Time (h)')
ylabel ('Cyto SRX tot (uM)')
set(gca,'fontsize',font_size);
set(get(gca,'XLabel'),'FontSize',font_size);
set(get(gca,'YLabel'),'FontSize',font_size);
pbaspect([pb_ratio 1 1])
hold on


%% --------------- Generate final_fig. 10B and 10C ---------------------------------------------------------------- %%
% Caculate normalized distance in each variable

H2O2_mito = y3(:,2);
H2O2_mito_max = max(H2O2_mito);
dist_H2O2_mito = (H2O2_mito(2:end) - H2O2_mito(1:end-1))/H2O2_mito_max;

PRXSO2H_tot = y3(:,7) + y3(:,10);
PRXSO2H_tot_max = max(PRXSO2H_tot);
dist_PRXSO2H_tot = (PRXSO2H_tot(2:end) - PRXSO2H_tot(1:end-1))/PRXSO2H_tot_max;

SRX_mito_tot = (y3(:,6) + y3(:,7));
SRX_mito_tot_max = max(SRX_mito_tot);
dist_SRX_mito_tot = (SRX_mito_tot(2:end) - SRX_mito_tot(1:end-1))/SRX_mito_tot_max;


    % k0 = 68 (default value)
    figure (1107)
    % Overlay XPP-aut H2O2 vs. SRX_tot bifurcation
    [num, txt] = xlsread('ModelC_SRX_saddle_node_bifurcation_k0_default.xlsx', 1);
    column_num = length(txt);
    
    % Obtain branch index: 
    branch_1_index = find(num(:, column_num) == 1); 
    branch_2_index = find(num(:, column_num) == 2);
    branch_3_index = find(num(:, column_num) == 3);
    branch_4_index = find(num(:, column_num) == 4);

    % Obtain branch data: 
    branch_1 = num(branch_1_index, [1, 2]);
    branch_2 = num(branch_2_index, [1, 2]);
    branch_3 = num(branch_3_index, [1, 2]);
    branch_4 = num(branch_4_index, [1, 2]);

    plot(branch_4(:, 1), branch_4(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    hold on
    plot(branch_1(:, 1), branch_1(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    plot(branch_2(:, 1), branch_2(:, 2), ':', 'Color', [51,102,255]/255, 'LineWidth', 3);
    plot(branch_3(:, 1), branch_3(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    pbaspect([1.1 1 1])
    xlabel(txt(1))
    ylabel(txt(2))
%     xlim([0, 5.5])
%     ylim([0, 0.3])
    box on


    % k0 = 57.917 (HB1 value)
    figure (1107)
    % Overlay XPP-aut H2O2 vs. SRX_tot bifurcation
    [num, txt] = xlsread('ModelC_SRX_saddle_node_bifurcation_k0_HB1.xlsx', 1);
    column_num = length(txt);
    
    % Obtain branch index: 
    branch_1_index = find(num(:, column_num) == 1); 
    branch_2_index = find(num(:, column_num) == 2);
    branch_3_index = find(num(:, column_num) == 3);
    branch_4_index = find(num(:, column_num) == 4);
    
    % Obtain branch data: 
    branch_1 = num(branch_1_index, [1, 2]);
    branch_2 = num(branch_2_index, [1, 2]);
    branch_3 = num(branch_3_index, [1, 2]);
    branch_4 = num(branch_4_index, [1, 2]);

    plot(branch_4(:, 1), branch_4(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    hold on
    plot(branch_1(:, 1), branch_1(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    plot(branch_2(:, 1), branch_2(:, 2), ':', 'Color', [51,102,255]/255, 'LineWidth', 3);
    plot(branch_3(:, 1), branch_3(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    pbaspect([1.1 1 1])
    xlabel(txt(1))
    ylabel(txt(2))
%     xlim([0, 5.5])
%     ylim([0, 0.3])
    box on


    % k0 = 84.281 (HB2 value)
    figure (1107)
    % Overlay XPP-aut H2O2 vs. SRX_tot bifurcation
    [num, txt] = xlsread('ModelC_SRX_saddle_node_bifurcation_k0_HB2.xlsx', 1);
    column_num = length(txt);
    
    % Obtain branch index: 
    branch_1_index = find(num(:, column_num) == 1); 
    branch_2_index = find(num(:, column_num) == 2);
    branch_3_index = find(num(:, column_num) == 3);
    branch_4_index = find(num(:, column_num) == 4);
    
    % Obtain branch data: 
    branch_1 = num(branch_1_index, [1, 2]);
    branch_2 = num(branch_2_index, [1, 2]);
    branch_3 = num(branch_3_index, [1, 2]);
    branch_4 = num(branch_4_index, [1, 2]);
   
    plot(branch_4(:, 1), branch_4(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    hold on
    plot(branch_1(:, 1), branch_1(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    plot(branch_2(:, 1), branch_2(:, 2), ':', 'Color', [51,102,255]/255, 'LineWidth', 3);
    plot(branch_3(:, 1), branch_3(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    pbaspect([1.1 1 1])
    xlabel(txt(1))
    ylabel(txt(2))
%     xlim([0, 5.5])
%     ylim([0, 0.3])
    box on




% Overlay Mito SRX tot null curve by changing Mito H2O2;  

% ****** needs to set H2O2_mito ODE to zero in ModelC_ode.m file  ******
% ****** needs to set H2O2_mito ODE to zero in ModelC_ode.m file  ******
% ****** needs to set H2O2_mito ODE to zero in ModelC_ode.m file  ******

init = default_init;
param = default_param;

Mito_SRX_tot_vector = [];
Mito_H2O2_vector = [];

init.H2O2_mito = 0.001;
tspan1 = [0:100:3600*24*100]; 

while init.H2O2_mito <= 10

        Mito_H2O2_vector = [Mito_H2O2_vector, init.H2O2_mito];
        y0 = cell2mat(struct2cell(init));
        [t1,y1]=ode15s('ModelC_ode', tspan1, y0, [], param); % needs to set H2O2_mito ODE to zero in ModelC_ode.m file;

        Mito_SRX_tot_vector  =  [Mito_SRX_tot_vector, y1(end,6) + y1(end,7)]; 
        init.H2O2_mito =  init.H2O2_mito * 1.05;
 
end

figure(1107)
plot(Mito_SRX_tot_vector , Mito_H2O2_vector, 'Color', [51,102,255]/255, 'LineWidth', 3);
% plot(Mito_H2O2_vector, Mito_SRX_tot_vector, 'Color', [51,102,255]/255, 'LineWidth', 3);
% xlim([0, 1.5])
% ylim([0, 0.3])

% Overlay trajectory with velocity
    dist = sqrt(dist_H2O2_mito.^2 + dist_SRX_mito_tot.^2);
    % dist = sqrt(dist_H2O2_mito.^2);
    % dist = sqrt(dist_SRX_mito_tot.^2);
    
    dist_max = max(dist);
    norm_dist = dist/dist_max; % nomalized distance
    norm_dist_log10 = log10(norm_dist);
    
    figure(1107)
    scatter(SRX_mito_tot(1:end-1), H2O2_mito(1:end-1), 50, norm_dist, 'filled') 
%   scatter(SRX_mito_tot(1:end-1), H2O2_mito(1:end-1), 50, norm_dist_log10, 'filled') 
    xlabel ('Mito SRX tot (uM)')
    ylabel ('Mito H2O2 (uM)')
    colorbar
    box on
    hold on
    xlim([0.3 1.5])
    ylim([0.04, 0.25])
    pbaspect([1 1 0.9])

% Figure 10C, % 3D-plot
dist = sqrt(dist_H2O2_mito.^2 + dist_PRXSO2H_tot.^2 + dist_SRX_mito_tot.^2);
dist_max = max(dist);
norm_dist = dist/dist_max; % nomalized distance
norm_dist_log10 = log10(norm_dist);

figure (1108)
scatter3(SRX_mito_tot(1:end-1), H2O2_mito(1:end-1), PRXSO2H_tot(1:end-1), 50, norm_dist, 'filled') 
% scatter3(SRX_mito_tot(1:end-1), H2O2_mito(1:end-1), PRXSO2H_tot(1:end-1), 50, norm_dist_log10, 'filled')  
pbaspect([1 1 0.9])
xlabel ('Mito SRX tot (uM)')
ylabel ('Mito H2O2 (uM)')
zlabel ('PRXSO2H tot (uM)')
colorbar
box on

% ！！！ reset H2O2_mito ODE in ModelC_ode.m file;
% ！！！ reset H2O2_mito ODE in ModelC_ode.m file;
% ！！！ reset H2O2_mito ODE in ModelC_ode.m file;

%





%% --------------- Generate final_fig. 11A for Supercritical_hopf_bifurcation_Mito_H2O2_vs_k0.xlsx ---------------- %%
[num, txt] = xlsread('Supercritical_hopf_bifurcation_Mito_H2O2_vs_k0.xlsx');

% Obtain branch index: 
branch_1_index = find(num(:, 4) == 1); 
branch_2_index = find(num(:, 4) == 2);
branch_3_index = find(num(:, 4) == 3);
branch_4_index = find(num(:, 4) == 4);
    
% Obtain branch data: 
branch_1 = num(branch_1_index, [1:3]);
branch_2 = num(branch_2_index, [1:3]);
branch_3 = num(branch_3_index, [1:3]);
branch_4 = num(branch_4_index, [1:3]);

font_size = 24;
% Plot bifurcation: 

figure(1101)
% 1: Stable fixed point, to left hopf bifurcation point;
plot(branch_1(:, 1), branch_1(:, 2), 'Color', [255,50,50]/255, 'LineWidth', 3);
hold on

% 1: Stable fixed point, to hopf bifurcation point;
plot(branch_4(:, 1), branch_4(:, 2), 'Color', [255,50,50]/255, 'LineWidth', 3);

% 2: Unstable fixed point
plot(branch_2(:, 1), branch_2(:, 2), '--', 'Color', [20,20,20]/255, 'LineWidth', 3);

% 3: Stable limit cycle
plot(branch_3(1:20:end, 1), branch_3(1:20:end, 2:3), 'o', 'MarkerFaceColor', [119,172,48]/255, 'Color', [119,172,48]/255, 'LineWidth', 3);

% xlim([24, 38])
% ylim([0.004, 1])
pbaspect([2.5 1 1])
set(gca,'fontsize',font_size);
xlabel(txt(1))
ylabel(txt(2))
box off


%% --------------- Generate final_fig. 11B for Supercritical_hopf_bifurcation_Mito_H2O2_vs_k0.xlsx period --------- %%

[num, txt] = xlsread('Supercritical_Period_vs_k0.xlsx');

font_size = 24;
figure(1104)
plot(num(1:2:end,1), num(1:2:end,2)/3600, 'Color', [119,172,48]/255, 'LineWidth', 3);

% xlim([10, 26])
% ylim([0.004, 1])
pbaspect([2.5 1 1])
set(gca,'fontsize',font_size);
xlabel(txt(1))
ylabel('Period (h)')
box off

%

%% --------------- Generate final_fig. 11C for sensitivity analysis ----------------------------------------------- %%
percent_change = 0.05; % (5%)
tspan = [0:100:3600*24*10]; 

% Default magnitude and period using H2O2_mito;
  default_param_2 = param; % default parameter for supercritical

    % run to stable oscillation
    tspan0 = [0:100:3600*24*100]; 
    y0 = cell2mat(struct2cell(init));
    [t2,y2] = ode15s('ModelC_ode', tspan0, y0, [], param);
    
    % start stable oscillation
    y0 = y2(end,:); 
    options = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [t3,y3] = ode15s('ModelC_ode', tspan, y0, options, param);

    H2O2_mito = y3(:,2);
    sign_changes_1 = sign(H2O2_mito(2:end) - H2O2_mito(1:end-1));
    sign_changes_2 = sign_changes_1(2:end) - sign_changes_1(1:end-1);
    
    index_last_peak = find(sign_changes_2 == -2, 1,'last') + 1;
    index_last_trough = find(sign_changes_2 == 2, 1,'last') + 1;
    
    % Magnitude;
    mag_ratio_0 = H2O2_mito(index_last_peak) / H2O2_mito(index_last_trough);
    mag_delta_0 = H2O2_mito(index_last_peak) - H2O2_mito(index_last_trough);

    % Period;
    index_last_two_peaks = find(sign_changes_2 == -2, 2,'last') + 1;
    period_0 = (t3(index_last_two_peaks(2)) - t3(index_last_two_peaks(1))) / 3600; % period in hours;

%     figure(11001)
%     plot(t3/3600, H2O2_mito)
%     hold on

%     figure(1101)
%     plot(t1/3600, H2O2_mito,'LineWidth', line_width)
%     set(gca, 'xtick', 0:24:240)
%     xlabel ('Time (h)')
%     ylabel ('Mito H2O2 (uM)')
%     hold on

% Increase or decrease parameter values;
param_cell = struct2cell(default_param_2);
param_names = fieldnames(default_param_2);

for i = 1:length(param_cell)
    i
    param = default_param_2;

    % Increased parameter value: magnitude and period using H2O2_mito;
        eval(strcat('param.', param_names{i}, '=', num2str(param_cell{i}), '*(1 + percent_change);'));
        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t3,y3] = ode15s('ModelC_ode', tspan, y0, options, param);
        H2O2_mito = y3(:,2);
        sign_changes_1 = sign(H2O2_mito(2:end) - H2O2_mito(1:end-1));
        sign_changes_2 = sign_changes_1(2:end) - sign_changes_1(1:end-1);
        
        index_last_peak = find(sign_changes_2 == -2, 1,'last') + 1;
        index_last_trough = find(sign_changes_2 == 2, 1,'last') + 1;       

        % Magnitude;
        mag_ratio_plus(i) = H2O2_mito(index_last_peak) / H2O2_mito(index_last_trough);
        mag_delta_plus(i) = H2O2_mito(index_last_peak) - H2O2_mito(index_last_trough);
        
        % Period;
        index_last_two_peaks = find(sign_changes_2 == -2, 2,'last') + 1;
        period_plus(i) = (t3(index_last_two_peaks(2)) - t3(index_last_two_peaks(1))) / 3600; % period in hours; 
    
%         plot(t3/3600, H2O2_mito)

    % Decreased parameter value: magnitude and period using H2O2_mito;
        eval(strcat('param.', param_names{i}, '=', num2str(param_cell{i}), '*(1 - percent_change);'));
        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t3,y3] = ode15s('ModelC_ode', tspan, y0, options, param);
        H2O2_mito = y3(:,2);

        sign_changes_1 = sign(H2O2_mito(2:end) - H2O2_mito(1:end-1));
        sign_changes_2 = sign_changes_1(2:end) - sign_changes_1(1:end-1);
        
        index_last_peak = find(sign_changes_2 == -2, 1,'last') + 1;
        index_last_trough = find(sign_changes_2 == 2, 1,'last') + 1;
        
        % Magnitude;
        mag_ratio_minus(i) = H2O2_mito(index_last_peak) / H2O2_mito(index_last_trough);
        mag_delta_minus(i) = H2O2_mito(index_last_peak) - H2O2_mito(index_last_trough);
        
        % Period;
        index_last_two_peaks = find(sign_changes_2 == -2, 2,'last') + 1;
        period_minus(i) = (t3(index_last_two_peaks(2)) - t3(index_last_two_peaks(1))) / 3600; % period in hours;
   
%         plot(t3/3600, H2O2_mito)
        
        % Caculate sensitivity coefficient (SC)

        SC_mag_ratio(i) = mean([(mag_ratio_plus(i) - mag_ratio_0) / mag_ratio_0 / percent_change,  (mag_ratio_minus(i) - mag_ratio_0) / mag_ratio_0 / (-percent_change)]);
    
        SC_mag_delta(i) = mean([(mag_delta_plus(i) - mag_delta_0) / mag_delta_0 / percent_change,  (mag_delta_minus(i) - mag_delta_0) / mag_delta_0 / (-percent_change)]);
    
        SC_period(i) = mean([(period_plus(i) - period_0) / period_0 / percent_change,  (period_minus(i) - period_0) / period_0 / (-percent_change)]);
end


% Tornado plot;
% mag_ratio sensitivity 
[~,idx] = sort(abs(SC_mag_ratio), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_mag_ratio = param_names(idx);
figure(1104)
barh(SC_mag_ratio(idx));
yticks([1:1:length(param_cell)]);
xlabel('sensitivity coefficient (mag ratio)');
yticklabels(param_names_mag_ratio);
pbaspect([0.55 1 1])


% mag_delta sensitivity 
[~,idx] = sort(abs(SC_mag_delta), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_mag_delta = param_names(idx);
figure(1105)
barh(SC_mag_delta(idx));
yticks([1:1:length(param_cell)]);
xlabel('sensitivity coefficient (mag delta)')
yticklabels(param_names_mag_delta);
pbaspect([0.55 1 1])

% period sensitivity 
[~,idx] = sort(abs(SC_period), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_period = param_names(idx);
figure(1106)
barh(SC_period(idx));
yticks([1:1:length(param_cell)])
xlabel('sensitivity coefficient (period)')
yticklabels(param_names_period);
pbaspect([0.55 1 1])

%

