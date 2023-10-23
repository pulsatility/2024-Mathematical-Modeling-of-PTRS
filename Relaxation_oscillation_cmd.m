%% Matlab code of Relaxation Oscillation Model in "Origins of Bistability and Circadian Oscillation of Cellular Hydrogen Peroxide and Hyperoxidation of Peroxiredoxin" %%
%Concentration unit: uM, time unit: second
close all
clear all
clc
%

%% --------------- PARAMETERS ---------------------------------------------------------------------------- %%
param.k0 = 30;
param.k1 = 30;
param.k21 = 10;
param.k22 = 0.21;
param.k3 = 0.052;
param.k4f = 1.4e-3;
param.k4b = 1e-3;
param.k4c = 3e-3;
param.k4d = 4.825e-5;
param.k5 = 100;
param.k6 = 1;
param.k7 = 30;
param.k8 = 8.15818e-4;
param.k9 = 4.825e-5;
param.k10 = 0.0104212;
param.k11 = 0;
param.k12f = 1.92541e-4;
param.k12b = 0;
param.k13 = 4.825e-5;
param.TRX = 10;
param.k14 = 1.92541e-4;
param.k15 = 4.825e-5;
param.k16 = 3.85e-4;
param.HSP90 = 1;
param.Vratio = 1;
param.PRXtot = 100;
default_param = param;
%

%% --------------- INITIAL CONDITION --------------------------------------------------------------------------------------------- %%
init.H2O2cyto = 0;
init.H2O2mito = 0.01;
init.SRXcyto = 35;
init.SRXOH = 0.1;
init.HSPSRXOH = 0.008;
init.SRXmito = 0.1;
init.PS = 0;
init.PRX	= 100;
init.PRXSOH = 0;
init.PRXSO2H = 0;
default_init = init;
%

%% --------------- Run to stable oscillation ------------------------------------------------------------------------------------- %%
tspan0 = [0:100:3600*24*100]; 
y0 = cell2mat(struct2cell(init));
[t,y] = ode15s('Oscillation_ode', tspan0, y0, [], param); 
%

%% --------------- Generate Fig. 8 - Time course and trajectories for relaxation oscillator -------------------------------------- %%
y0 = y(end,:); 
param = default_param;
tspan = [0:10:3600*24*5]; 
[t1,y1] = ode15s('Oscillation_ode', tspan, y0, [], param);
pb_ratio = 1.8;
font_size = 24;
x_lim = 24*5;
line_width = 2;

% Figure 8A
% H2O2mito
figure(801)
plot(t1/3600, y1(:,2),'LineWidth', line_width)
set(gca, 'xtick', 0:24:240)
xlim([0,x_lim])
xlabel ('Time (h)')
ylabel ('Mito H2O2 (uM)')
set(gca,'fontsize',font_size);
set(get(gca,'XLabel'),'FontSize',font_size);
set(get(gca,'YLabel'),'FontSize',font_size);
pbaspect([pb_ratio 1 1])
hold on

% H2O2cyto
figure(802)
plot(t1/3600, y1(:,1),'LineWidth', line_width)
set(gca, 'xtick', 0:24:240)
xlim([0,x_lim])
xlabel ('Time (h)')
ylabel ('Cyto H2O2 (uM)')
set(gca,'fontsize',font_size);
set(get(gca,'XLabel'),'FontSize',font_size);
set(get(gca,'YLabel'),'FontSize',font_size);
pbaspect([pb_ratio 1 1])
ax=gca; ax.YAxis.Exponent = -2;
hold on

% PRXSO2Htot = PRXSO2H + PS;
figure(803)
plot(t1/3600, y1(:,7) + y1(:,10),'LineWidth', line_width)
set(gca, 'xtick', 0:24:240)
xlim([0,x_lim])
xlabel ('Time (h)')
ylabel ('PRXSO2Htot (uM)')
set(gca,'fontsize',font_size);
set(get(gca,'XLabel'),'FontSize',font_size);
set(get(gca,'YLabel'),'FontSize',font_size);
pbaspect([pb_ratio 1 1])
hold on

% Non-PRXSO2H = PRX + PRXSOH + PRXSS
figure(804)
plot(t1/3600, param.PRXtot - (y1(:,7) + y1(:,10)),'LineWidth', line_width)
set(gca, 'xtick', 0:24:240)
xlim([0,x_lim])
xlabel ('Time (h)')
ylabel ('Non-PRXSO2H (uM)')
set(gca,'fontsize',font_size);
set(get(gca,'XLabel'),'FontSize',font_size);
set(get(gca,'YLabel'),'FontSize',font_size);
pbaspect([pb_ratio 1 1])
hold on

% Mito SRXtot = SRXmito + PS
figure(805)
plot(t1/3600, y1(:,6) + y1(:,7),'LineWidth', line_width)
set(gca, 'xtick', 0:24:240)
xlim([0,x_lim])
xlabel ('Time (h)')
ylabel ('Mito SRXtot (uM)')
set(gca,'fontsize',font_size);
set(get(gca,'XLabel'),'FontSize',font_size);
set(get(gca,'YLabel'),'FontSize',font_size);
pbaspect([pb_ratio 1 1])
hold on

% Cyto SRXtot = SRXcyto + SRXOH + HSPSRXOH;
figure(806)
plot(t1/3600, y1(:,3) + y1(:,4) + y1(:,5),'LineWidth', line_width)
set(gca, 'xtick', 0:24:240)
xlim([0,x_lim])
xlabel ('Time (h)')
ylabel ('Cyto SRXtot (uM)')
set(gca,'fontsize',font_size);
set(get(gca,'XLabel'),'FontSize',font_size);
set(get(gca,'YLabel'),'FontSize',font_size);
pbaspect([pb_ratio 1 1])
hold on

% Figure 8B and 8C
% Caculate normalized distance in each variable
H2O2mito = y1(:,2);
H2O2mito_max = max(H2O2mito);
dist_H2O2mito = (H2O2mito(2:end) - H2O2mito(1:end-1))/H2O2mito_max;

PRXSO2H_tot = y1(:,7) + y1(:,10);
PRXSO2H_tot_max = max(PRXSO2H_tot);
dist_PRXSO2H_tot = (PRXSO2H_tot(2:end) - PRXSO2H_tot(1:end-1))/PRXSO2H_tot_max;

SRXmito_tot = (y1(:,6) + y1(:,7));
SRXmito_tot_max = max(SRXmito_tot);
dist_SRXmito_tot = (SRXmito_tot(2:end) - SRXmito_tot(1:end-1))/SRXmito_tot_max;


% Figure 8B
figure (807)
% Overlay H2O2 null curve vs. Mito SRXtot from XPP-AUT results
[num, txt] = xlsread('XPP-AUT\Relaxation_oscillator_H2O2mito_null_curve.xlsx', 5);
column_num = length(txt);

% Obtain branch index:
branch_1_index = find(num(:, column_num) == 1);
branch_2_index = find(num(:, column_num) == 2);
branch_3_index = find(num(:, column_num) == 3);

% Obtain branch data:
branch_1 = num(branch_1_index, [1, 2]);
branch_2 = num(branch_2_index, [1, 2]);
branch_3 = num(branch_3_index, [1, 2]);

plot(branch_1(:, 1), branch_1(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
hold on
plot(branch_2(:, 1), branch_2(:, 2), ':', 'Color', [51,102,255]/255, 'LineWidth', 3);
plot(branch_3(:, 1), branch_3(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
pbaspect([1.1 1 1])
xlabel(txt(1))
ylabel(txt(2))
xlim([0, 5.5])
ylim([0, 0.3])
box on

% Overlay Mito SRXtot null curve by changing H2O2mito;  
% Needs to set H2O2mito ODE to zero in Oscillation_ode.m file;
init = default_init;
param = default_param;
Mito_SRX_tot_vector = [];
Mito_H2O2_vector = [];
init.H2O2mito = 0.001;
tspan1 = [0:100:3600*24*10]; 

while init.H2O2mito <= 10
    Mito_H2O2_vector = [Mito_H2O2_vector, init.H2O2mito];
    y0 = cell2mat(struct2cell(init));
    
    % needs to set H2O2mito ODE to zero in Oscillation_ode.m file;
    [t1,y1]=ode15s('Oscillation_ode', tspan1, y0, [], param);
    Mito_SRX_tot_vector  =  [Mito_SRX_tot_vector, y1(end,6) + y1(end,7)];
    init.H2O2mito =  init.H2O2mito * 1.1
    
    %         Visually check if steady state is reached
    %         figure(110)
    %         plot(t1, y1(:,6) + y1(:,7))
    %         xlabel ('Time')
    %         ylabel ('Mito SRXtot')
    %         hold on
end
% Remember to reset H2O2mito ODE in Oscillation_ode.m file;


figure(807)
plot(Mito_SRX_tot_vector, Mito_H2O2_vector, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 3);
xlim([0, 5.5])
ylim([0, 0.3])


% Overlay trajectory with velocity
dist = sqrt(dist_H2O2mito.^2 + dist_SRXmito_tot.^2);
% dist = sqrt(dist_H2O2mito.^2);
% dist = sqrt(dist_SRXmito_tot.^2);

dist_max = max(dist);
norm_dist = dist/dist_max; % nomalized distance
norm_dist_log10 = log10(norm_dist);

figure(807)
scatter(SRXmito_tot(1:end-1), H2O2mito(1:end-1), 50, norm_dist, 'filled')
%   scatter(SRXmito_tot(1:end-1), H2O2mito(1:end-1), 50, norm_dist_log10, 'filled')
xlabel ('Mito SRXtot (uM)')
ylabel ('Mito H2O2 (uM)')
colorbar
box on
hold on


% Figure 8C
dist = sqrt(dist_H2O2mito.^2 + dist_PRXSO2H_tot.^2 + dist_SRXmito_tot.^2);
dist_max = max(dist);
norm_dist = dist/dist_max; % nomalized distance
norm_dist_log10 = log10(norm_dist);

figure (808)
scatter3(SRXmito_tot(1:end-1), H2O2mito(1:end-1), PRXSO2H_tot(1:end-1), 50, norm_dist, 'filled') 
% scatter3(SRXmito_tot(1:end-1), H2O2mito(1:end-1), PRXSO2H_tot(1:end-1), 50, norm_dist_log10, 'filled')  
pbaspect([1 1 0.9])
xlabel ('Mito SRXtot (uM)')
ylabel ('Mito H2O2 (uM)')
zlabel ('PRXSO2Htot (uM)')
colorbar
box on
%

%% --------------- Generate Fig. 9A - Amplitude of subcritical hopf bifurcation for relaxation oscillator from XPP-AUT results --- %%
[num, txt] = xlsread('XPP-AUT\Relaxation_oscillator_amplitude_subcritical_bifurcation_H2O2mito_vs_k0.xlsx');

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

figure(901)
font_size = 24;

% 1: Stable fixed point. Index 1943 is the left hopf bifurcation point;
semilogy(branch_1(1:1943, 1), branch_1(1:1943, 2), 'Color', [255,50,50]/255, 'LineWidth', 3);
hold on

% 1: Stable fixed point. Index 1944 is the right hopf bifurcation point;
semilogy(branch_1(1944:end, 1), branch_1(1944:end, 2), 'Color', [255,50,50]/255, 'LineWidth', 3);

% 2: Unstable fixed point
semilogy(branch_2(:, 1), branch_2(:, 2), '--', 'Color', [20,20,20]/255, 'LineWidth', 3);

% 3: Stable limit cycle
semilogy(branch_3(1:3:end, 1), branch_3(1:3:end, 2:3), 'o', 'MarkerFaceColor', [119,172,48]/255, 'Color', [119,172,48]/255, 'LineWidth', 3);

% 4: Unstable limit cycle
semilogy(branch_4(1:5:end, 1), branch_4(1:5:end, 2:3),  'o', 'Color', [51,102,255]/255, 'LineWidth', 1);

xlim([24, 38])
ylim([0.004, 1])
pbaspect([2.5 1 1])
set(gca,'fontsize',font_size);
xlabel('k0 (uM/S)')
ylabel('Mito H2O2 (uM)')
box off
%

%% --------------- Generate Fig. 9B - Period of subcritical hopf bifurcation for relaxation oscillator from XPP-AUT results  ----- %%
[num, txt] = xlsread('XPP-AUT\Relaxation_oscillator_period_subcritical_bifurcation_H2O2mito_vs_k0.xlsx');

% Obtain branch index: 
branch_1_index = find(num(:, 3) == 1); 
branch_2_index = find(num(:, 3) == 2);
branch_3_index = find(num(:, 3) == 3);
branch_4_index = find(num(:, 3) == 4);
    
% Obtain branch data: 
branch_1 = num(branch_1_index, [1:2]);
branch_2 = num(branch_2_index, [1:2]);
branch_3 = num(branch_3_index, [1:2]);
branch_4 = num(branch_4_index, [1:2]);

figure(902)
font_size = 24;

% 1: Stable fixed point. Index 1943 is the left hopf bifurcation point;
% semilogy(branch_1(1:1943, 1), branch_1(1:1943, 2)/3600, 'Color', [255,50,50]/255, 'LineWidth', 3);
hold on

% 1: Stable fixed point. Index 1944 is the right hopf bifurcation point;
% semilogy(branch_1(1944:end, 1), branch_1(1944:end, 2)/3600, 'Color', [255,50,50]/255, 'LineWidth', 3);

% 2: Unstable fixed point
% semilogy(branch_2(:, 1), branch_2(:, 2)/3600, '--', 'Color', [20,20,20]/255, 'LineWidth', 3);

% 3: Stable limit cycle
semilogy(branch_3(1:end, 1), branch_3(1:end, 2)/3600, 'o', 'MarkerFaceColor', [119,172,48]/255, 'Color', [119,172,48]/255, 'LineWidth', 3);

% 4: Unstable limit cycle
semilogy(branch_4(1:end, 1), branch_4(1:end, 2)/3600,  'o', 'Color', [51,102,255]/255, 'LineWidth', 1);

xlim([24, 38])
% ylim([0.004, 1])
pbaspect([2.5 1 1])
set(gca,'fontsize',font_size);
xlabel(txt(1))
ylabel('Period (h)')
box off
%

%% --------------- Generate Fig. 9C - Sensitivity analysis ----------------------------------------------------------------------- %%
percent_change = 0.01;
y0 = y(end,:); 
tspan = [0:100:3600*24*10]; 

% Default magnitude and period of H2O2mito;
param = default_param;
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t1,y1] = ode15s('Oscillation_ode', tspan, y0, options, param);
H2O2mito = y1(:,2);
sign_changes_1 = sign(H2O2mito(2:end) - H2O2mito(1:end-1));
sign_changes_2 = sign_changes_1(2:end) - sign_changes_1(1:end-1);

index_last_peak = find(sign_changes_2 == -2, 1,'last') + 1;
index_last_trough = find(sign_changes_2 == 2, 1,'last') + 1;

% Magnitude;
mag_ratio_0 = H2O2mito(index_last_peak) / H2O2mito(index_last_trough);
mag_delta_0 = H2O2mito(index_last_peak) - H2O2mito(index_last_trough);

% Period;
index_last_two_peaks = find(sign_changes_2 == -2, 2,'last') + 1;
period_0 = (t1(index_last_two_peaks(2)) - t1(index_last_two_peaks(1))) / 3600; % period in hours;

%     figure(9031)
%     plot(t1/3600, H2O2mito,'LineWidth', line_width)
%     set(gca, 'xtick', 0:24:240)
%     xlabel ('Time (h)')
%     ylabel ('Mito H2O2 (uM)')
%     hold on


% Increase or decrease parameter values;
param_names = fieldnames(default_param);
param_names([16,18]) = []; %Exclude k11 and k12b because their default values are zero 

param_cell = struct2cell(default_param);
param_cell([16,18]) = []; %Exclude k11 and k12b because their default values are zero 

for i = 1:length(param_cell)
    i
    param = default_param;

    % Increased parameter value: magnitude and period of H2O2mito;
        eval(strcat('param.', param_names{i}, '=', num2str(param_cell{i}), '*(1 + percent_change);'));
        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t1,y1] = ode15s('Oscillation_ode', tspan, y0, options, param);
        H2O2mito = y1(:,2);
        sign_changes_1 = sign(H2O2mito(2:end) - H2O2mito(1:end-1));
        sign_changes_2 = sign_changes_1(2:end) - sign_changes_1(1:end-1);
        
        index_last_peak = find(sign_changes_2 == -2, 1,'last') + 1;
        index_last_trough = find(sign_changes_2 == 2, 1,'last') + 1;
        
        % Magnitude;
        mag_ratio_plus(i) = H2O2mito(index_last_peak) / H2O2mito(index_last_trough);
        mag_delta_plus(i) = H2O2mito(index_last_peak) - H2O2mito(index_last_trough);
        
        % Period;
        index_last_two_peaks = find(sign_changes_2 == -2, 2,'last') + 1;
        period_plus(i) = (t1(index_last_two_peaks(2)) - t1(index_last_two_peaks(1))) / 3600; % period in hours;
    

    % Decreased parameter value: magnitude and period of H2O2mito;
        eval(strcat('param.', param_names{i}, '=', num2str(param_cell{i}), '*(1 - percent_change);'));
        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t1,y1] = ode15s('Oscillation_ode', tspan, y0, options, param);
        H2O2mito = y1(:,2);
        sign_changes_1 = sign(H2O2mito(2:end) - H2O2mito(1:end-1));
        sign_changes_2 = sign_changes_1(2:end) - sign_changes_1(1:end-1);
        
        index_last_peak = find(sign_changes_2 == -2, 1,'last') + 1;
        index_last_trough = find(sign_changes_2 == 2, 1,'last') + 1;
        
        % Magnitude;
        mag_ratio_minus(i) = H2O2mito(index_last_peak) / H2O2mito(index_last_trough);
        mag_delta_minus(i) = H2O2mito(index_last_peak) - H2O2mito(index_last_trough);
        
        % Period;
        index_last_two_peaks = find(sign_changes_2 == -2, 2,'last') + 1;
        period_minus(i) = (t1(index_last_two_peaks(2)) - t1(index_last_two_peaks(1))) / 3600; % period in hours;
 
    
    % Caculate sensitivity coefficient (SC)
        SC_mag_ratio(i) = mean([(mag_ratio_plus(i) - mag_ratio_0) / mag_ratio_0 / percent_change,  (mag_ratio_minus(i) - mag_ratio_0) / mag_ratio_0 / (-percent_change)]);
    
        SC_mag_delta(i) = mean([(mag_delta_plus(i) - mag_delta_0) / mag_delta_0 / percent_change,  (mag_delta_minus(i) - mag_delta_0) / mag_delta_0 / (-percent_change)]);
    
        SC_period(i) = mean([(period_plus(i) - period_0) / period_0 / percent_change,  (period_minus(i) - period_0) / period_0 / (-percent_change)]);
end


% Tornado plot;
% Relative amplitude sensitivity
[~,idx] = sort(abs(SC_mag_ratio), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_mag_ratio = param_names(idx);
figure(903)
barh(SC_mag_ratio(idx));
ylim([0,length(param_cell)+1]);
yticks([1:1:length(param_cell)])
xlabel('sensitivity coefficient (relative amplitude)');
yticklabels(param_names_mag_ratio);
pbaspect([0.75 1 1])

% Absolute amplitude sensitivity 
[~,idx] = sort(abs(SC_mag_delta), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_mag_delta = param_names(idx);
figure(904)
barh(SC_mag_delta(idx));
ylim([0,length(param_cell)+1]);
yticks([1:1:length(param_cell)])
xlabel('sensitivity coefficient (absolute amplitude)')
yticklabels(param_names_mag_delta);
pbaspect([0.75 1 1])

% period sensitivity 
[~,idx] = sort(abs(SC_period), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_period = param_names(idx);
figure(905)
barh(SC_period(idx));
ylim([0,length(param_cell)+1]);
yticks([1:1:length(param_cell)])
xlabel('sensitivity coefficient (period)')
yticklabels(param_names_period);
pbaspect([0.75 1 1])
%
