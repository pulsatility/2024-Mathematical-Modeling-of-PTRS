%% --------------- Generate Fig. 9A and 9B --------------------------------------------- %%
% P3TRS module with k2c = 60, k5 = 110, TRtot = 5.52, favoring bistability
% Figure 9A
% k0 = 10, 20, 40, 80, 160 (number of sheets = 5)
for i = 1 : 1 : 5 
    [num, txt] = xlsread('Bifurcation_SRX_versus_H2O2_k2c_60&k5_110&TRtot_5.52', i);
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

    figure(901)
    loglog(branch_1(:, 1), branch_1(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    hold on
    loglog(branch_2(:, 1), branch_2(:, 2), ':', 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_3(:, 1), branch_3(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_4(:, 1), branch_4(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);

    pbaspect([1.2 0.9 1])
    xlabel(txt(1))
    ylabel(txt(2))
    box on

end

% Figure 9B
% k0 = 10, 20, 40, 80, 160 (number of sheets = 5)
for i = 1 : 1 : 5
    [num, txt] = xlsread('Bifurcation_SRX_versus_PRXSO2H_k2c_60&k5_110&TRtot_5.52', i);
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

    figure(902)
    loglog(branch_1(:, 1), branch_1(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    hold on
    loglog(branch_2(:, 1), branch_2(:, 2), ':', 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_3(:, 1), branch_3(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_4(:, 1), branch_4(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);

    pbaspect([1.2 0.9 1])
    xlabel(txt(1))
    ylabel(txt(2))
    box on

end

%

%% --------------- Generate Fig. S6A and S6B ------------------------------------------- %%
% P3TRS module with k2c = 30, k5 = 220, TRtot = 1.38, favoring ultrasensitivity
% Figure S6A
% k0 = 10, 20, 40, 80, 160 (number of sheets = 5)
for i = 1 : 1 : 5 
    [num, txt] = xlsread('Bifurcation_SRX_versus_H2O2_k2c_30&k5_220&TRtot_1.38', i);
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

    figure(6001)
    loglog(branch_1(:, 1), branch_1(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    hold on
    loglog(branch_2(:, 1), branch_2(:, 2), ':', 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_3(:, 1), branch_3(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_4(:, 1), branch_4(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    
    pbaspect([1.2 0.9 1])
    xlabel(txt(1))
    ylabel(txt(2))
    box on

end

% Figure S6B
% k0 = 10, 20, 40, 80, 160 (number of sheets = 5)
for i = 1 : 1 : 5
    [num, txt] = xlsread('Bifurcation_SRX_versus_PRXSO2H_k2c_30&k5_220&TRtot_1.38', i);
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

    figure(6002)
    loglog(branch_1(:, 1), branch_1(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    hold on
    loglog(branch_2(:, 1), branch_2(:, 2), ':', 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_3(:, 1), branch_3(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_4(:, 1), branch_4(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);

    pbaspect([1.2 0.9 1])
    xlabel(txt(1))
    ylabel(txt(2))
    box on

end

%% --------------- Matlab code for Fig. 9C, 9D, S6C, and S6D --------------------------- %%
% Concentration unit: uM, time unit: second
close all
clear all
clc
tspan = [0:100:1000000];

% -------------- PARAMETERS --------------
% P3TRS module with k2c = 60, k5 = 110, k6 = 100, k7 = 3000, TRtot = 5.52
param.k0 = 10;
param.k1 = 20;
param.k2a = 22;
param.k2b = 1;
param.k2c = 60;
param.k3 = 0.012;
param.k4f = 0.0014;
param.k4b = 0.001;
param.k4c = 0.006;
param.k5 = 110;
param.PRXtot = 100;
param.TRXtot = 30;
param.TRtot = 5.52;
param.Km2c = 3;
param.V_ratio = 0.3;
param.k6 = 100;
param.k7 = 3000;
param.SRXtot = 0.3;
default_param = param;

% --------------- INITIAL CONDITION -------------
init.H2O2 = 0.01;
init.H2O2_cyto = 0;
init.PRXSH = param.PRXtot;
init.PRXSOH = 0;
init.PRXSS = 0;
init.PRXSO2H = 0; 
init.TRXSS = 0;
default_init = init;

%

%% --------------- Generate Fig. 9C and 9D --------------------------------------------- %%
y0 = cell2mat(struct2cell(default_init));
param = default_param;
param.k0 = 40; % run with different param.k0 values (40, 80, 160);
param.k2c = 60;
param.k5 = 110;
param.TRtot = 5.52; 

tspan0 = [0 : 1 : 3600*100]; 
tspan1 = [0 : 1 : 3600*100]; 
tspan2 = [3600*100 : 1 : 3600*300]; 
tspan3 = [3600*300 : 1 : 3600*400]; 

% (0) run to steady state;
param.SRXtot = 5;
[t,y] = ode15s('P3TRS_ode', tspan0, y0, [], param); 
y1_0 = y(end,:);  

% run with different param.SRXtot values;
SRXtot_vector = [5, 0.5, 0.25, 0.1, 0.05, 0.01, 0.005, 0];
num_steps = length(SRXtot_vector); 

for j = 1:1:num_steps 
    j
    % (1) run between tspan1;
    param.SRXtot = 5;
    [t,y1] = ode15s('P3TRS_ode', tspan1, y1_0, [], param); 
    y2_0 = y1(end,:);  

    % (2) run with different param.SRXtot values between tspan2;
    param.SRXtot = SRXtot_vector(j);
    [t,y2] = ode15s('P3TRS_ode', tspan2, y2_0, [], param); 
    y3_0 = y2(end,:);  

    % (3) run between tspan3;
    param.SRXtot = 5;
    [t,y3] = ode15s('P3TRS_ode', tspan3, y3_0, [], param); 
    y4_0 = y3(end,:);  

    % Plotting time course;
    t123 = [tspan1, tspan2, tspan3];
    y123 = [y1; y2; y3];

    pb_ratio = 1.8;
    font_size = 22;
    line_width = 3;

    % H2O2
    figure(903)
    plot(t123/3600, y123(:,1), 'LineWidth', line_width)
    xlabel ('Time (h)')
    ylabel ('H2O2')
    pbaspect([pb_ratio 1 1])
    set(gca,'fontsize', font_size);
    hold on  

    % PRXSO2H
    figure(904)
    plot(t123/3600, y123(:,6), 'LineWidth', line_width)
    xlabel ('Time (h)')
    ylabel ('PRXSO2H')
    pbaspect([pb_ratio 1 1])
    set(gca,'fontsize', font_size);
    hold on  
    
end

%

%% --------------- Generate Fig. S6C and S6D ------------------------------------------- %%
y0 = cell2mat(struct2cell(default_init));
param = default_param;
param.k0 = 40; % run with different param.k0 values (40, 80, 160);
param.k2c = 30;
param.k5 = 220;
param.TRtot = 1.38; 

tspan0 = [0 : 1 : 3600*100]; 
tspan1 = [0 : 1 : 3600*100]; 
tspan2 = [3600*100 : 1 : 3600*300]; 
tspan3 = [3600*300 : 1 : 3600*400]; 

% (0) run to steady state;
param.SRXtot = 5;
[t,y] = ode15s('P3TRS_ode', tspan0, y0, [], param); 
y1_0 = y(end,:);  

% run with different param.SRXtot values;
SRXtot_vector = [5, 0.5, 0.25, 0.1, 0.05, 0.01, 0.005, 0];
num_steps = length(SRXtot_vector); 

for j = 1:1:num_steps 
    j
    % (1) 
    param.SRXtot = 5;
    [t,y1] = ode15s('P3TRS_ode', tspan1, y1_0, [], param); 
    y2_0 = y1(end,:);  

    % (2) run with different param.SRXtot values 
    param.SRXtot = SRXtot_vector(j);
    [t,y2] = ode15s('P3TRS_ode', tspan2, y2_0, [], param); 
    y3_0 = y2(end,:);  

    % (3) 
    param.SRXtot = 5;
    [t,y3] = ode15s('P3TRS_ode', tspan3, y3_0, [], param); 
    y4_0 = y3(end,:);  

    % Plotting time course;
    t123 = [tspan1, tspan2, tspan3];
    y123 = [y1; y2; y3];

    pb_ratio = 1.8;
    font_size = 22;
    line_width = 3;

    % H2O2
    figure(6003)
    plot(t123/3600, y123(:,1), 'LineWidth', line_width)
    xlabel ('Time (h)')
    ylabel ('H2O2')
    pbaspect([pb_ratio 1 1])
    set(gca,'fontsize',font_size);
    hold on  

    % PRXSO2H
    figure(6004)
    plot(t123/3600, y123(:,6), 'LineWidth', line_width)
    xlabel ('Time (h)')
    ylabel ('PRXSO2H')
    pbaspect([pb_ratio 1 1])
    set(gca,'fontsize',font_size);
    hold on  
    
end

%