%% --------------- Matlab codes of Ultrasensitivity Model --------------------------------------------------- %%
% Concentration unit: uM, time unit: second
close all
clear all
clc
tspan = [0:100:1000000];

%

%% --------------- PARAMETERS ------------------------------------------------------------------------------- %%
param.k0 = 10;
param.k1 = 50;
param.k2a = 10;
param.k2b = 1;
param.k2c = 30;         
param.Km2c = 3;    
param.k3 = 0.002;
param.k4f = 0.0014;
param.k4b = 0.001;
param.k4c = 0.006;
param.k5 = 220;        
param.PRXtot = 100;
param.SRXtot = 0.6;     
param.TRXtot = 30;      
param.TRtot = 1.38;     
param.PRXSO2H_switch = 1;
param.TRXSS_switch = 1; 
param.PRXSS_switch = 1; 
param.PS_switch = 1; 
param.PRXSH_switch = 1;
param.H2O2_switch = 1; 
param.k3_switch = 1;
default_param = param;

%

%% --------------- INITIAL CONDITION ------------------------------------------------------------------------ %%
init.H2O2 = 0.01;
init.PRXSH = param.PRXtot;
init.PRXSOH = 0;
init.PRXSS = 0;
init.PRXSO2H = 0; 
init.TRXSS = 0;
default_init = init;

%

%% --------------- Generate Fig. 2A-2I - Ultrasensitivity in response to k0 --------------------------------- %%
tspan = [0:100:1000000]; 
init = default_init;
param = default_param;
param.k0 = 0.1;
% param.k3 = 0;
increment = 1.01; 

k0_vector = [];
H2O2_vector = []; 
PRXSH_vector = [];
PRXSOH_vector = [];
PRXSS_vector = [];
PRXSO2H_vector = [];
PS_vector = [];
PRXSO2Htot_vector = [];
TRXSH_vector = [];
TRXSS_vector = [];

LRC_H2O2 = [];
LRC_PRXSH = [];
LRC_PRXSOH = [];
LRC_PRXSS = [];
LRC_PRXSO2H = [];
LRC_PS = [];
LRC_PRXSO2Htot = [];
LRC_TRXSH = [];
LRC_TRXSS = [];

while param.k0 <= 1000
        y0 = cell2mat(struct2cell(init));
        [t,y] = ode15s('PTRS_ode', tspan, y0, [], param);
        k0_vector = [k0_vector, param.k0]; 
        H2O2_vector = [H2O2_vector, y(end,1)];
        PRXSH_vector = [PRXSH_vector, y(end,2)];
        PRXSOH_vector = [PRXSOH_vector, y(end,3)];
        PRXSS_vector = [PRXSS_vector, y(end,4)];
        PRXSO2H_vector = [PRXSO2H_vector, y(end,5)];
        PS_vector = [PS_vector, param.PRXtot - y(end,2) - y(end,3) - y(end,4) - y(end,5)];  
        PRXSO2Htot_vector = [PRXSO2Htot_vector, param.PRXtot - y(end,2) - y(end,3) - y(end,4)];
        TRXSH_vector = [TRXSH_vector, param.TRXtot - y(end,6)];
        TRXSS_vector = [TRXSS_vector, y(end,6)];
        param.k0 = param.k0 * increment 
%     Visually check if steady state is reached
%     figure(1010)
%     plot(t, y(:,1))
%     xlabel ('Time')
%     ylabel ('H2O2')
%     hold on
end          

% calculate and plot local response coefficient (LRC) for H2O2 null curve
for j = 1:1:length(H2O2_vector)-1  

        delta_k0 = k0_vector(j+1) - k0_vector(j);
        PerInc_k0 = delta_k0 / k0_vector(j); 
        
        % H2O2
        delta_H2O2 = H2O2_vector(j+1) - H2O2_vector(j);
        PerInc_H2O2 = delta_H2O2 / H2O2_vector(j);
        LRC_H2O2(j) = PerInc_H2O2 / PerInc_k0;
        
        % PRXSH
        delta_PRXSH = PRXSH_vector(j+1) - PRXSH_vector(j);
        PerInc_PRXSH = delta_PRXSH / PRXSH_vector(j);
        LRC_PRXSH(j) = PerInc_PRXSH / PerInc_k0;
        
        % PRXSOH
        delta_PRXSOH = PRXSOH_vector(j+1) - PRXSOH_vector(j);
        PerInc_PRXSOH = delta_PRXSOH / PRXSOH_vector(j);
        LRC_PRXSOH(j) = PerInc_PRXSOH / PerInc_k0;
        
        % PRXSS
        delta_PRXSS = PRXSS_vector(j+1) - PRXSS_vector(j);
        PerInc_PRXSS = delta_PRXSS / PRXSS_vector(j);
        LRC_PRXSS(j) = PerInc_PRXSS / PerInc_k0;
        
        % PRXSO2H
        delta_PRXSO2H = PRXSO2H_vector(j+1) - PRXSO2H_vector(j);
        PerInc_PRXSO2H = delta_PRXSO2H / PRXSO2H_vector(j);
        LRC_PRXSO2H(j) = PerInc_PRXSO2H / PerInc_k0;

        % PS
        delta_PS = PS_vector(j+1) - PS_vector(j);
        PerInc_PS = delta_PS / PS_vector(j);
        LRC_PS(j) = PerInc_PS / PerInc_k0;

        % PRXSO2Htot
        delta_PRXSO2Htot = PRXSO2Htot_vector(j+1) - PRXSO2Htot_vector(j);
        PerInc_PRXSO2Htot = delta_PRXSO2Htot / PRXSO2Htot_vector(j);
        LRC_PRXSO2Htot(j) = PerInc_PRXSO2Htot / PerInc_k0;

        % TRXSH
        delta_TRXSH = TRXSH_vector(j+1) - TRXSH_vector(j);
        PerInc_TRXSH = delta_TRXSH / TRXSH_vector(j);
        LRC_TRXSH(j) = PerInc_TRXSH / PerInc_k0;
        
        % TRXSS
        delta_TRXSS = TRXSS_vector(j+1) - TRXSS_vector(j);
        PerInc_TRXSS = delta_TRXSS / TRXSS_vector(j);
        LRC_TRXSS(j) = PerInc_TRXSS / PerInc_k0;
end

% H2O2
figure(201)
yyaxis left
loglog(k0_vector, H2O2_vector, 'Color', [51,102,255]/255, 'LineWidth', 3);
hold on  
xlim([1, 1000])
xticks([1, 10, 100, 1000])
ylim([1e-4, 10])
yticks([1e-4, 1e-3, 1e-2, 1e-1, 1, 10])
set(gca,'YColor', [0.02,0.26,0.98]);
xlabel ('k0')
ylabel ('H2O2')      
yyaxis right
semilogx(k0_vector(1:(end-1)), LRC_H2O2, ':', 'Color', [100,100,100]/255, 'LineWidth', 2);
set(gca,'YColor','k');
ylim([0, 16])
ylabel('LRC');
hold on
box on
pbaspect([1.1 1 1])
set(gca,'Fontsize', 16);

% PRXSH
figure(202)
yyaxis left
loglog(k0_vector, PRXSH_vector, 'Color', [51,102,255]/255, 'LineWidth', 3);
hold on  
xlim([1, 1000])
xticks([1, 10, 100, 1000])
ylim([1e-2 200])
yticks([1e-2, 1e-1, 1, 10, 100])
set(gca,'YColor', [0.02,0.26,0.98]);
xlabel ('k0')
ylabel ('PRXSH')      
yyaxis right
semilogx(k0_vector(1:(end-1)), LRC_PRXSH, ':', 'Color', [100,100,100]/255, 'LineWidth', 2);
ylim([-15 2])
set(gca,'YColor','k');
ylabel('LRC');
hold on
box on
pbaspect([1.1 1 1])
set(gca,'Fontsize', 16);

% PRXSOH
figure(203)
yyaxis left
loglog(k0_vector, PRXSOH_vector, 'Color', [51,102,255]/255, 'LineWidth', 3);
hold on  
xlim([1, 1000])
xticks([1, 10, 100, 1000])
ylim([1e-1, 10])
yticks([1e-1, 1, 10])
set(gca,'YColor', [0.02,0.26,0.98]);
xlabel ('k0')
ylabel ('PRXSOH')      
yyaxis right
semilogx(k0_vector(1:(end-1)), LRC_PRXSOH, ':', 'Color', [100,100,100]/255, 'LineWidth', 2);
ylim([-2 1.5])
set(gca,'YColor','k');
ylabel('LRC');
hold on
box on
pbaspect([1.1 1 1])
set(gca,'Fontsize', 16);

% PRXSS
figure(204)
yyaxis left
loglog(k0_vector, PRXSS_vector, 'Color', [51,102,255]/255, 'LineWidth', 3);
hold on  
xlim([1, 1000])
xticks([1, 10, 100, 1000])
ylim([1e-2, 200])
yticks([1e-2, 1e-1, 1, 10, 100])
set(gca,'YColor', [0.02,0.26,0.98]);
xlabel ('k0')
ylabel ('PRXSS')      
yyaxis right
semilogx(k0_vector(1:(end-1)), LRC_PRXSS, ':', 'Color', [100,100,100]/255, 'LineWidth', 2);
set(gca,'YColor','k');
ylabel('LRC');
hold on
box on
pbaspect([1.1 1 1])
set(gca,'Fontsize', 16);

% PRXSO2H
figure(205)
yyaxis left
loglog(k0_vector, PRXSO2H_vector, 'Color', [51,102,255]/255, 'LineWidth', 3);
hold on  
xlim([1, 1000])
xticks([1, 10, 100, 1000])
ylim([1e-4, 200])
yticks([1e-4, 1e-2, 1, 100])
set(gca,'YColor', [0.02,0.26,0.98]);
xlabel ('k0')
ylabel ('PRXSO2H')      
yyaxis right
semilogx(k0_vector(1:(end-1)), LRC_PRXSO2H, ':', 'Color', [100,100,100]/255, 'LineWidth', 2);
ylim([-1, 20])
set(gca,'YColor','k');
ylabel('LRC');
hold on
box on
pbaspect([1.1 1 1])
set(gca,'Fontsize', 16);

% PS
figure(206)
yyaxis left
loglog(k0_vector, PS_vector, 'Color', [51,102,255]/255, 'LineWidth', 3);
hold on  
xlim([1, 1000])
xticks([1, 10, 100, 1000])
ylim([1e-5, 1])
yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
set(gca,'YColor', [0.02,0.26,0.98]);
xlabel ('k0')
ylabel ('PS')      
yyaxis right
semilogx(k0_vector(1:(end-1)), LRC_PS, ':', 'Color', [100,100,100]/255, 'LineWidth', 2);
ylim([-1, 20])
set(gca,'YColor','k');
ylabel('LRC');
hold on
box on
pbaspect([1.1 1 1])
set(gca,'Fontsize', 16);

% PRXSO2Htot
figure(207)
yyaxis left
loglog(k0_vector, PRXSO2Htot_vector, 'Color', [51,102,255]/255, 'LineWidth', 3);
hold on  
xlim([1, 1000])
xticks([1, 10, 100, 1000])
ylim([1e-4, 200])
yticks([1e-4, 1e-2, 1, 100])
set(gca,'YColor', [0.02,0.26,0.98]);
xlabel ('k0')
ylabel ('PRXSO2Htot')      
yyaxis right
semilogx(k0_vector(1:(end-1)), LRC_PRXSO2Htot, ':', 'Color', [100,100,100]/255, 'LineWidth', 2);
ylim([-1, 20])
set(gca,'YColor','k');
ylabel('LRC');
hold on
box on
pbaspect([1.1 1 1])
set(gca,'Fontsize', 16);

% TRXSH
figure(208)
yyaxis left
loglog(k0_vector, TRXSH_vector, 'Color', [51,102,255]/255, 'LineWidth', 3);
hold on  
xlim([1, 1000])
xticks([1, 10, 100, 1000])
ylim([1e-1, 100])
yticks([1e-1, 1, 10, 100])
set(gca,'YColor', [0.02,0.26,0.98]);
xlabel ('k0')
ylabel ('TRXSH')      
yyaxis right
semilogx(k0_vector(1:(end-1)), LRC_TRXSH, ':', 'Color', [100,100,100]/255, 'LineWidth', 2);
ylim([-50, 150])
set(gca,'YColor','k');
ylabel('LRC');
hold on
box on
pbaspect([1.1 1 1])
set(gca,'Fontsize', 16);

% TRXSS
figure(209)
yyaxis left
loglog(k0_vector, TRXSS_vector, 'Color', [51,102,255]/255, 'LineWidth', 3);
hold on  
xlim([1, 1000])
xticks([1, 10, 100, 1000])
ylim([1e-1, 100])
yticks([1e-1, 1, 10, 100])
set(gca,'YColor', [0.02,0.26,0.98]);
xlabel ('k0')
ylabel ('TRXSS')      
yyaxis right
semilogx(k0_vector(1:(end-1)), LRC_TRXSS, ':', 'Color', [100,100,100]/255, 'LineWidth', 2);
ylim([-15, 10])
set(gca,'YColor','k');
ylabel('LRC');
hold on
box on
pbaspect([1.1 1 1])
set(gca,'Fontsize', 16);

%

%% --------------- Generate Fig. 3A - Mechanism for ultrasensitivity ---------------------------------------- %%
param = default_param;
H2O2_vector = [0.0036, 0.0072, 0.0144];
param.H2O2_switch = 0; % set to 0 to clamp H2O2  
param.TRXSS_switch = 0; % set to 0 to clamp TRXSS

figure(301)
m = tiledlayout(1,1);
ax1 = axes(m);

for i = 1 : 1 : 3;
    i
    TRXSS_vector = [];
    Flux_k2b_vector = [];
    Flux_k2c_vector = [];
    init.H2O2 = H2O2_vector(i);
    init.TRXSS = 0;
    increment = 0.1;

    while init.TRXSS <= 30.0 + increment/2 
        TRXSS_vector = [TRXSS_vector, init.TRXSS];
        y0 = cell2mat(struct2cell(init));   
        [t,y] = ode15s('PTRS_ode', tspan, y0, [], param);

        Flux_k2b = param.k2b * (param.TRXtot - y(end, 6)) * y(end, 4);     
        Flux_k2b_vector =  [Flux_k2b_vector, Flux_k2b];

        Flux_k2c = param.k2c*0.5*(param.TRtot + y(end, 6) + param.Km2c - ((param.TRtot + y(end, 6) + param.Km2c)^2 - 4*param.TRtot*y(end, 6))^0.5);
        Flux_k2c_vector =  [Flux_k2c_vector, Flux_k2c];     

        init.TRXSS = init.TRXSS + increment;
    end

    plot(ax1, TRXSS_vector, Flux_k2b_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3, 'DisplayName', ['Flux k2b (H2O2 = ', num2str(init.H2O2),')']);
    ax1.XAxisLocation = 'top';
    xticks([0, 10, 20, 30])
    yticks([])
    xticklabels({'30', '20', '10', '0'})
    xlabel ('TRXSH (uM)')
    xlim([0, 30])
    ylim([0, 70])
    hold on
%     legend
    pbaspect([1.3 1 1])
    set(gca,'Fontsize', 18);

end  

ax2 = axes(m);
plot(ax2, TRXSS_vector, Flux_k2c_vector, 'Color', [255,50,50]/255, 'LineWidth', 3, 'DisplayName', 'Flux k2c');
ax2.XAxisLocation = 'bottom';
ax2.Color = 'none'; % make the second background transparent
xticks([0, 10, 20, 30])
xticklabels({'0', '10', '20', '30'})
xlabel ('TRXSS (uM)')
ylabel ('Flux (uM/S)')
xlim([0, 30])
ylim([0, 70])
pbaspect([1.3 1 1])
set(gca,'Fontsize', 18);
% legend

%

%% --------------- Generate Fig. 3B - Mechanism for ultrasensitivity ---------------------------------------- %%
init = default_init;
param = default_param;
param.H2O2_switch = 0; % set to 0 to clamp H2O2  
H2O2_vector = [];
H2O2_LRC_vector = [];
TRXSH_vector = [];   
TRXSS_vector = [];
LRC_TRXSH = [];
LRC_TRXSS = [];
increment = 1.01;
init.H2O2 = 1e-4;  

while   init.H2O2 <= 10;        
        H2O2_vector = [H2O2_vector, init.H2O2];    
        y0 = cell2mat(struct2cell(init));    
        [t,y] = ode15s('PTRS_ode', tspan, y0, [], param);     
        TRXSH_vector = [TRXSH_vector, y(end,6)];  
        TRXSS_vector = [TRXSS_vector, param.TRXtot - y(end,6)];    
        H2O2_LRC_vector = [H2O2_LRC_vector, y(end,1)];
        init.H2O2 =  init.H2O2 * increment;
end

% calculate and plot local response coefficient (LRC) for H2O2 null curve
for j = 1:1:length(H2O2_LRC_vector)-1  
        delta_H2O2 = H2O2_LRC_vector(j+1) - H2O2_LRC_vector(j);
        PerInc_H2O2 = delta_H2O2 / H2O2_LRC_vector(j);      
        % TRXSH
        delta_TRXSH = TRXSH_vector(j+1) - TRXSH_vector(j);
        PerInc_TRXSH = delta_TRXSH / TRXSH_vector(j);
        LRC_TRXSH(j) = PerInc_TRXSH / PerInc_H2O2;   
        % TRXSS
        delta_TRXSS = TRXSS_vector(j+1) - TRXSS_vector(j);
        PerInc_TRXSS = delta_TRXSS / TRXSS_vector(j);
        LRC_TRXSS(j) = PerInc_TRXSS / PerInc_H2O2;
end 

figure(302)
yyaxis left
set(gca,'YColor','k');
loglog(H2O2_vector, TRXSH_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3, 'DisplayName', 'TRXSH');
hold on
loglog(H2O2_vector, TRXSS_vector, '-', 'Color', [255,50,50]/255, 'LineWidth', 3, 'DisplayName', 'TRXSS');
xlabel ('H2O2 (uM)')
ylabel ('Concentration (uM)')
xlim([1e-4, 10])
xticks([1e-4, 1e-3, 1e-2, 1e-1, 1, 10])
ylim([0.1, 100])
legend

yyaxis right
semilogx(H2O2_vector(1:(end-1)), LRC_TRXSH, ':', 'Color', [0.47,0.67,0.19], 'LineWidth', 2, 'DisplayName', 'LRC TRXSH');
semilogx(H2O2_vector(1:(end-1)), LRC_TRXSS, ':', 'Color', [255,50,50]/255, 'LineWidth', 2, 'DisplayName', 'LRC TRXSS');
set(gca,'YColor','k');
ylim([-15, 65])
ylabel('LRC');
hold on
box on
% daspect([11.538461538461538, 35, 1])
pbaspect([1.3 1 1])
set(gca,'Fontsize', 22);
legend


xline([0.0036 0.0072 0.0144], ':', 'LineWidth', 1.5, 'HandleVisibility', 'off')

%

%% --------------- Generate Fig. 3C - Mechanism for ultrasensitivity ---------------------------------------- %%
init = default_init;
param = default_param;
param.PS_switch = 0; % set to 0 to clamp PS
param.H2O2_switch = 0; % set to 0 to clamp H2O2  
param.PRXSS_switch = 0; % set to 0 to clamp PRXSS 
param.PRXSO2H_switch = 0; % set to 0 to clamp PRXSO2H

H2O2_vector = [0.0075, 0.015, 0.030];

figure(303)
m = tiledlayout(1,1);
ax1 = axes(m);
for i = 1 : 1 : length(H2O2_vector);
    i
    init.H2O2 = H2O2_vector(i);
    init.PRXSS = 1e-1;
    PRXSS_vector = [];
    PRXSH_vector = [];
    PRXSOH_vector = [];
    Flux_k2a_vector = [];
    Flux_k2b_vector = [];
    increment = 1;
    while   init.PRXSS <= param.PRXtot + increment/2
            init.PRXSH = param.PRXtot - init.PRXSS - init.PRXSOH;
            PRXSS_vector = [PRXSS_vector, init.PRXSS];    

            y0 = cell2mat(struct2cell(init));    
            [t,y] = ode15s('PTRS_ode', tspan, y0, [], param);    

            % When PRXSS is clamped, without considering PRXSO2H and PS, PRXSH and PRXSOH need to meet the following two conditions:

            % PRXSH + PRXSOH = PRXtot - PRXSS                    (1)
            % k1 * H2O2 * PRXSH = k2a * PRXSOH                   (2)
            % So, PRXSH = PRXtot - PRXSS - PRXSOH          
            % So, k1 * H2O2 * (PRXtot - PRXSS - PRXSOH) = k2a * PRXSOH
            % So, PRXtot - PRXSS - PRXSOH = k2a * PRXSOH / (k1 * H2O2)
            % So, PRXtot - PRXSS = k2a * PRXSOH / (k1 * H2O2) + PRXSOH
            % So, PRXtot - PRXSS = PRXSOH * (k2a / (k1 * H2O2) + 1)  
            % So, PRXSOH = (PRXtot - PRXSS) / (k2a / (k1 * H2O2) + 1) 
            % So, Flux_k2a = k2a * (PRXtot - PRXSS) / (k2a / (k1 * H2O2) + 1)  

            Flux_k2a = param.k2a * (param.PRXtot - y(end, 4)) / (param.k2a / (param.k1 * y(end, 1)) + 1); 
%             Flux_k2a = param.k2a * y(end, 3); 
            Flux_k2a_vector = [Flux_k2a_vector, Flux_k2a];
         
            PRXSH = y(end, 2);   
            PRXSH_vector = [PRXSH_vector, PRXSH];

            PRXSOH = y(end, 3);  
            PRXSOH_vector = [PRXSOH_vector, PRXSOH ];

            Flux_k2b = param.k2b * (param.TRXtot - y(end, 6)) * y(end, 4);     
            Flux_k2b_vector = [Flux_k2b_vector, Flux_k2b];  

            init.PRXSS =  init.PRXSS + increment;
    end

    plot(ax1, PRXSS_vector, Flux_k2a_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3, 'DisplayName', ['Flux k2a (H2O2 = ', num2str(init.H2O2),')']);
    hold on
    ax1.XAxisLocation = 'top';
    xlim([0, 100])
    ylim([0, 100])
    xticks([0:25:100])
    xlabel ('PRXSH + PRXSOH (uM)')
    yticks([])
    xticklabels({'100', '75', '50', '25', '0'})
    pbaspect([1.3 1 1])
    set(gca,'Fontsize', 18);
%     legend
end

ax2 = axes(m);
plot(ax2, PRXSS_vector, Flux_k2b_vector, 'Color', [255,50,50]/255, 'LineWidth', 3, 'DisplayName', 'Flux k2b');
hold on
ax2.XAxisLocation = 'bottom';
ax2.Color = 'none'; % make the second background transparent
xticklabels({'0', '25', '50', '75', '100'})
xlabel ('PRXSS (uM)')
ylabel ('Flux (uM/S)')
xlim([0, 100])
ylim([0, 100])
xticks([0:25:100])
yticks([0:25:100])
pbaspect([1.3 1 1])
pbaspect([1.3 1 1])
set(gca,'Fontsize', 18);
% legend

%
%% --------------- Generate Fig. 3D - Mechanism for ultrasensitivity ---------------------------------------- %%
init = default_init;
param = default_param;
param.H2O2_switch = 0; % set to 0 to clamp H2O2  
H2O2_vector = []; 
H2O2_LRC_vector = [];
% Flux_k2b_vector = [];
% Flux_k2c_vector = [];
PRXSH_vector = [];
PRXSS_vector = [];
PRXSOH_vector = [];
LRC_PRXSH = [];
LRC_PRXSS = [];
LRC_PRXSOH = [];

increment = 1.01;
init.H2O2 = 1e-4;
while   init.H2O2 <= 10;     
        H2O2_vector = [H2O2_vector, init.H2O2];    
        y0 = cell2mat(struct2cell(init));    
        [t,y] = ode15s('PTRS_ode', tspan, y0, [], param);    
%         Flux_k2b = param.k2b * (param.TRXtot - y(end, 6)) * y(end, 4);     
%         Flux_k2b_vector =  [Flux_k2b_vector, Flux_k2b];
% 
%         Flux_k2c = param.k2c*0.5*(param.TRtot + y(end, 6) + param.Km2c - ((param.TRtot + y(end, 6) + param.Km2c)^2 - 4*param.TRtot*y(end, 6))^0.5);
%         Flux_k2c_vector =  [Flux_k2c_vector, Flux_k2c];   
        H2O2_LRC_vector = [H2O2_LRC_vector, y(end,1)];
        PRXSH_vector = [PRXSH_vector, y(end,2)];   
        PRXSOH_vector = [PRXSOH_vector, y(end,3)];   
        PRXSS_vector = [PRXSS_vector, y(end,4)];  
        init.H2O2 =  init.H2O2 * increment;
end

% calculate and plot local response coefficient (LRC) for H2O2 null curve
for j = 1:1:length(H2O2_LRC_vector)-1  
        delta_H2O2 = H2O2_LRC_vector(j+1) - H2O2_LRC_vector(j);
        PerInc_H2O2 = delta_H2O2 / H2O2_LRC_vector(j);   
        % PRXSH
        delta_PRXSH = PRXSH_vector(j+1) - PRXSH_vector(j);
        PerInc_PRXSH = delta_PRXSH / PRXSH_vector(j);
        LRC_PRXSH(j) = PerInc_PRXSH / PerInc_H2O2;
        % PRXSS
        delta_PRXSS = PRXSS_vector(j+1) - PRXSS_vector(j);
        PerInc_PRXSS = delta_PRXSS / PRXSS_vector(j);
        LRC_PRXSS(j) = PerInc_PRXSS / PerInc_H2O2;
        % PRXSOH
        delta_PRXSOH = PRXSOH_vector(j+1) - PRXSOH_vector(j);
        PerInc_PRXSOH = delta_PRXSOH / PRXSOH_vector(j);
        LRC_PRXSOH(j) = PerInc_PRXSOH / PerInc_H2O2;
end 


figure(304)
yyaxis left
set(gca,'YColor','k');
loglog(H2O2_vector, PRXSH_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3, 'DisplayName', 'PRXSH');
hold on
loglog(H2O2_vector, PRXSOH_vector, '-', 'Color', [0,0.45,0.74], 'LineWidth', 3, 'DisplayName', 'PRXSOH');
loglog(H2O2_vector, PRXSS_vector, '-', 'Color', [255,50,50]/255, 'LineWidth', 3, 'DisplayName', 'PRXSS');
xlabel ('H2O2 (uM)')
ylabel ('Concentration (uM)')
xlim([1e-4, 10])
xticks([1e-4, 1e-3, 1e-2, 1e-1, 1, 10])
ylim([0.1, 200])
yticks([0.1, 1, 10, 100])
legend

yyaxis right
semilogx(H2O2_vector(1:(end-1)), LRC_PRXSH, ':', 'Color', [0.47,0.67,0.19], 'LineWidth', 2, 'DisplayName', 'LRC PRXSH');
semilogx(H2O2_vector(1:(end-1)), LRC_PRXSOH, ':', 'Color', [0,0.45,0.74], 'LineWidth', 2, 'DisplayName', 'LRC PRXSOH');
semilogx(H2O2_vector(1:(end-1)), LRC_PRXSS, ':', 'Color', [255,50,50]/255, 'LineWidth', 2, 'DisplayName', 'LRC PRXSS');
set(gca,'YColor','k');
ylim([-40, 15])
yticks([-30, -15, 0, 15])
ylabel('LRC');
hold on
box on
% daspect([11.538461538461538, 35, 1])
pbaspect([1.3 1 1])
set(gca,'Fontsize', 22);
legend

xline([0.0075 0.015 0.03], ':', 'LineWidth', 1.5, 'HandleVisibility', 'off')

%
%% --------------- Generate Fig. 3E - Mechanism for ultrasensitivity ---------------------------------------- %%
% may not need a separate ode file;
init = default_init;
init.PRXSH = 0;
param = default_param;
param.H2O2_switch = 0; % set to 0 to clamp H2O2
param.PRXSH_switch = 0; % set to 0 to clamp PRXSH 
param.PRXSO2H_switch = 0; % set to 0 to clamp PRXSO2H  
H2O2_vector = [0.22, 0.44, 0.88];

figure(305)
m = tiledlayout(1,1);
ax1 = axes(m);

for i = 1 : 1 : length(H2O2_vector);
    i
    init.H2O2 = H2O2_vector(i);
    increment = 0.2;
    init.PRXSO2H = 1e-2;
    PRXSO2H_vector = [];
    Flux_k3_vector = [];
    Flux_k4_vector = [];

    while   init.PRXSO2H <= 100
            PRXSO2H_vector = [PRXSO2H_vector, init.PRXSO2H];  
            y0 = cell2mat(struct2cell(init));    
            [t,y] = ode15s('PTRS_ode_Fig3E', tspan, y0, [], param);    
    
    %         flux_k3 = k3*H2O2*PRXSOH
    %         flux_k4 = k4c*PS
    %         PS_switch * (PRXtot - PRXSH - PRXSOH - PRXSS - PRXSO2H);

            % When PRXSO2H is clamped, the following conditions need to be met:
            % k4f * PRXSO2H * SRX - k4b * PS = k4c * PS  (1)
            % SRX + PS = SRXtot                          (2)
            % So, SRX = SRXtot - PS
            % So, (k4c + k4b) * PS = k4f * PRXSO2H * (SRXtot - PS)
            % So, (k4c + k4b) * PS = k4f * PRXSO2H * SRXtot - k4f * PRXSO2H * PS
            % So, (k4c + k4b + k4f * PRXSO2H) * PS = k4f * PRXSO2H * SRXtot
            % So, PS = k4f * PRXSO2H * SRXtot / (k4c + k4b + k4f * PRXSO2H) 
            % So, Flux_k4c = k4c * PRXSO2H * SRXtot / ((k4c + k4b)/k4f + PRXSO2H)), which is Michaelis-Menten


            Flux_k3 = param.k3 * y(end, 1) * y(end, 3);     
            Flux_k3_vector =  [Flux_k3_vector, Flux_k3];

            Flux_k4 = param.k4c * init.PRXSO2H * param.SRXtot / ((param.k4c + param.k4b)/param.k4f + init.PRXSO2H);
            Flux_k4_vector =  [Flux_k4_vector, Flux_k4]; 

            init.PRXSO2H =  init.PRXSO2H + increment;
    end

    plot(ax1, PRXSO2H_vector, Flux_k3_vector*1000, 'Color', [0.47,0.67,0.19], 'LineWidth', 3, 'DisplayName', ['Flux k2b (H2O2 = ', num2str(init.H2O2),')']);
    ax1.XAxisLocation = 'top';
    xticks([0:25:100])
    yticks([])
    xticklabels({'100', '75', '50', '25', '0'})
    xlabel ('PRXtot-PRXSO2H (uM)')
    xlim([0, 100])
    ylim([0, 8])
    hold on
    pbaspect([1.3 1 1])
    set(gca,'Fontsize', 18);

end

ax2 = axes(m);
plot(ax2, PRXSO2H_vector, Flux_k4_vector*1000, 'Color', [255,50,50]/255, 'LineWidth', 3, 'DisplayName', 'Flux k2c');
ax2.XAxisLocation = 'bottom';
ax2.Color = 'none'; % make the second background transparent
xticks([0:25:100])
xticklabels({'0', '25', '50', '75', '100'})
xlabel ('PRXSO2H (uM)')
ylabel ('Flux (0.001 uM/S)')
xlim([0, 100])
ylim([0, 8])
pbaspect([1.3 1 1])
set(gca,'Fontsize', 18);


%
%% --------------- Generate Fig. 3F - Mechanism for ultrasensitivity ---------------------------------------- %%
init = default_init;
param = default_param;
param.H2O2_switch = 0; % set to 0 to clamp H2O2  
H2O2_vector = []; 
H2O2_LRC_vector = [];
% Flux_k2b_vector = [];
% Flux_k2c_vector = [];
PRXSH_vector = [];
PS_vector = [];
PRXSO2H_vector = [];
LRC_PRXSH = [];
LRC_PS = [];
LRC_PRXSO2H = [];

increment = 1.01;
init.H2O2 = 1e-4;
while   init.H2O2 <= 10;     
        H2O2_vector = [H2O2_vector, init.H2O2];    
        y0 = cell2mat(struct2cell(init));    
        [t,y] = ode15s('PTRS_ode', tspan, y0, [], param);    
%         Flux_k2b = param.k2b * (param.TRXtot - y(end, 6)) * y(end, 4);     
%         Flux_k2b_vector =  [Flux_k2b_vector, Flux_k2b];
% 
%         Flux_k2c = param.k2c*0.5*(param.TRtot + y(end, 6) + param.Km2c - ((param.TRtot + y(end, 6) + param.Km2c)^2 - 4*param.TRtot*y(end, 6))^0.5);
%         Flux_k2c_vector =  [Flux_k2c_vector, Flux_k2c];   
        H2O2_LRC_vector = [H2O2_LRC_vector, y(end,1)];
        PRXSH_vector = [PRXSH_vector, y(end,2)];   
        PRXSO2H_vector = [PRXSO2H_vector, y(end,5)];   
        PS_vector = [PS_vector, param.PRXtot - y(end,2) - y(end,3) - y(end,4) - y(end,5)];  
        init.H2O2 =  init.H2O2 * increment;
end

% calculate and plot local response coefficient (LRC) for H2O2 null curve
for j = 1:1:length(H2O2_LRC_vector)-1  
        delta_H2O2 = H2O2_LRC_vector(j+1) - H2O2_LRC_vector(j);
        PerInc_H2O2 = delta_H2O2 / H2O2_LRC_vector(j);   
        % PRXSH
        delta_PRXSH = PRXSH_vector(j+1) - PRXSH_vector(j);
        PerInc_PRXSH = delta_PRXSH / PRXSH_vector(j);
        LRC_PRXSH(j) = PerInc_PRXSH / PerInc_H2O2;
        % PS
        delta_PS = PS_vector(j+1) - PS_vector(j);
        PerInc_PS = delta_PS / PS_vector(j);
        LRC_PS(j) = PerInc_PS / PerInc_H2O2;
        % PRXSO2H
        delta_PRXSO2H = PRXSO2H_vector(j+1) - PRXSO2H_vector(j);
        PerInc_PRXSO2H = delta_PRXSO2H / PRXSO2H_vector(j);
        LRC_PRXSO2H(j) = PerInc_PRXSO2H / PerInc_H2O2;
end 

figure(306)
yyaxis left
set(gca,'YColor','k');
% loglog(H2O2_vector, PRXSH_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3, 'DisplayName', 'PRXSH');
loglog(H2O2_vector, PRXSO2H_vector, '-', 'Color', [0,0.45,0.74], 'LineWidth', 3, 'DisplayName', 'PRXSO2H');
hold on
loglog(H2O2_vector, PS_vector, '-', 'Color', [255,50,50]/255, 'LineWidth', 3, 'DisplayName', 'PS');
xlabel ('H2O2 (uM)')
ylabel ('Concentration (uM)')
xlim([1e-3, 10])
xticks([1e-3, 1e-2, 1e-1, 1, 10])
ylim([0.0001, 200])
yticks([0.0001, 0.01, 1, 100])
legend

yyaxis right
% semilogx(H2O2_vector(1:(end-1)), LRC_PRXSH, ':', 'Color', [0.47,0.67,0.19], 'LineWidth', 2, 'DisplayName', 'LRC PRXSH');
semilogx(H2O2_vector(1:(end-1)), LRC_PRXSO2H, ':', 'Color', [0,0.45,0.74], 'LineWidth', 2, 'DisplayName', 'LRC PRXSO2H');
semilogx(H2O2_vector(1:(end-1)), LRC_PS, ':', 'Color', [255,50,50]/255, 'LineWidth', 2, 'DisplayName', 'LRC PS');
set(gca,'YColor','k');
ylim([-1, 15])
yticks([0, 5, 10, 15])
ylabel('LRC');
hold on
box on
% daspect([11.538461538461538, 35, 1])
pbaspect([1.3 1 1])
set(gca,'Fontsize', 22);
legend

xline([0.22 0.44 0.88], ':', 'LineWidth', 1.5, 'HandleVisibility', 'off')

%


%% --------------- Generate Fig. 4 and S1 - Sensitivity analysis -------------------------------------------- %%
tspan = [0:100:1000000]; 
init = default_init;
param = default_param;
param.k0 = 0.1;
increment = 1.01; % for k0 increment

% -------------  STEP 1: calculate default LRC_H2O2, PRXSH, PRXSOH, PRXSS, PRXSO2Htot, TRXSH, TRXSS;
k0_default_vector = [];
H2O2_default_vector = []; 
PRXSH_default_vector = [];
PRXSOH_default_vector = [];
PRXSS_default_vector = [];
PRXSO2Htot_default_vector = [];
TRXSH_default_vector = [];
TRXSS_default_vector = [];
LRC_H2O2_default = [];
LRC_PRXSH_default = [];
LRC_PRXSOH_default = [];
LRC_PRXSS_default = [];
LRC_PRXSO2Htot_default = [];
LRC_TRXSH_default = [];
LRC_TRXSS_default = [];

while param.k0 <= 1000
        y0 = cell2mat(struct2cell(init));
        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t,y] = ode15s('PTRS_ode', tspan, y0, options, param);
        k0_default_vector = [k0_default_vector, param.k0]; 
        H2O2_default_vector = [H2O2_default_vector, y(end,1)];
        PRXSH_default_vector = [PRXSH_default_vector, y(end,2)];
        PRXSOH_default_vector = [PRXSOH_default_vector, y(end,3)];
        PRXSS_default_vector = [PRXSS_default_vector, y(end,4)]; 
        PRXSO2Htot_default_vector = [PRXSO2Htot_default_vector, param.PRXtot - y(end,2) - y(end,3) - y(end,4)];
        TRXSH_default_vector = [TRXSH_default_vector, param.TRXtot - y(end,6)];
        TRXSS_default_vector = [TRXSS_default_vector, y(end,6)];
        param.k0 = param.k0 * increment;
end          

% calculate and plot local response coefficient (LRC) for H2O2 null curve
for j = 1:1:length(H2O2_default_vector)-1  

        delta_k0_default = k0_default_vector(j+1) - k0_default_vector(j);
        PerInc_k0_default = delta_k0_default / k0_default_vector(j); 
        
        % H2O2_default
        delta_H2O2_default = H2O2_default_vector(j+1) - H2O2_default_vector(j);
        PerInc_H2O2_default = delta_H2O2_default / H2O2_default_vector(j);
        LRC_H2O2_default(j) = PerInc_H2O2_default / PerInc_k0_default;
        
        % PRXSH_default
        delta_PRXSH_default = PRXSH_default_vector(j+1) - PRXSH_default_vector(j);
        PerInc_PRXSH_default = delta_PRXSH_default / PRXSH_default_vector(j);
        LRC_PRXSH_default(j) = PerInc_PRXSH_default / PerInc_k0_default;
        
        % PRXSOH_default
        delta_PRXSOH_default = PRXSOH_default_vector(j+1) - PRXSOH_default_vector(j);
        PerInc_PRXSOH_default = delta_PRXSOH_default / PRXSOH_default_vector(j);
        LRC_PRXSOH_default(j) = PerInc_PRXSOH_default / PerInc_k0_default;
        
        % PRXSS_default
        delta_PRXSS_default = PRXSS_default_vector(j+1) - PRXSS_default_vector(j);
        PerInc_PRXSS_default = delta_PRXSS_default / PRXSS_default_vector(j);
        LRC_PRXSS_default(j) = PerInc_PRXSS_default / PerInc_k0_default;
       
        % PRXSO2Htot_default
        delta_PRXSO2Htot_default = PRXSO2Htot_default_vector(j+1) - PRXSO2Htot_default_vector(j);
        PerInc_PRXSO2Htot_default = delta_PRXSO2Htot_default / PRXSO2Htot_default_vector(j);
        LRC_PRXSO2Htot_default(j) = PerInc_PRXSO2Htot_default / PerInc_k0_default;

        % TRXSH_default
        delta_TRXSH_default = TRXSH_default_vector(j+1) - TRXSH_default_vector(j);
        PerInc_TRXSH_default = delta_TRXSH_default / TRXSH_default_vector(j);
        LRC_TRXSH_default(j) = PerInc_TRXSH_default / PerInc_k0_default;
        
        % TRXSS_default
        delta_TRXSS_default = TRXSS_default_vector(j+1) - TRXSS_default_vector(j);
        PerInc_TRXSS_default = delta_TRXSS_default / TRXSS_default_vector(j);
        LRC_TRXSS_default(j) = PerInc_TRXSS_default / PerInc_k0_default;
end

% Default values: H2O2, PRXSH, PRXSOH, PRXSS, PRXSO2Htot, TRXSH, TRXSS;
% ----- 01. H2O2_two_peaks_default -----
% findpeak
[LRC_H2O2_peaks_default, idx_specially_for_H2O2_default] = findpeaks(LRC_H2O2_default);
LRC_H2O2_firstpeak_default = LRC_H2O2_peaks_default(1);
LRC_H2O2_secondpeak_default = LRC_H2O2_peaks_default(2);
% idx_specially_for_H2O2
k0_at_LRC_H2O2_peaks_default = k0_default_vector(idx_specially_for_H2O2_default);
k0_threshold_H2O2_firstpeak_default = k0_at_LRC_H2O2_peaks_default(1);
k0_threshold_H2O2_secondpeak_default = k0_at_LRC_H2O2_peaks_default(2);

% ----- 02. PRXSH_two_troughs_default -----
% findtrough
% [LRC_trough_temp, idx_trough]  = findpeaks(-data);
% LRC_trough_true = - LRC_trough_temp;
[LRC_PRXSH_troughs_temp_default, idx_specially_for_PRXSH_default] = findpeaks( - LRC_PRXSH_default);
LRC_PRXSH_troughs_true_default = - LRC_PRXSH_troughs_temp_default;
LRC_PRXSH_firsttrough_default = LRC_PRXSH_troughs_true_default(1);
LRC_PRXSH_secondtrough_default = LRC_PRXSH_troughs_true_default(2);
% idx_specially_for_PRXSH
k0_at_LRC_PRXSH_troughs_default = k0_default_vector(idx_specially_for_PRXSH_default);
k0_threshold_PRXSH_firsttrough_default = k0_at_LRC_PRXSH_troughs_default(1);
k0_threshold_PRXSH_secondtrough_default = k0_at_LRC_PRXSH_troughs_default(2);

% ----- 03. PRXSOH_only_one_trough_default ----- 
% find trough(minmum)
% [LRC_trough idx_trough]  = min(data);
[LRC_PRXSOH_trough_default, idx_specially_for_PRXSOH_default] = min(LRC_PRXSOH_default);
% idx_specially_for_PRXSOH
k0_at_LRC_PRXSOH_trough_default = k0_default_vector(idx_specially_for_PRXSOH_default);
k0_threshold_PRXSOH_trough_default = k0_at_LRC_PRXSOH_trough_default;

% ----- 04. PRXSS_peak_default and PRXSS_trough_default ----- 
% PRXSS_peak_default:
[LRC_PRXSS_peak_default, idx_specially_for_PRXSS_max_default] = max(LRC_PRXSS_default);
% idx_specially_for_PRXSS_max
k0_at_LRC_PRXSS_peak_default = k0_default_vector(idx_specially_for_PRXSS_max_default);
k0_threshold_PRXSS_peak_default = k0_at_LRC_PRXSS_peak_default;
% PRXSS_trough_default
[LRC_PRXSS_trough_default, idx_specially_for_PRXSS_min_default] = min(LRC_PRXSS_default);
% idx_specially_for_PRXSS_min
k0_at_LRC_PRXSS_trough_default = k0_default_vector(idx_specially_for_PRXSS_min_default);
k0_threshold_PRXSS_trough_default = k0_at_LRC_PRXSS_trough_default;

% ----- 05. PRXSO2Htot_two_peaks_default ----- 
% findpeak
[LRC_PRXSO2Htot_peaks_default, idx_specially_for_PRXSO2Htot_default] = findpeaks(LRC_PRXSO2Htot_default);
LRC_PRXSO2Htot_peaks_default_2 = sort(LRC_PRXSO2Htot_peaks_default, 'descend'); % Has to sort to obtain the true peaks because there are minor peaks due to numerical imprecision
LRC_PRXSO2Htot_firstpeak_default = LRC_PRXSO2Htot_peaks_default_2(2);
LRC_PRXSO2Htot_secondpeak_default = LRC_PRXSO2Htot_peaks_default_2(1);
% idx_specially_for_PRXSO2Htot
idx_specially_for_PRXSO2Htot_default_2 = sort(idx_specially_for_PRXSO2Htot_default, 'descend'); % Has to sort to obtain the true peaks because there are minor peaks due to numerical imprecision
idx_specially_for_PRXSO2Htot_firstpeak_default = idx_specially_for_PRXSO2Htot_default_2(2);
idx_specially_for_PRXSO2Htot_secondpeak_default = idx_specially_for_PRXSO2Htot_default_2(1);
k0_at_LRC_PRXSO2Htot_firstpeak_default = k0_default_vector(idx_specially_for_PRXSO2Htot_firstpeak_default);
k0_at_LRC_PRXSO2Htot_secondpeak_default = k0_default_vector(idx_specially_for_PRXSO2Htot_secondpeak_default);
k0_threshold_PRXSO2Htot_firstpeak_default = k0_at_LRC_PRXSO2Htot_firstpeak_default;
k0_threshold_PRXSO2Htot_secondpeak_default = k0_at_LRC_PRXSO2Htot_secondpeak_default;

% ----- 06. TRXSH_trough_default and TRXSH_peak_default ----- 
% TRXSH_trough_default
[LRC_TRXSH_trough_default, idx_specially_for_TRXSH_min_default] = min(LRC_TRXSH_default);
% idx_specially_for_TRXSH_min
k0_at_LRC_TRXSH_trough_default = k0_default_vector(idx_specially_for_TRXSH_min_default);
k0_threshold_TRXSH_trough_default = k0_at_LRC_TRXSH_trough_default;
% TRXSH_peak_default:
[LRC_TRXSH_peak_default, idx_specially_for_TRXSH_max_default] = max(LRC_TRXSH_default);
% idx_specially_for_TRXSH_max
k0_at_LRC_TRXSH_peak_default = k0_default_vector(idx_specially_for_TRXSH_max_default);
k0_threshold_TRXSH_peak_default = k0_at_LRC_TRXSH_peak_default;

% ----- 07. TRXSS_peak_default and TRXSS_trough_default ----- 
% TRXSS_peak_default:
[LRC_TRXSS_peak_default, idx_specially_for_TRXSS_max_default] = max(LRC_TRXSS_default);
% idx_specially_for_TRXSS_max
k0_at_LRC_TRXSS_peak_default = k0_default_vector(idx_specially_for_TRXSS_max_default);
k0_threshold_TRXSS_peak_default = k0_at_LRC_TRXSS_peak_default;
% TRXSS_trough_default
[LRC_TRXSS_trough_default, idx_specially_for_TRXSS_min_default] = min(LRC_TRXSS_default);
% idx_specially_for_TRXSS_min
k0_at_LRC_TRXSS_trough_default = k0_default_vector(idx_specially_for_TRXSS_min_default);
k0_threshold_TRXSS_trough_default = k0_at_LRC_TRXSS_trough_default;





% ------------ STEP 2: increase or decrease parameter values;
param_cell = struct2cell(default_param);
param_cell_test = param_cell (2:15);
param_names = fieldnames(default_param);
param_names_test = param_names (2:15);
init = default_init;

percent_change = 0.01; % for parameter change 

for i = 1:1:length(param_cell_test)
    i
% ------ Increased parameter value;

    k0_increased_vector = [];
    H2O2_increased_vector = []; 
    PRXSH_increased_vector = [];
    PRXSOH_increased_vector = [];
    PRXSS_increased_vector = [];
    PRXSO2Htot_increased_vector = [];
    TRXSH_increased_vector = [];
    TRXSS_increased_vector = [];
    LRC_H2O2_increased = [];
    LRC_PRXSH_increased = [];
    LRC_PRXSOH_increased = [];
    LRC_PRXSS_increased = [];
    LRC_PRXSO2Htot_increased = [];
    LRC_TRXSH_increased = [];
    LRC_TRXSS_increased = [];

    param = default_param;
    param.k0 = 0.1;
    eval(strcat('param.', param_names_test{i}, '=', num2str(param_cell_test{i}), '*(1 + percent_change);'));

    while param.k0 <= 1000
        y0 = cell2mat(struct2cell(init));
        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t,y] = ode15s('PTRS_ode', tspan, y0, options, param);
        k0_increased_vector = [k0_increased_vector, param.k0]; 
        H2O2_increased_vector = [H2O2_increased_vector, y(end,1)];
        PRXSH_increased_vector = [PRXSH_increased_vector, y(end,2)];
        PRXSOH_increased_vector = [PRXSOH_increased_vector, y(end,3)];
        PRXSS_increased_vector = [PRXSS_increased_vector, y(end,4)];
        PRXSO2Htot_increased_vector = [PRXSO2Htot_increased_vector, param.PRXtot - y(end,2) - y(end,3) - y(end,4)];
        TRXSH_increased_vector = [TRXSH_increased_vector, param.TRXtot - y(end,6)];
        TRXSS_increased_vector = [TRXSS_increased_vector, y(end,6)];
        param.k0 = param.k0 * increment; 
    end          

    % calculate and plot local response coefficient (LRC) for H2O2
    for k = 1:1:length(H2O2_increased_vector)-1  
    
        delta_k0_increased = k0_increased_vector(k+1) - k0_increased_vector(k);
        PerInc_k0_increased = delta_k0_increased / k0_increased_vector(k); 
        
        % H2O2_increased
        delta_H2O2_increased = H2O2_increased_vector(k+1) - H2O2_increased_vector(k);
        PerInc_H2O2_increased = delta_H2O2_increased / H2O2_increased_vector(k);
        LRC_H2O2_increased(k) = PerInc_H2O2_increased / PerInc_k0_increased;
        
        % PRXSH_increased
        delta_PRXSH_increased = PRXSH_increased_vector(k+1) - PRXSH_increased_vector(k);
        PerInc_PRXSH_increased = delta_PRXSH_increased / PRXSH_increased_vector(k);
        LRC_PRXSH_increased(k) = PerInc_PRXSH_increased / PerInc_k0_increased;
        
        % PRXSOH_increased
        delta_PRXSOH_increased = PRXSOH_increased_vector(k+1) - PRXSOH_increased_vector(k);
        PerInc_PRXSOH_increased = delta_PRXSOH_increased / PRXSOH_increased_vector(k);
        LRC_PRXSOH_increased(k) = PerInc_PRXSOH_increased / PerInc_k0_increased;
        
        % PRXSS_increased
        delta_PRXSS_increased = PRXSS_increased_vector(k+1) - PRXSS_increased_vector(k);
        PerInc_PRXSS_increased = delta_PRXSS_increased / PRXSS_increased_vector(k);
        LRC_PRXSS_increased(k) = PerInc_PRXSS_increased / PerInc_k0_increased;
        
        % PRXSO2Htot_increased
        delta_PRXSO2Htot_increased = PRXSO2Htot_increased_vector(k+1) - PRXSO2Htot_increased_vector(k);
        PerInc_PRXSO2Htot_increased = delta_PRXSO2Htot_increased / PRXSO2Htot_increased_vector(k);
        LRC_PRXSO2Htot_increased(k) = PerInc_PRXSO2Htot_increased / PerInc_k0_increased;

        % TRXSH_increased
        delta_TRXSH_increased = TRXSH_increased_vector(k+1) - TRXSH_increased_vector(k);
        PerInc_TRXSH_increased = delta_TRXSH_increased / TRXSH_increased_vector(k);
        LRC_TRXSH_increased(k) = PerInc_TRXSH_increased / PerInc_k0_increased;
        
        % TRXSS_increased
        delta_TRXSS_increased = TRXSS_increased_vector(k+1) - TRXSS_increased_vector(k);
        PerInc_TRXSS_increased = delta_TRXSS_increased / TRXSS_increased_vector(k);
        LRC_TRXSS_increased(k) = PerInc_TRXSS_increased / PerInc_k0_increased;

    end
        

    % increased values: H2O2, PRXSH, PRXSOH, PRXSS, PRXSO2Htot, TRXSH, TRXSS;
    % ----- 01. H2O2_two_peaks_increased -----
    % findpeak
    [LRC_H2O2_peaks_increased, idx_specially_for_H2O2_increased] = findpeaks(LRC_H2O2_increased);
    LRC_H2O2_firstpeak_increased(i) = LRC_H2O2_peaks_increased(1);
    LRC_H2O2_secondpeak_increased(i) = LRC_H2O2_peaks_increased(2);
    % idx_specially_for_H2O2
    k0_at_LRC_H2O2_peaks_increased = k0_increased_vector(idx_specially_for_H2O2_increased);
    k0_threshold_H2O2_firstpeak_increased(i) = k0_at_LRC_H2O2_peaks_increased(1); 
    k0_threshold_H2O2_secondpeak_increased(i) = k0_at_LRC_H2O2_peaks_increased(2); 
    
    % ----- 02. PRXSH_two_troughs_increased -----
    % findtrough
    % [LRC_trough_temp, idx_trough]  = findpeaks(-data);
    % LRC_trough_true = - LRC_trough_temp;
    [LRC_PRXSH_troughs_temp_increased, idx_specially_for_PRXSH_increased] = findpeaks( - LRC_PRXSH_increased);
    LRC_PRXSH_troughs_true_increased = - LRC_PRXSH_troughs_temp_increased;
    LRC_PRXSH_firsttrough_increased(i) = LRC_PRXSH_troughs_true_increased(1);
    LRC_PRXSH_secondtrough_increased(i) = LRC_PRXSH_troughs_true_increased(2);
    % idx_specially_for_PRXSH
    k0_at_LRC_PRXSH_troughs_increased = k0_increased_vector(idx_specially_for_PRXSH_increased);
    k0_threshold_PRXSH_firsttrough_increased(i) = k0_at_LRC_PRXSH_troughs_increased(1);
    k0_threshold_PRXSH_secondtrough_increased(i) = k0_at_LRC_PRXSH_troughs_increased(2); 
    
    % ----- 03. PRXSOH_only_one_trough_increased ----- 
    % find trough(minmum)
    % [LRC_trough idx_trough] = min(data);
    [LRC_PRXSOH_trough_increased, idx_specially_for_PRXSOH_increased] = min(LRC_PRXSOH_increased);
    LRC_PRXSOH_trough_increased(i) = LRC_PRXSOH_trough_increased;
    % idx_specially_for_PRXSOH
    k0_at_LRC_PRXSOH_trough_increased = k0_increased_vector(idx_specially_for_PRXSOH_increased);
    k0_threshold_PRXSOH_trough_increased(i) = k0_at_LRC_PRXSOH_trough_increased;
    
    % ----- 04. PRXSS_peak_increased and PRXSS_trough_increased ----- 
    % PRXSS_peak_increased:
    [LRC_PRXSS_peak_increased, idx_specially_for_PRXSS_max_increased] = max(LRC_PRXSS_increased);
    LRC_PRXSS_peak_increased(i) = LRC_PRXSS_peak_increased;
    % idx_specially_for_PRXSS_max
    k0_at_LRC_PRXSS_peak_increased = k0_increased_vector(idx_specially_for_PRXSS_max_increased);
    k0_threshold_PRXSS_peak_increased(i) = k0_at_LRC_PRXSS_peak_increased;
    % PRXSS_trough_increased
    [LRC_PRXSS_trough_increased, idx_specially_for_PRXSS_min_increased] = min(LRC_PRXSS_increased);
    LRC_PRXSS_trough_increased(i) = LRC_PRXSS_trough_increased;
    % idx_specially_for_PRXSS_min
    k0_at_LRC_PRXSS_trough_increased = k0_increased_vector(idx_specially_for_PRXSS_min_increased);
    k0_threshold_PRXSS_trough_increased(i) = k0_at_LRC_PRXSS_trough_increased;
    
    % ----- 05. PRXSO2Htot_two_peaks_increased ----- 
    % findpeak
    [LRC_PRXSO2Htot_peaks_increased, idx_specially_for_PRXSO2Htot_increased] = findpeaks(LRC_PRXSO2Htot_increased);
    LRC_PRXSO2Htot_peaks_increased_2 = sort(LRC_PRXSO2Htot_peaks_increased, 'descend');
    LRC_PRXSO2Htot_firstpeak_increased(i)  = LRC_PRXSO2Htot_peaks_increased_2(2);
    LRC_PRXSO2Htot_secondpeak_increased(i)  = LRC_PRXSO2Htot_peaks_increased_2(1);
    % idx_specially_for_PRXSO2Htot
    idx_specially_for_PRXSO2Htot_increased_2 = sort(idx_specially_for_PRXSO2Htot_increased, 'descend');
    idx_specially_for_PRXSO2Htot_firstpeak_increased = idx_specially_for_PRXSO2Htot_increased_2(2);
    idx_specially_for_PRXSO2Htot_secondpeak_increased = idx_specially_for_PRXSO2Htot_increased_2(1);
    k0_at_LRC_PRXSO2Htot_firstpeak_increased = k0_increased_vector(idx_specially_for_PRXSO2Htot_firstpeak_increased);
    k0_at_LRC_PRXSO2Htot_secondpeak_increased = k0_increased_vector(idx_specially_for_PRXSO2Htot_secondpeak_increased);
    k0_threshold_PRXSO2Htot_firstpeak_increased(i) = k0_at_LRC_PRXSO2Htot_firstpeak_increased;
    k0_threshold_PRXSO2Htot_secondpeak_increased(i) = k0_at_LRC_PRXSO2Htot_secondpeak_increased;

    % ----- 06. TRXSH_trough_increased and TRXSH_peak_increased ----- 
    % TRXSH_trough_increased
    [LRC_TRXSH_trough_increased, idx_specially_for_TRXSH_min_increased] = min(LRC_TRXSH_increased);
    LRC_TRXSH_trough_increased(i) = LRC_TRXSH_trough_increased;
    % idx_specially_for_TRXSH_min
    k0_at_LRC_TRXSH_trough_increased = k0_increased_vector(idx_specially_for_TRXSH_min_increased);
    k0_threshold_TRXSH_trough_increased(i) = k0_at_LRC_TRXSH_trough_increased;
    % TRXSH_peak_increased:
    [LRC_TRXSH_peak_increased, idx_specially_for_TRXSH_max_increased] = max(LRC_TRXSH_increased);
    LRC_TRXSH_peak_increased(i) = LRC_TRXSH_peak_increased;
    % idx_specially_for_TRXSH_max
    k0_at_LRC_TRXSH_peak_increased = k0_increased_vector(idx_specially_for_TRXSH_max_increased);
    k0_threshold_TRXSH_peak_increased(i) = k0_at_LRC_TRXSH_peak_increased;
    
    % ----- 07. TRXSS_peak_increased and TRXSS_trough_increased ----- 
    % TRXSS_peak_increased:
    [LRC_TRXSS_peak_increased, idx_specially_for_TRXSS_max_increased] = max(LRC_TRXSS_increased);
    LRC_TRXSS_peak_increased(i) = LRC_TRXSS_peak_increased;
    % idx_specially_for_TRXSS_max
    k0_at_LRC_TRXSS_peak_increased = k0_increased_vector(idx_specially_for_TRXSS_max_increased);
    k0_threshold_TRXSS_peak_increased(i) = k0_at_LRC_TRXSS_peak_increased;
    % TRXSS_trough_increased
    [LRC_TRXSS_trough_increased, idx_specially_for_TRXSS_min_increased] = min(LRC_TRXSS_increased);
    LRC_TRXSS_trough_increased(i) = LRC_TRXSS_trough_increased;
    % idx_specially_for_TRXSS_min
    k0_at_LRC_TRXSS_trough_increased = k0_increased_vector(idx_specially_for_TRXSS_min_increased);
    k0_threshold_TRXSS_trough_increased(i) = k0_at_LRC_TRXSS_trough_increased;
    

 % ---------- Decreased parameter value;

    k0_decreased_vector = [];
    H2O2_decreased_vector = []; 
    PRXSH_decreased_vector = [];
    PRXSOH_decreased_vector = [];
    PRXSS_decreased_vector = [];
    PRXSO2Htot_decreased_vector = [];
    TRXSH_decreased_vector = [];
    TRXSS_decreased_vector = [];
    LRC_H2O2_decreased = [];
    LRC_PRXSH_decreased = [];
    LRC_PRXSOH_decreased = [];
    LRC_PRXSS_decreased = [];
    LRC_PRXSO2Htot_decreased = [];
    LRC_TRXSH_decreased = [];
    LRC_TRXSS_decreased = [];

    param = default_param;
    param.k0 = 0.1;
    eval(strcat('param.', param_names_test{i}, '=', num2str(param_cell_test{i}), '*(1 - percent_change);'));

    while param.k0 <= 1000
        y0 = cell2mat(struct2cell(init));
        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t,y] = ode15s('PTRS_ode', tspan, y0, options, param);
        k0_decreased_vector = [k0_decreased_vector, param.k0]; 
        H2O2_decreased_vector = [H2O2_decreased_vector, y(end,1)];
        PRXSH_decreased_vector = [PRXSH_decreased_vector, y(end,2)];
        PRXSOH_decreased_vector = [PRXSOH_decreased_vector, y(end,3)];
        PRXSS_decreased_vector = [PRXSS_decreased_vector, y(end,4)];
        PRXSO2Htot_decreased_vector = [PRXSO2Htot_decreased_vector, param.PRXtot - y(end,2) - y(end,3) - y(end,4)];
        TRXSH_decreased_vector = [TRXSH_decreased_vector, param.TRXtot - y(end,6)];
        TRXSS_decreased_vector = [TRXSS_decreased_vector, y(end,6)];
        param.k0 = param.k0 * increment; 
    end          


    % calculate and plot local response coefficient (LRC) for H2O2
    for M = 1:1:length(H2O2_decreased_vector)-1  
    
        delta_k0_decreased = k0_decreased_vector(M+1) - k0_decreased_vector(M);
        PerInc_k0_decreased = delta_k0_decreased / k0_decreased_vector(M); 
        
        % H2O2_decreased
        delta_H2O2_decreased = H2O2_decreased_vector(M+1) - H2O2_decreased_vector(M);
        PerInc_H2O2_decreased = delta_H2O2_decreased / H2O2_decreased_vector(M);
        LRC_H2O2_decreased(M) = PerInc_H2O2_decreased / PerInc_k0_decreased;
        
        % PRXSH_decreased
        delta_PRXSH_decreased = PRXSH_decreased_vector(M+1) - PRXSH_decreased_vector(M);
        PerInc_PRXSH_decreased = delta_PRXSH_decreased / PRXSH_decreased_vector(M);
        LRC_PRXSH_decreased(M) = PerInc_PRXSH_decreased / PerInc_k0_decreased;
        
        % PRXSOH_decreased
        delta_PRXSOH_decreased = PRXSOH_decreased_vector(M+1) - PRXSOH_decreased_vector(M);
        PerInc_PRXSOH_decreased = delta_PRXSOH_decreased / PRXSOH_decreased_vector(M);
        LRC_PRXSOH_decreased(M) = PerInc_PRXSOH_decreased / PerInc_k0_decreased;
        
        % PRXSS_decreased
        delta_PRXSS_decreased = PRXSS_decreased_vector(M+1) - PRXSS_decreased_vector(M);
        PerInc_PRXSS_decreased = delta_PRXSS_decreased / PRXSS_decreased_vector(M);
        LRC_PRXSS_decreased(M) = PerInc_PRXSS_decreased / PerInc_k0_decreased;
        
        % PRXSO2Htot_decreased
        delta_PRXSO2Htot_decreased = PRXSO2Htot_decreased_vector(M+1) - PRXSO2Htot_decreased_vector(M);
        PerInc_PRXSO2Htot_decreased = delta_PRXSO2Htot_decreased / PRXSO2Htot_decreased_vector(M);
        LRC_PRXSO2Htot_decreased(M) = PerInc_PRXSO2Htot_decreased / PerInc_k0_decreased;

        % TRXSH_decreased
        delta_TRXSH_decreased = TRXSH_decreased_vector(M+1) - TRXSH_decreased_vector(M);
        PerInc_TRXSH_decreased = delta_TRXSH_decreased / TRXSH_decreased_vector(M);
        LRC_TRXSH_decreased(M) = PerInc_TRXSH_decreased / PerInc_k0_decreased;
        
        % TRXSS_decreased
        delta_TRXSS_decreased = TRXSS_decreased_vector(M+1) - TRXSS_decreased_vector(M);
        PerInc_TRXSS_decreased = delta_TRXSS_decreased / TRXSS_decreased_vector(M);
        LRC_TRXSS_decreased(M) = PerInc_TRXSS_decreased / PerInc_k0_decreased;
    
    end

    % decreased values: H2O2, PRXSH, PRXSOH, PRXSS, PRXSO2Htot, TRXSH, TRXSS;
    % ----- 01. H2O2_two_peaks_decreased -----
    % findpeak
    [LRC_H2O2_peaks_decreased, idx_specially_for_H2O2_decreased] = findpeaks(LRC_H2O2_decreased);
    LRC_H2O2_firstpeak_decreased(i) = LRC_H2O2_peaks_decreased(1);
    LRC_H2O2_secondpeak_decreased(i) = LRC_H2O2_peaks_decreased(2);
    % idx_specially_for_H2O2
    k0_at_LRC_H2O2_peaks_decreased = k0_decreased_vector(idx_specially_for_H2O2_decreased);
    k0_threshold_H2O2_firstpeak_decreased(i) = k0_at_LRC_H2O2_peaks_decreased(1);
    k0_threshold_H2O2_secondpeak_decreased(i) = k0_at_LRC_H2O2_peaks_decreased(2); 
    
    % ----- 02. PRXSH_two_troughs_decreased -----
    % findtrough
    % [LRC_trough_temp, idx_trough]  = findpeaks(-data);
    % LRC_trough_true = - LRC_trough_temp;
    [LRC_PRXSH_troughs_temp_decreased, idx_specially_for_PRXSH_decreased] = findpeaks( - LRC_PRXSH_decreased);
    LRC_PRXSH_troughs_true_decreased = - LRC_PRXSH_troughs_temp_decreased;
    LRC_PRXSH_firsttrough_decreased(i) = LRC_PRXSH_troughs_true_decreased(1);
    LRC_PRXSH_secondtrough_decreased(i) = LRC_PRXSH_troughs_true_decreased(2);
    % idx_specially_for_PRXSH
    k0_at_LRC_PRXSH_troughs_decreased = k0_decreased_vector(idx_specially_for_PRXSH_decreased);
    k0_threshold_PRXSH_firsttrough_decreased(i) = k0_at_LRC_PRXSH_troughs_decreased(1);
    k0_threshold_PRXSH_secondtrough_decreased(i) = k0_at_LRC_PRXSH_troughs_decreased(2);
    
    % ----- 03. PRXSOH_only_one_trough_decreased ----- 
    % find trough(minmum)
    % [LRC_trough idx_trough]  = min(data);
    [LRC_PRXSOH_trough_decreased, idx_specially_for_PRXSOH_decreased] = min(LRC_PRXSOH_decreased);
    LRC_PRXSOH_trough_decreased(i) = LRC_PRXSOH_trough_decreased;
    % idx_specially_for_PRXSOH
    k0_at_LRC_PRXSOH_trough_decreased = k0_decreased_vector(idx_specially_for_PRXSOH_decreased);
    k0_threshold_PRXSOH_trough_decreased(i) = k0_at_LRC_PRXSOH_trough_decreased;  
    
    % ----- 04. PRXSS_peak_decreased and PRXSS_trough_decreased ----- 
    % PRXSS_peak_decreased:
    [LRC_PRXSS_peak_decreased, idx_specially_for_PRXSS_max_decreased] = max(LRC_PRXSS_decreased);
    LRC_PRXSS_peak_decreased(i) = LRC_PRXSS_peak_decreased;
    % idx_specially_for_PRXSS_max
    k0_at_LRC_PRXSS_peak_decreased = k0_decreased_vector(idx_specially_for_PRXSS_max_decreased);
    k0_threshold_PRXSS_peak_decreased(i) = k0_at_LRC_PRXSS_peak_decreased;
    % PRXSS_trough_decreased
    [LRC_PRXSS_trough_decreased, idx_specially_for_PRXSS_min_decreased] = min(LRC_PRXSS_decreased);
    LRC_PRXSS_trough_decreased(i) = LRC_PRXSS_trough_decreased;
    % idx_specially_for_PRXSS_min
    k0_at_LRC_PRXSS_trough_decreased = k0_decreased_vector(idx_specially_for_PRXSS_min_decreased);
    k0_threshold_PRXSS_trough_decreased(i) = k0_at_LRC_PRXSS_trough_decreased;
    
    % ----- 05. PRXSO2Htot_two_peaks_decreased ----- 
    % findpeak
    [LRC_PRXSO2Htot_peaks_decreased, idx_specially_for_PRXSO2Htot_decreased] = findpeaks(LRC_PRXSO2Htot_decreased);
    LRC_PRXSO2Htot_peaks_decreased_2 = sort(LRC_PRXSO2Htot_peaks_decreased, 'descend');
    LRC_PRXSO2Htot_firstpeak_decreased(i)  = LRC_PRXSO2Htot_peaks_decreased_2(2);
    LRC_PRXSO2Htot_secondpeak_decreased(i)  = LRC_PRXSO2Htot_peaks_decreased_2(1);
    % idx_specially_for_PRXSO2Htot
    idx_specially_for_PRXSO2Htot_decreased_2 = sort(idx_specially_for_PRXSO2Htot_decreased, 'descend');
    idx_specially_for_PRXSO2Htot_firstpeak_decreased = idx_specially_for_PRXSO2Htot_decreased_2(2);
    idx_specially_for_PRXSO2Htot_secondpeak_decreased = idx_specially_for_PRXSO2Htot_decreased_2(1);
    k0_at_LRC_PRXSO2Htot_firstpeak_decreased = k0_decreased_vector(idx_specially_for_PRXSO2Htot_firstpeak_decreased);
    k0_at_LRC_PRXSO2Htot_secondpeak_decreased = k0_decreased_vector(idx_specially_for_PRXSO2Htot_secondpeak_decreased);
    k0_threshold_PRXSO2Htot_firstpeak_decreased(i) = k0_at_LRC_PRXSO2Htot_firstpeak_decreased;
    k0_threshold_PRXSO2Htot_secondpeak_decreased(i) = k0_at_LRC_PRXSO2Htot_secondpeak_decreased;
    
    % ----- 06. TRXSH_trough_decreased and TRXSH_peak_decreased ----- 
    % TRXSH_trough_decreased
    [LRC_TRXSH_trough_decreased, idx_specially_for_TRXSH_min_decreased] = min(LRC_TRXSH_decreased);
    LRC_TRXSH_trough_decreased(i) = LRC_TRXSH_trough_decreased;
    % idx_specially_for_TRXSH_min
    k0_at_LRC_TRXSH_trough_decreased = k0_decreased_vector(idx_specially_for_TRXSH_min_decreased);
    k0_threshold_TRXSH_trough_decreased(i) = k0_at_LRC_TRXSH_trough_decreased;
    % TRXSH_peak_decreased:
    [LRC_TRXSH_peak_decreased, idx_specially_for_TRXSH_max_decreased] = max(LRC_TRXSH_decreased);
    LRC_TRXSH_peak_decreased(i) = LRC_TRXSH_peak_decreased;
    % idx_specially_for_TRXSH_max
    k0_at_LRC_TRXSH_peak_decreased = k0_decreased_vector(idx_specially_for_TRXSH_max_decreased);
    k0_threshold_TRXSH_peak_decreased(i) = k0_at_LRC_TRXSH_peak_decreased;
    
    % ----- 07. TRXSS_peak_decreased and TRXSS_trough_decreased ----- 
    % TRXSS_peak_decreased:
    [LRC_TRXSS_peak_decreased, idx_specially_for_TRXSS_max_decreased] = max(LRC_TRXSS_decreased);
    LRC_TRXSS_peak_decreased(i) = LRC_TRXSS_peak_decreased;
    % idx_specially_for_TRXSS_max
    k0_at_LRC_TRXSS_peak_decreased = k0_decreased_vector(idx_specially_for_TRXSS_max_decreased);
    k0_threshold_TRXSS_peak_decreased(i) = k0_at_LRC_TRXSS_peak_decreased;
    % TRXSS_trough_decreased
    [LRC_TRXSS_trough_decreased, idx_specially_for_TRXSS_min_decreased] = min(LRC_TRXSS_decreased);
    LRC_TRXSS_trough_decreased(i) = LRC_TRXSS_trough_decreased;
    % idx_specially_for_TRXSS_min
    k0_at_LRC_TRXSS_trough_decreased = k0_decreased_vector(idx_specially_for_TRXSS_min_decreased);
    k0_threshold_TRXSS_trough_decreased(i) = k0_at_LRC_TRXSS_trough_decreased;


% ------ Step3: Caculate sensitivity coefficient (SC) : H2O2, PRXSH, PRXSOH, PRXSS, PRXSO2Htot, TRXSH, TRXSS;  

    % H2O2_firstpeak:
    delta_LRC_H2O2_firstpeak(i) = mean([(LRC_H2O2_firstpeak_increased(i) - LRC_H2O2_firstpeak_default) / LRC_H2O2_firstpeak_default / percent_change,  (LRC_H2O2_firstpeak_decreased(i) - LRC_H2O2_firstpeak_default) / LRC_H2O2_firstpeak_default / (-percent_change)]);
    delta_k0_threshold_H2O2_firstpeak(i) = mean([(k0_threshold_H2O2_firstpeak_increased(i) - k0_threshold_H2O2_firstpeak_default) / k0_threshold_H2O2_firstpeak_default / percent_change,  (k0_threshold_H2O2_firstpeak_decreased(i) - k0_threshold_H2O2_firstpeak_default) / k0_threshold_H2O2_firstpeak_default / (-percent_change)]);
    % H2O2_secondpeak:
    delta_LRC_H2O2_secondpeak(i) = mean([(LRC_H2O2_secondpeak_increased(i) - LRC_H2O2_secondpeak_default) / LRC_H2O2_secondpeak_default / percent_change,  (LRC_H2O2_secondpeak_decreased(i) - LRC_H2O2_secondpeak_default) / LRC_H2O2_secondpeak_default / (-percent_change)]);
    delta_k0_threshold_H2O2_secondpeak(i) = mean([(k0_threshold_H2O2_secondpeak_increased(i) - k0_threshold_H2O2_secondpeak_default) / k0_threshold_H2O2_secondpeak_default / percent_change,  (k0_threshold_H2O2_secondpeak_decreased(i) - k0_threshold_H2O2_secondpeak_default) / k0_threshold_H2O2_secondpeak_default / (-percent_change)]);

    % PRXSH_firsttrough: 
    delta_LRC_PRXSH_firsttrough(i) = mean([(LRC_PRXSH_firsttrough_increased(i) - LRC_PRXSH_firsttrough_default) / LRC_PRXSH_firsttrough_default / percent_change,  (LRC_PRXSH_firsttrough_decreased(i) - LRC_PRXSH_firsttrough_default) / LRC_PRXSH_firsttrough_default / (-percent_change)]);
    delta_k0_threshold_PRXSH_firsttrough(i) = mean([(k0_threshold_PRXSH_firsttrough_increased(i) - k0_threshold_PRXSH_firsttrough_default) / k0_threshold_PRXSH_firsttrough_default / percent_change,  (k0_threshold_PRXSH_firsttrough_decreased(i) - k0_threshold_PRXSH_firsttrough_default) / k0_threshold_PRXSH_firsttrough_default / (-percent_change)]);
    % PRXSH_secondtrough:
    delta_LRC_PRXSH_secondtrough(i) = mean([(LRC_PRXSH_secondtrough_increased(i) - LRC_PRXSH_secondtrough_default) / LRC_PRXSH_secondtrough_default / percent_change,  (LRC_PRXSH_secondtrough_decreased(i) - LRC_PRXSH_secondtrough_default) / LRC_PRXSH_secondtrough_default / (-percent_change)]);
    delta_k0_threshold_PRXSH_secondtrough(i) = mean([(k0_threshold_PRXSH_secondtrough_increased(i) - k0_threshold_PRXSH_secondtrough_default) / k0_threshold_PRXSH_secondtrough_default / percent_change,  (k0_threshold_PRXSH_secondtrough_decreased(i) - k0_threshold_PRXSH_secondtrough_default) / k0_threshold_PRXSH_secondtrough_default / (-percent_change)]);

    % PRXSOH_trough: 
    delta_LRC_PRXSOH_trough(i) = mean([(LRC_PRXSOH_trough_increased(i) - LRC_PRXSOH_trough_default) / LRC_PRXSOH_trough_default / percent_change,  (LRC_PRXSOH_trough_decreased(i) - LRC_PRXSOH_trough_default) / LRC_PRXSOH_trough_default / (-percent_change)]);
    delta_k0_threshold_PRXSOH_trough(i) = mean([(k0_threshold_PRXSOH_trough_increased(i) - k0_threshold_PRXSOH_trough_default) / k0_threshold_PRXSOH_trough_default / percent_change,  (k0_threshold_PRXSOH_trough_decreased(i) - k0_threshold_PRXSOH_trough_default) / k0_threshold_PRXSOH_trough_default / (-percent_change)]);

    % PRXSS_peak: 
    delta_LRC_PRXSS_peak(i) = mean([(LRC_PRXSS_peak_increased(i) - LRC_PRXSS_peak_default) / LRC_PRXSS_peak_default / percent_change,  (LRC_PRXSS_peak_decreased(i) - LRC_PRXSS_peak_default) / LRC_PRXSS_peak_default / (-percent_change)]);
    delta_k0_threshold_PRXSS_peak(i) = mean([(k0_threshold_PRXSS_peak_increased(i) - k0_threshold_PRXSS_peak_default) / k0_threshold_PRXSS_peak_default / percent_change,  (k0_threshold_PRXSS_peak_decreased(i) - k0_threshold_PRXSS_peak_default) / k0_threshold_PRXSS_peak_default / (-percent_change)]);

    % PRXSS_trough:
    delta_LRC_PRXSS_trough(i) = mean([(LRC_PRXSS_trough_increased(i) - LRC_PRXSS_trough_default) / LRC_PRXSS_trough_default / percent_change,  (LRC_PRXSS_trough_decreased(i) - LRC_PRXSS_trough_default) / LRC_PRXSS_trough_default / (-percent_change)]);
    delta_k0_threshold_PRXSS_trough(i) = mean([(k0_threshold_PRXSS_trough_increased(i) - k0_threshold_PRXSS_trough_default) / k0_threshold_PRXSS_trough_default / percent_change,  (k0_threshold_PRXSS_trough_decreased(i) - k0_threshold_PRXSS_trough_default) / k0_threshold_PRXSS_trough_default / (-percent_change)]);

    % PRXSO2Htot_firstpeak:
    delta_LRC_PRXSO2Htot_firstpeak(i) = mean([(LRC_PRXSO2Htot_firstpeak_increased(i) - LRC_PRXSO2Htot_firstpeak_default) / LRC_PRXSO2Htot_firstpeak_default / percent_change,  (LRC_PRXSO2Htot_firstpeak_decreased(i) - LRC_PRXSO2Htot_firstpeak_default) / LRC_PRXSO2Htot_firstpeak_default / (-percent_change)]);
    delta_k0_threshold_PRXSO2Htot_firstpeak(i) = mean([(k0_threshold_PRXSO2Htot_firstpeak_increased(i) - k0_threshold_PRXSO2Htot_firstpeak_default) / k0_threshold_PRXSO2Htot_firstpeak_default / percent_change,  (k0_threshold_PRXSO2Htot_firstpeak_decreased(i) - k0_threshold_PRXSO2Htot_firstpeak_default) / k0_threshold_PRXSO2Htot_firstpeak_default / (-percent_change)]);

    % PRXSO2Htot_secondpeak:
    delta_LRC_PRXSO2Htot_secondpeak(i) = mean([(LRC_PRXSO2Htot_secondpeak_increased(i) - LRC_PRXSO2Htot_secondpeak_default) / LRC_PRXSO2Htot_secondpeak_default / percent_change,  (LRC_PRXSO2Htot_secondpeak_decreased(i) - LRC_PRXSO2Htot_secondpeak_default) / LRC_PRXSO2Htot_secondpeak_default / (-percent_change)]);
    delta_k0_threshold_PRXSO2Htot_secondpeak(i) = mean([(k0_threshold_PRXSO2Htot_secondpeak_increased(i) - k0_threshold_PRXSO2Htot_secondpeak_default) / k0_threshold_PRXSO2Htot_secondpeak_default / percent_change,  (k0_threshold_PRXSO2Htot_secondpeak_decreased(i) - k0_threshold_PRXSO2Htot_secondpeak_default) / k0_threshold_PRXSO2Htot_secondpeak_default / (-percent_change)]);

    % TRXSH_trough:
    delta_LRC_TRXSH_trough(i) = mean([(LRC_TRXSH_trough_increased(i) - LRC_TRXSH_trough_default) / LRC_TRXSH_trough_default / percent_change,  (LRC_TRXSH_trough_decreased(i) - LRC_TRXSH_trough_default) / LRC_TRXSH_trough_default / (-percent_change)]);
    delta_k0_threshold_TRXSH_trough(i) = mean([(k0_threshold_TRXSH_trough_increased(i) - k0_threshold_TRXSH_trough_default) / k0_threshold_TRXSH_trough_default / percent_change,  (k0_threshold_TRXSH_trough_decreased(i) - k0_threshold_TRXSH_trough_default) / k0_threshold_TRXSH_trough_default / (-percent_change)]);

    % TRXSH_peak:
    delta_LRC_TRXSH_peak(i) = mean([(LRC_TRXSH_peak_increased(i) - LRC_TRXSH_peak_default) / LRC_TRXSH_peak_default / percent_change,  (LRC_TRXSH_peak_decreased(i) - LRC_TRXSH_peak_default) / LRC_TRXSH_peak_default / (-percent_change)]);
    delta_k0_threshold_TRXSH_peak(i) = mean([(k0_threshold_TRXSH_peak_increased(i) - k0_threshold_TRXSH_peak_default) / k0_threshold_TRXSH_peak_default / percent_change,  (k0_threshold_TRXSH_peak_decreased(i) - k0_threshold_TRXSH_peak_default) / k0_threshold_TRXSH_peak_default / (-percent_change)]);

    % TRXSS_peak_decreased:
    delta_LRC_TRXSS_peak(i) = mean([(LRC_TRXSS_peak_increased(i) - LRC_TRXSS_peak_default) / LRC_TRXSS_peak_default / percent_change,  (LRC_TRXSS_peak_decreased(i) - LRC_TRXSS_peak_default) / LRC_TRXSS_peak_default / (-percent_change)]);
    delta_k0_threshold_TRXSS_peak(i) = mean([(k0_threshold_TRXSS_peak_increased(i) - k0_threshold_TRXSS_peak_default) / k0_threshold_TRXSS_peak_default / percent_change,  (k0_threshold_TRXSS_peak_decreased(i) - k0_threshold_TRXSS_peak_default) / k0_threshold_TRXSS_peak_default / (-percent_change)]);

    % TRXSS_trough_decreased
    delta_LRC_TRXSS_trough(i) = mean([(LRC_TRXSS_trough_increased(i) - LRC_TRXSS_trough_default) / LRC_TRXSS_trough_default / percent_change,  (LRC_TRXSS_trough_decreased(i) - LRC_TRXSS_trough_default) / LRC_TRXSS_trough_default / (-percent_change)]);
    delta_k0_threshold_TRXSS_trough(i) = mean([(k0_threshold_TRXSS_trough_increased(i) - k0_threshold_TRXSS_trough_default) / k0_threshold_TRXSS_trough_default / percent_change,  (k0_threshold_TRXSS_trough_decreased(i) - k0_threshold_TRXSS_trough_default) / k0_threshold_TRXSS_trough_default / (-percent_change)]);

end


% -------- STEP 4: tornado plot;
% H2O2
figure(4011)
[~,idx] = sort(abs(delta_LRC_H2O2_firstpeak), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_H2O2_firstpeak = param_names_test(idx);
barh(delta_LRC_H2O2_firstpeak(idx));
xlabel('sensitivity coefficient (delta LRC H2O2 firstpeak)');
yticklabels(param_names_delta_LRC_H2O2_firstpeak);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-1.2,1.2]);
xticks([-1 -0.5 0 0.5 1])
set(gca,'Fontsize', 13);
figure(40111)
[~,idx] = sort(abs(delta_k0_threshold_H2O2_firstpeak), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_H2O2_firstpeak = param_names_test(idx);
barh(delta_k0_threshold_H2O2_firstpeak(idx));
xlabel('sensitivity coefficient (k0 threshold H2O2 firstpeak)');
yticklabels(param_names_delta_LRC_H2O2_firstpeak);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-1.2,1.2]);
xticks([-1 -0.5 0 0.5 1])
set(gca,'Fontsize', 13);

figure(4012)
[~,idx] = sort(abs(delta_LRC_H2O2_secondpeak), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_H2O2_secondpeak = param_names_test(idx);
barh(delta_LRC_H2O2_secondpeak(idx));
xlabel('sensitivity coefficient (delta LRC H2O2 secondpeak)');
yticklabels(param_names_delta_LRC_H2O2_secondpeak);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-2,2]);
xticks([-2 -1 0 1 2])
set(gca,'Fontsize', 13);
figure(40122)
[~,idx] = sort(abs(delta_k0_threshold_H2O2_secondpeak), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_H2O2_secondpeak = param_names_test(idx);
barh(delta_k0_threshold_H2O2_secondpeak(idx));
xlabel('sensitivity coefficient (k0 threshold H2O2 secondpeak)');
yticklabels(param_names_delta_LRC_H2O2_secondpeak);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-1.2,1.2]);
xticks([-1 -0.5 0 0.5 1])
set(gca,'Fontsize', 13);


% PRXSH
figure(4021)
[~,idx] = sort(abs(delta_LRC_PRXSH_firsttrough), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_PRXSH_firsttrough = param_names_test(idx);
barh(delta_LRC_PRXSH_firsttrough(idx));
xlabel('sensitivity coefficient (delta LRC PRXSH firsttrough)');
yticklabels(param_names_delta_LRC_PRXSH_firsttrough);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-1.2,1.2]);
xticks([-1 -0.5 0 0.5 1])
set(gca,'Fontsize', 13);
figure(40211)
[~,idx] = sort(abs(delta_k0_threshold_PRXSH_firsttrough), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_PRXSH_firsttrough = param_names_test(idx);
barh(delta_k0_threshold_PRXSH_firsttrough(idx));
xlabel('sensitivity coefficient (k0_threshold PRXSH firsttrough)');
yticklabels(param_names_delta_LRC_PRXSH_firsttrough);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-1.2,1.2]);
xticks([-1 -0.5 0 0.5 1])
set(gca,'Fontsize', 13);

figure(4022)
[~,idx] = sort(abs(delta_LRC_PRXSH_secondtrough), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_PRXSH_secondtrough = param_names_test(idx);
barh(delta_LRC_PRXSH_secondtrough(idx));
xlabel('sensitivity coefficient (delta LRC PRXSH secondtrough)');
yticklabels(param_names_delta_LRC_PRXSH_secondtrough);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-2,2]);
xticks([-2 -1 0 1 2])
set(gca,'Fontsize', 13);
figure(40222)
[~,idx] = sort(abs(delta_k0_threshold_PRXSH_secondtrough), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_PRXSH_secondtrough = param_names_test(idx);
barh(delta_k0_threshold_PRXSH_secondtrough(idx));
xlabel('sensitivity coefficient (k0_threshold PRXSH secondtrough)');
yticklabels(param_names_delta_LRC_PRXSH_secondtrough);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-2,2]);
xticks([-2 -1 0 1 2])
set(gca,'Fontsize', 13);


% PRXSOH
figure(4031)
[~,idx] = sort(abs(delta_LRC_PRXSOH_trough), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_PRXSOH_trough = param_names_test(idx);
barh(delta_LRC_PRXSOH_trough(idx));
xlabel('sensitivity coefficient (delta LRC PRXSOH trough)');
yticklabels(param_names_delta_LRC_PRXSOH_trough);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-2,2]);
xticks([-2 -1 0 1 2])
set(gca,'Fontsize', 13);
figure(40311)
[~,idx] = sort(abs(delta_k0_threshold_PRXSOH_trough), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_PRXSOH_trough = param_names_test(idx);
barh(delta_k0_threshold_PRXSOH_trough(idx));
xlabel('sensitivity coefficient (k0_threshold PRXSOH trough)');
yticklabels(param_names_delta_LRC_PRXSOH_trough);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-2,2]);
xticks([-2 -1 0 1 2])
set(gca,'Fontsize', 13);


% PRXSS
figure(4041)
[~,idx] = sort(abs(delta_LRC_PRXSS_peak), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_PRXSS_peak = param_names_test(idx);
barh(delta_LRC_PRXSS_peak(idx));
xlabel('sensitivity coefficient (delta LRC PRXSS firstpeak)');
yticklabels(param_names_delta_LRC_PRXSS_peak);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-2.5,2.5]);
xticks([-2 -1 0 1 2])
set(gca,'Fontsize', 13);
figure(40411)
[~,idx] = sort(abs(delta_k0_threshold_PRXSS_peak), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_PRXSS_peak = param_names_test(idx);
barh(delta_k0_threshold_PRXSS_peak(idx));
xlabel('sensitivity coefficient (k0_threshold PRXSS firstpeak)');
yticklabels(param_names_delta_LRC_PRXSS_peak);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-2.5,2.5]);
xticks([-2 -1 0 1 2])
set(gca,'Fontsize', 13);

figure(4042)
[~,idx] = sort(abs(delta_LRC_PRXSS_trough), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_PRXSS_trough = param_names_test(idx);
barh(delta_LRC_PRXSS_trough(idx));
xlabel('sensitivity coefficient (delta LRC PRXSS secondtrough)');
yticklabels(param_names_delta_LRC_PRXSS_trough);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-2.8,2.8]);
xticks([-2 -1 0 1 2])
set(gca,'Fontsize', 13);
figure(40422)
[~,idx] = sort(abs(delta_k0_threshold_PRXSS_trough), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_PRXSS_trough = param_names_test(idx);
barh(delta_k0_threshold_PRXSS_trough(idx));
xlabel('sensitivity coefficient (k0_threshold PRXSS secondtrough)');
yticklabels(param_names_delta_LRC_PRXSS_trough);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-2.8,2.8]);
xticks([-2 -1 0 1 2])
set(gca,'Fontsize', 13);


% PRXSO2Htot
figure(4051)
[~,idx] = sort(abs(delta_LRC_PRXSO2Htot_firstpeak), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_PRXSO2Htot_firstpeak = param_names_test(idx);
barh(delta_LRC_PRXSO2Htot_firstpeak(idx));
xlabel('sensitivity coefficient (delta LRC PRXSO2Htot firstpeak)');
yticklabels(param_names_delta_LRC_PRXSO2Htot_firstpeak);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-1.2,1.2])
xticks([-1 -0.5 0 0.5 1])
set(gca,'Fontsize', 13);
% PRXSO2Htot
figure(40511)
[~,idx] = sort(abs(delta_k0_threshold_PRXSO2Htot_firstpeak), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_PRXSO2Htot_firstpeak = param_names_test(idx);
barh(delta_k0_threshold_PRXSO2Htot_firstpeak(idx));
xlabel('sensitivity coefficient (k0_threshold PRXSO2Htot firstpeak)');
yticklabels(param_names_delta_LRC_PRXSO2Htot_firstpeak);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-1.2,1.2])
xticks([-1 -0.5 0 0.5 1])
set(gca,'Fontsize', 13);

figure(4052)
[~,idx] = sort(abs(delta_LRC_PRXSO2Htot_secondpeak), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_PRXSO2Htot_secondpeak = param_names_test(idx);
barh(delta_LRC_PRXSO2Htot_secondpeak(idx));
xlabel('sensitivity coefficient (delta LRC PRXSO2Htot secondpeak)');
yticklabels(param_names_delta_LRC_PRXSO2Htot_secondpeak);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-2,2])
xticks([-2 -1 0 1 2])
set(gca,'Fontsize', 13);
figure(40522)
[~,idx] = sort(abs(delta_k0_threshold_PRXSO2Htot_secondpeak), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_PRXSO2Htot_secondpeak = param_names_test(idx);
barh(delta_k0_threshold_PRXSO2Htot_secondpeak(idx));
xlabel('sensitivity coefficient (k0_threshold PRXSO2Htot secondpeak)');
yticklabels(param_names_delta_LRC_PRXSO2Htot_secondpeak);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-2,2])
xticks([-2 -1 0 1 2])
set(gca,'Fontsize', 13);


% TRXSH
figure(4061)
[~,idx] = sort(abs(delta_LRC_TRXSH_trough), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_TRXSH_trough = param_names_test(idx);
barh(delta_LRC_TRXSH_trough(idx));
xlabel('sensitivity coefficient (delta LRC TRXSH firsttrough)');
yticklabels(param_names_delta_LRC_TRXSH_trough);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-2,2])
xticks([-2 -1 0 1 2])
set(gca,'Fontsize', 13);
figure(40611)
[~,idx] = sort(abs(delta_k0_threshold_TRXSH_trough), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_TRXSH_trough = param_names_test(idx);
barh(delta_k0_threshold_TRXSH_trough(idx));
xlabel('sensitivity coefficient (k0_threshold TRXSH firsttrough)');
yticklabels(param_names_delta_LRC_TRXSH_trough);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-2,2])
xticks([-2 -1 0 1 2])
set(gca,'Fontsize', 13);

figure(4062)
[~,idx] = sort(abs(delta_LRC_TRXSH_peak), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_TRXSH_peak = param_names_test(idx);
barh(delta_LRC_TRXSH_peak(idx));
xlabel('sensitivity coefficient (delta LRC TRXSH secondpeak)');
yticklabels(param_names_delta_LRC_TRXSH_peak);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-6,6]);
xticks([-6 -3 0 3 6])
set(gca,'Fontsize', 13);
figure(40622)
[~,idx] = sort(abs(delta_k0_threshold_TRXSH_peak), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_TRXSH_peak = param_names_test(idx);
barh(delta_k0_threshold_TRXSH_peak(idx));
xlabel('sensitivity coefficient (k0_threshold TRXSH secondpeak)');
yticklabels(param_names_delta_LRC_TRXSH_peak);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-6,6]);
xticks([-6 -3 0 3 6])
set(gca,'Fontsize', 13);


% TRXSS
figure(4071)
[~,idx] = sort(abs(delta_LRC_TRXSS_peak), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_TRXSS_peak = param_names_test(idx);
barh(delta_LRC_TRXSS_peak(idx));
xlabel('sensitivity coefficient (delta LRC TRXSS firstpeak)');
% xlabel(['$(\mathrm{\overline{u}})(m/s)$'],'interpreter','latex')
% xticklabels([-1 -0.5 0 0.5 1])
set(gca,'FontAngle', 'normal')
yticklabels(param_names_delta_LRC_TRXSS_peak);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-1,1])
xticks([-1 -0.5 0 0.5 1])
set(gca,'Fontsize', 13);
figure(40711)
[~,idx] = sort(abs(delta_k0_threshold_TRXSS_peak), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_TRXSS_peak = param_names_test(idx);
barh(delta_k0_threshold_TRXSS_peak(idx));
xlabel('sensitivity coefficient (k0_threshold TRXSS firstpeak)');
% xlabel(['$(\mathrm{\overline{u}})(m/s)$'],'interpreter','latex')
% xticklabels([-1 -0.5 0 0.5 1])
set(gca,'FontAngle', 'normal')
yticklabels(param_names_delta_LRC_TRXSS_peak);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-1,1])
xticks([-1 -0.5 0 0.5 1])
set(gca,'Fontsize', 13);

figure(4072)
[~,idx] = sort(abs(delta_LRC_TRXSS_trough), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_TRXSS_trough = param_names_test(idx);
barh(delta_LRC_TRXSS_trough(idx));
xlabel('sensitivity coefficient (delta LRC TRXSS secondtroughh)');
yticklabels(param_names_delta_LRC_TRXSS_trough);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-2,2])
xticks([-2 -1 0 1 2])
set(gca,'Fontsize', 13);
figure(40722)
[~,idx] = sort(abs(delta_k0_threshold_TRXSS_trough), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_delta_LRC_TRXSS_trough = param_names_test(idx);
barh(delta_k0_threshold_TRXSS_trough(idx));
xlabel('sensitivity coefficient (delta k0_threshold TRXSS secondtroughh)');
yticklabels(param_names_delta_LRC_TRXSS_trough);
pbaspect([0.6*0.57 1.2*0.61 1])
xlim([-2,2])
xticks([-2 -1 0 1 2])
set(gca,'Fontsize', 13);

%

%% --------------- Generate Fig. S2A-S2F - Time course for varying k0 --------------------------------------- %%
y0 = cell2mat(struct2cell(default_init));
param = default_param;

tspan0 = [0 : 1 : 3600*10]; 
tspan1 = [0 : 1 : 100]; 
tspan2 = [100 : 1 : 300]; 
tspan3 = [300 : 1 : 400]; 
tspan4 = [400 : 1 : 10000]; 

% (0) run to steady state;
param.k0 = 10;
[t,y] = ode15s('PTRS_ode', tspan0, y0, [], param); 
y1_0 = y(end,:);  

% run with different k0 values(20:10:80);
num_steps = 7; 
k0_vector = [];
for j = 1:1:num_steps 
    j
    % (1) run for 100 seconds with k0 = 10;
    param.k0 = 10;
    [t,y1] = ode15s('PTRS_ode', tspan1, y1_0, [], param); 
    y2_0 = y1(end,:);  

    % (2) run with different k0 values(20:10:80) between 100-300 seconds;
    param.k0 = 10 + 10*j;
    k0_vector = [k0_vector, param.k0];
    [t,y2] = ode15s('PTRS_ode', tspan2, y2_0, [], param); 
    y3_0 = y2(end,:);  

    % (3) run for 300-400 seconds with k0 = 10;
    param.k0 = 10;
    [t,y3] = ode15s('PTRS_ode', tspan3, y3_0, [], param); 
    y4_0 = y3(end,:);  

    % (4) run for 400-10000 seconds with k0 = 10 for PRXSO2H to reach steady state;
    param.k0 = 10;
    [t,y4] = ode15s('PTRS_ode', tspan4, y4_0, [], param); 
    
    % Plotting time course;
    t123 = [tspan1, tspan2, tspan3];
    y123 = [y1; y2; y3];

    t1234 = [tspan1, tspan2, tspan3, tspan4];
    y1234 = [y1; y2; y3; y4];

    pb_ratio = 1.8;
    font_size = 22;
    line_width = 3;

    % H2O2, PRXSH, PRXSOH, PRXSS, PRXSO2H, TRXSS
    figure(2001)
    semilogy(t123, y123(:,1), 'LineWidth', line_width)
    xlabel ('Time (s)')
    ylabel ('H2O2')
    pbaspect([pb_ratio 1 1])
    set(gca,'fontsize',font_size);
    hold on  

    figure(2002)
    semilogy(t123, y123(:,2), 'LineWidth', line_width)
    xlabel ('Time (s)')
    ylabel ('PRXSH')
    pbaspect([pb_ratio 1 1])
    set(gca,'fontsize',font_size);
    hold on

    figure(2003)
    semilogy(t123, y123(:,3), 'LineWidth', line_width)
    xlabel ('Time (s)')
    ylabel ('PRXSOH')
    pbaspect([pb_ratio 1 1])
    set(gca,'fontsize',font_size);
    hold on

    figure(2004)
    semilogy(t123, y123(:,4), 'LineWidth', line_width)
    xlabel ('Time (s)')
    ylabel ('PRXSS')
    pbaspect([pb_ratio 1 1])
    set(gca,'fontsize',font_size);
    hold on

    figure(2005)
    semilogy(t1234, y1234(:,5), 'LineWidth', line_width)
    xlabel ('Time (s)')
    ylabel ('PRXSO2H')
    pbaspect([pb_ratio 1 1])
    set(gca,'fontsize',font_size);
    hold on

    figure(2006)
    semilogy(t123, y123(:,6), 'LineWidth', line_width)
    xlabel ('Time (s)')
    ylabel ('TRXSS')
    pbaspect([pb_ratio 1 1])
    set(gca,'fontsize',font_size);
    hold on
 
end

%





























