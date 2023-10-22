%% Matlab code of Bistability Model in "Origins of Bistability and Circadian Oscillation of Cellular Hydrogen Peroxide and Hyperoxidation of Peroxiredoxin" %%
%Concentration unit: uM, time unit: second
close all
clear all
clc
tspan = [0:100:1000000];
%

%% --------------- PARAMETERS ------------------------------------------------------------------------------------------------ %%
param.k0 = 10;
param.k1 = 30;
param.k21 = 10;
param.k22 = 0.21;
param.k3 = 0.012;
param.k4f = 0.0014;
param.k4b = 0.001;
param.k4c = 0.003;
param.k5 = 50;
param.PRXtot = 100;
param.SRXtot = 0.5;
param.TRX = 10;
param.H2O2_null_curve_switch = 1; % set to 0 to generate H2O2 null curve
param.PRXSO2H_null_curve_switch = 1; % set to 0 to generate PRXSO2H null curve
param.k3_switch = 1; % set to 0 to remove H2O2 from k3 step
default_param = param;
%

%% --------------- INITIAL CONDITION ----------------------------------------------------------------------------------------- %%
init.H2O2 = 0.01;
init.PRX = param.PRXtot;
init.PRXSOH = 0;
init.PRXSS = 0;
init.PRXSO2H = 0; 
default_init = init;
%

%% --------------- Generate Fig. 2 - Saddle-node bifurcation from XPP-AUT results ------------------------------------------------------------------- %%
[num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_k0.xlsx'); %Bifurcation_k0.xlsx contains bufucation results obtained from XPP-AUT
column_num = length(txt);

% Obtain branch index: 
branch_1_index = find(num(:, column_num) == 1); 
branch_2_index = find(num(:, column_num) == 2);
branch_3_index = find(num(:, column_num) == 3);
for i = 1 : 1 : column_num - 2
    % Obtain branch data: 
    branch_1 = num(branch_1_index, [1, 1+i]);
    branch_2 = num(branch_2_index, [1, 1+i]);
    branch_3= num(branch_3_index, [1, 1+i]);
    % Plot bifurcation: 
    figure(200+i)
    loglog(branch_1(:, 1), branch_1(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    hold on
    loglog(branch_2(:, 1), branch_2(:, 2), ':', 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_3(:, 1), branch_3(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    xlim([1, 1000])
    pbaspect([1.5 1 1])
    xlabel(txt(1))
    ylabel(txt(1+i))
    box on
end
%

%% --------------- Generate Fig. 3A - H2O2 null curves ------------------------------------------------------------------------- %%
init = default_init;
param = default_param;
param.H2O2_null_curve_switch = 0;
k0_vector = [15.5763, 27, 44.6811];
LRC_H2O2 = [];
increment = 1.01;

for i = 1 : 1 : 3
        param.k0 = k0_vector(i);
        PRXSO2H_vector = [];
        H2O2_vector = [];
        init.PRXSO2H = 0.1;
        
                while init.PRXSO2H <= 1000
                    init.PRX = param.PRXtot - init.PRXSO2H;
                    PRXSO2H_vector = [PRXSO2H_vector, init.PRXSO2H];
                    y0 = cell2mat(struct2cell(init));
                    [t,y]=ode15s('Bistability_ode', tspan, y0, [], param);
                    H2O2_vector =  [H2O2_vector, y(end,1)];
                    init.PRXSO2H =  init.PRXSO2H * increment

                    % Visually check if steady state is reached
                    % figure(3010)
                    % plot(t, y(:,1))
                    % xlabel ('Time')
                    % ylabel ('H2O2')
                    % hold on
                end
                
        % find the first index where PRXSO2H >= 100 
        index_100 = find(PRXSO2H_vector >= 100, 1); 

        figure(301)
        loglog(PRXSO2H_vector(1 : index_100), H2O2_vector(1 : index_100),  'Color', [51,102,255]/255, 'LineWidth', 2);
        hold on  
        loglog(PRXSO2H_vector(index_100 : end), H2O2_vector(index_100 : end), ':', 'Color', [51,102,255]/255, 'LineWidth', 2);        
        xlabel ('PRXSO2H')
        ylabel ('H2O2')      
        pbaspect([1.3 1 1])
                          
        % calculate and plot local response coefficient (LRC) for H2O2 null curve     
        for j= 2:1:length(H2O2_vector)   % skipping the first concentration where H2O2 could be 0
            if j==length(H2O2_vector)    % stop when reaching the last concentration because no next concentration avaiable to calculate the delta
                break
            end
            delta_H2O2=H2O2_vector(j+1) - H2O2_vector(j);
            PerInc_H2O2=delta_H2O2 / H2O2_vector(j);

            delta_PRXSO2H=PRXSO2H_vector(j+1) - PRXSO2H_vector(j);
            PerInc_PRXSO2H=delta_PRXSO2H / PRXSO2H_vector(j);    

            LRC_H2O2(i, j-1) = PerInc_H2O2 / PerInc_PRXSO2H;
        end
        
        [LRC_H2O2_Max(i), I_PRXSO2H(i)] = max(LRC_H2O2(i,:));   
        PRXSO2H_value_for_LRC_H2O2_Max(i) = PRXSO2H_vector(1+I_PRXSO2H(i));

        figure(3011)
        lineH2O2 = plot(PRXSO2H_vector(2:(end-1)), LRC_H2O2(i,:), 'LineWidth', 2);
        xlabel('PRXSO2H');
        ylabel('LRC of H2O2');
        hold on
        box on
        pbaspect([1.3 1 1])
        %set(gca,'Fontsize',28);
end
%

%% --------------- Generate Fig. 3A - PRXSO2H null curve ----------------------------------------------------------------------- %%
init = default_init;
param = default_param;
param.PRXSO2H_null_curve_switch = 0;
PRXSO2H_vector = [];
H2O2_vector = [];
init.H2O2 = 0.001;
LRC_PRXSO2H = [];
increment = 1.01;

while init.H2O2 <= 100    
        init.PRX = param.PRXtot - init.PRXSO2H;
        H2O2_vector = [H2O2_vector, init.H2O2];
        y0 = cell2mat(struct2cell(init));
        [t,y]=ode15s('Bistability_ode', tspan, y0, [], param);
        PRXSO2H_vector =  [PRXSO2H_vector, y(end,5)];
        init.H2O2 =  init.H2O2 * increment

%         #Visually check if steady state is reached
%         figure(3012)
%         plot(t, y(:,5))
%         xlabel ('Time')
%         ylabel ('PRXSO2H')
%         hold on   
end

figure(301)
loglog(PRXSO2H_vector, H2O2_vector, 'Color', [255,50,50]/255, 'LineWidth', 2);
xlabel ('PRXSO2H')
ylabel ('H2O2')
hold on
xlim([0.1, 300])
ylim([0.001, 5])
pbaspect([1.3 1 1])
set(gca,'Fontsize',10);

% calculate and plot LRC for PRXSO2H null curve
for j= 2:1:length(PRXSO2H_vector)   % skipping the first concentration where PRXSO2H could be 0
    if j==length(PRXSO2H_vector)    % stop when reaching the last concentration becasue no next concentration avaiable to calculate the delta
        break
    end
    delta_PRXSO2H=PRXSO2H_vector(j+1) - PRXSO2H_vector(j);
    PerInc_PRXSO2H=delta_PRXSO2H / PRXSO2H_vector(j);    
    delta_H2O2=H2O2_vector(j+1) - H2O2_vector(j);
    PerInc_H2O2=delta_H2O2 / H2O2_vector(j);   
    LRC_PRXSO2H(j-1) = PerInc_PRXSO2H / PerInc_H2O2;  
end

[LRC_PRXSO2H_Max, I_H2O2] = max(LRC_PRXSO2H);
H2O2_value_for_LRC_PRXSO2H_Max = H2O2_vector(1+I_H2O2);

figure(3013)
linePRXSO2H = semilogx(H2O2_vector(2:(end-1)), LRC_PRXSO2H, 'LineWidth', 2);
xlabel('H2O2');
ylabel('LRC of PRXSO2H');
hold on
box on
pbaspect([1.3 1 1])
%set(gca,'Fontsize',28);
%

%% --------------- Generate Fig. 3B ------------------------------------------------------------------------------------------ %%
init = default_init;
param = default_param;
param.H2O2_null_curve_switch = 0;
PRXSO2H_vector = [];
H2O2_vector = [];
Total_nonPRXSO2H_vector = [];
Flux_k1_vector = [];
Flux_k3_vector = [];
Flux_k5_vector = [];
init.PRXSO2H = 1;
increment = 1.01;

        while init.PRXSO2H <= 1000
            init.PRX = param.PRXtot - init.PRXSO2H;
            PRXSO2H_vector = [PRXSO2H_vector, init.PRXSO2H];
            y0 = cell2mat(struct2cell(init));
      
            [t,y]=ode15s('Bistability_ode', tspan, y0, [], param);
            H2O2_vector =  [H2O2_vector, y(end,1)];

            Total_nonPRXSO2H =  y(end, 2) + y(end, 3) + y(end, 4);
            Total_nonPRXSO2H_vector =  [Total_nonPRXSO2H_vector, Total_nonPRXSO2H];

            Flux_k1 = param.k1 * y(end, 1) * y(end, 2);
            Flux_k1_vector =  [Flux_k1_vector, Flux_k1];

            Flux_k3 = param.k3 * y(end, 1) * y(end, 3);
            Flux_k3_vector =  [Flux_k3_vector, Flux_k3];

            Flux_k5  = param.k5 * y(end, 1);
            Flux_k5_vector =  [Flux_k5_vector, Flux_k5];

    
            if  init.PRXSO2H == 99.99999 % 99.99999 is used to achieve low values for Flux_k1 and Total_nonPRXSO2H to plot on log scale;
                init.PRXSO2H = param.PRXtot * increment
            else
                init.PRXSO2H =  init.PRXSO2H * increment
            end

            if  init.PRXSO2H >= param.PRXtot & init.PRXSO2H < (param.PRXtot * increment)
                init.PRXSO2H = 99.99999
            end
            
%             %Visually check if steady state is reached
%             figure(3020)
%             plot(t, y(:,1))
%             xlabel ('Time')
%             ylabel ('H2O2')
%             hold on
        end

% find the first index where PRXSO2H >= 100 
index_100 = find(PRXSO2H_vector >= 100, 1); 
figure(302)
xline(100)

yyaxis left
loglog(PRXSO2H_vector(1 : index_100), Flux_k1_vector(1 : index_100), '-r', 'LineWidth', 2);
hold on
loglog(PRXSO2H_vector(index_100 : end), Flux_k1_vector(index_100 : end), ':r', 'LineWidth', 2);

loglog(PRXSO2H_vector(1 : index_100), Flux_k3_vector(1 : index_100), '-y', 'LineWidth', 2);
loglog(PRXSO2H_vector(index_100 : end), Flux_k3_vector(index_100 : end), ':y', 'LineWidth', 2);

loglog(PRXSO2H_vector(1 : index_100), Flux_k5_vector(1 : index_100), '-g', 'LineWidth', 2);
loglog(PRXSO2H_vector(index_100 : end), Flux_k5_vector(index_100 : end), ':g', 'LineWidth', 2);

xlim([10, 200])
ylim([2e-5, 100])
ylabel ('Flux (uM/S)')
set(gca,'YColor','k');

yyaxis right
loglog(PRXSO2H_vector(1 : index_100), H2O2_vector(1 : index_100), '-', 'Color', [51,102,255]/255, 'LineWidth', 2);
loglog(PRXSO2H_vector(index_100 : end), H2O2_vector(index_100 : end), ':', 'Color', [51,102,255]/255, 'LineWidth', 2);

loglog(PRXSO2H_vector(1 : index_100), Total_nonPRXSO2H_vector(1 : index_100), '-m', 'LineWidth', 2);
loglog(PRXSO2H_vector(index_100 : end), Total_nonPRXSO2H_vector(index_100 : end), ':m', 'LineWidth', 2);

ylabel ('Concentration (uM)')
set(gca,'YColor','k');
xlim([10, 200])
ylim([1e-3, 200])

xlabel ('PRXSO2H')
pbaspect([1.3 1 1])
set(gca,'YScale','log')
%

%% --------------- Generate Fig. 3C ------------------------------------------------------------------------------------------ %%
init = default_init;
param = default_param;
param.PRXSO2H_null_curve_switch = 0;
increment = 1.01;

for i = 1:1:2
    if i == 1
       param.k3_switch = 1;
    else 
       param.k3_switch = 0;
    end

    PRXSO2H_vector = [];
    H2O2_vector = [];
    init.H2O2 = 1e-6;    

    while init.H2O2 <= 1        
            init.PRX = param.PRXtot - init.PRXSO2H;    
            H2O2_vector = [H2O2_vector, init.H2O2];   
            y0 = cell2mat(struct2cell(init));   
            [t,y]=ode15s('Bistability_ode', tspan, y0, [], param);
            PRXSO2H_vector =  [PRXSO2H_vector, y(end,5)];
            init.H2O2 =  init.H2O2 * increment
    end

    figure(303)
    loglog(H2O2_vector, PRXSO2H_vector, 'Color', [255,50,50]/255, 'LineWidth', 2);
    xlabel ('H2O2')
    ylabel ('PRXSO2H')
    hold on
    xlim([1e-6, 1])
    ylim([0.01, 100])
    pbaspect([1.3 1 1])
    set(gca,'Fontsize',10);
end
%

%% --------------- Generate Fig. 3D ------------------------------------------------------------------------------------------ %%
init = default_init;
param = default_param;
param.PRXSO2H_null_curve_switch = 0; 
k4f_vector = [0.2*default_param.k4f, default_param.k4f, 5*default_param.k4f];
LRC_PRXSO2H = [];
increment = 1.01;

for i = 1 : 1 : 3
    param.k4f = k4f_vector(i);
    PRXSO2H_vector = [];
    H2O2_vector = [];
    init.H2O2 = 1e-6;
    
    while init.H2O2 <= 1        
            init.PRX = param.PRXtot - init.PRXSO2H;  
            H2O2_vector = [H2O2_vector, init.H2O2];    
            y0 = cell2mat(struct2cell(init));    
            [t,y]=ode15s('Bistability_ode', tspan, y0, [], param);    
            PRXSO2H_vector =  [PRXSO2H_vector, y(end,5)];   
            init.H2O2 =  init.H2O2 * increment
    end

    figure(304)
    loglog(H2O2_vector, PRXSO2H_vector, 'Color', [255,50,50]/255, 'LineWidth', 2);
    xlabel ('H2O2')
    ylabel ('PRXSO2H')
    hold on
    xlim([1e-3, 0.1])
    ylim([0.05, 200])
    pbaspect([1.3 1 1])
    set(gca,'Fontsize',10);
  
    % calculate and plot LRC for PRXSO2H null curve
    for j= 2:1:length(PRXSO2H_vector)
        if j==length(PRXSO2H_vector)
            break
        end
        delta_PRXSO2H=PRXSO2H_vector(j+1) - PRXSO2H_vector(j);
        PerInc_PRXSO2H=delta_PRXSO2H / PRXSO2H_vector(j);

        delta_H2O2=H2O2_vector(j+1) - H2O2_vector(j);
        PerInc_H2O2=delta_H2O2 / H2O2_vector(j);

        LRC_PRXSO2H(i, j-1) = PerInc_PRXSO2H / PerInc_H2O2;
    end

    [LRC_PRXSO2H_Max(i), I_H2O2(i)] = max(LRC_PRXSO2H(i,:));
    H2O2_value_for_LRC_PRXSO2H_Max(i) = H2O2_vector(1+I_H2O2(i));       
end

%

%% --------------- Generate Fig. 4A and 4B ------------------------------------------------------------------------------------------ %%
% Fig. 4A is reproduced from Fig. 3A and Fig. 4B is reproduced from Fig. 2A.

%% --------------- Generate Fig. 4C ------------------------------------------------------------------------------------------ %%
init = default_init;
param = default_param;
param.PRXSO2H_null_curve_switch = 0;
H2O2_vector = [];
Flux_k1_vector = [];
Flux_k3_vector = [];
Flux_k5_vector = [];
Flux_total_removal_vector = [];
init.H2O2 = 0.001;
increment = 1.01;

while init.H2O2 <= 100    
        init.PRX = param.PRXtot - init.PRXSO2H;
        H2O2_vector = [H2O2_vector, init.H2O2];
        y0 = cell2mat(struct2cell(init));
        [t,y]=ode15s('Bistability_ode', tspan, y0, [], param);

        Flux_k1 = param.k1 * y(end, 1) * y(end, 2);
        Flux_k1_vector =  [Flux_k1_vector, Flux_k1];

        Flux_k3 = param.k3 * y(end, 1) * y(end, 3);
        Flux_k3_vector =  [Flux_k3_vector, Flux_k3];

        Flux_k5  = param.k5 * y(end, 1);
        Flux_k5_vector =  [Flux_k5_vector, Flux_k5];

        Flux_total_removal = Flux_k1 + Flux_k3 + Flux_k5;
        Flux_total_removal_vector =  [Flux_total_removal_vector, Flux_total_removal];

        init.H2O2 =  init.H2O2 * increment 
end

figure(403)
loglog(H2O2_vector, Flux_k1_vector, '-r', 'LineWidth', 2);
hold on
loglog(H2O2_vector, Flux_k3_vector, '-y', 'LineWidth', 2);
loglog(H2O2_vector, Flux_k5_vector, '-g', 'LineWidth', 2);
loglog(H2O2_vector, Flux_total_removal_vector, '-b', 'LineWidth', 2);
xlabel('H2O2 (uM)')
ylabel('Flux (uM/S)')
pbaspect([1.05 1 1])
yline(44.68)
yline(27)
yline(15.58)
%

%% --------------- Generate Fig. S1 --------------------------------------------------------- %%
% Fig. S1 is reproduced from Fig. 4.
%

%% --------------- Generate Fig. 5A left panel and S2A - 2-parameter bifurcation k1 vs. k0 --------------------------------------------------------- %%
%Expect several hours to complete
param = default_param;
% H2O2, PRX, PRXSOH, PRXSS, PRXSO2H
init_off = [0, 100, 0, 0, 0]; % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100]; % initial concentrations for ON state;
num_steps = 400; 
k1_step_size = 300/num_steps;
k0_step_size = 150/num_steps;  

for i=1:1:num_steps + 1;
    i
    param.k1 = 1e-6 + (i-1) * k1_step_size;
    k1(i) = param.k1;
    
    for j = 1:1:num_steps + 1;
        j;
        param.k0 = 1e-6 + (j-1) * k0_step_size;
        k0(j) = param.k0;      

        [t,y] = ode15s('Bistability_ode', tspan, init_off, [], param);
        H2O2_switch_on(i,j) = y(end,1);

        [t,y] = ode15s('Bistability_ode', tspan, init_on, [], param);
        H2O2_switch_off(i,j) = y(end,1);
    end    
end

% save Two_par_k0k1_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0k1_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0k1_k1.mat k1
% save Two_par_k0k1_k0.mat k0

% Figure 5A left panel
figure(5011)
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;
Bistability_mag = log10(H2O2_ratio);
surf(k0, k1, Bistability_mag');  
pbaspect([1.1 1 1])
xlim([1e-6,150])
ylim([1e-6,300])
xlabel('k0')
ylabel('k1')
zlabel('Bistability magnitude')
shading interp 
colorbar
hold on

% Overlay XPP-AUT results to confirm
[num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 1);    
column_num = length(txt);
row_num = length(num);
figure(5011)
z_offset = zeros(row_num,1) + max(max(H2O2_ratio));
scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
xlabel(txt(1))
ylabel(txt(2))
box on

% Figure S2A
figure(5013) 
mesh(k0(1:20:end), k1(1:20:end), H2O2_switch_on(1:20:end,1:20:end));
hold on
mesh(k0(1:20:end), k1(1:20:end), H2O2_switch_off(1:20:end,1:20:end));
xlabel('k0')
ylabel('k1')
zlabel('H2O2 (uM)')
colorbar
%
%% --------------- Generate Fig. 5A right panel - Flux analysis ---------------------------------------------------------------- %%
init = default_init;
param = default_param;
param.PRXSO2H_null_curve_switch = 0;
increment = 1.1;

k1_vector = 2.^[-2:1:2] * default_param.k1;

for i = 1 : 1 : length(k1_vector)
    param.k1 = k1_vector(i);    
    H2O2_vector = [];
    Flux_k1_vector = [];
    Flux_k3_vector = [];
    Flux_k5_vector = [];
    Flux_total_removal_vector = [];

    init.H2O2 = 0.001;
    while init.H2O2 <= 100 

        init.PRX = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('Bistability_ode', tspan, y0, [], param);

        Flux_k1 = param.k1 * y(end, 1) * y(end, 2);
        Flux_k1_vector =  [Flux_k1_vector, Flux_k1];

        Flux_k3 = param.k3 * y(end, 1) * y(end, 3);
        Flux_k3_vector =  [Flux_k3_vector, Flux_k3];

        Flux_k5  = param.k5 * y(end, 1);
        Flux_k5_vector =  [Flux_k5_vector, Flux_k5];

        Flux_total_removal = Flux_k1 + Flux_k3 + Flux_k5;
        Flux_total_removal_vector =  [Flux_total_removal_vector, Flux_total_removal];

        init.H2O2 =  init.H2O2 * increment;   
    end

    % Figure 5A right panel
    figure(5012)
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  
    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);
    
    xlim([1e-3,5])
    ylim([1e-2,1e3])
    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')
    pbaspect([1.3 1 1])   
end
%

%% --------------- Generate Fig. 5B left panel and S2B - 2-parameter bifurcation k21 vs. k0  -------------------------------------------------------- %%
%Expect several hours to complete
param = default_param;
% H2O2, PRX, PRXSOH, PRXSS, PRXSO2H
init_off = [0, 100, 0, 0, 0]; % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100]; % initial concentrations for ON state;
num_steps = 400; 
k21_step_size = 3500/num_steps;
k0_step_size = 300/num_steps;  

for i=1:1:num_steps + 1;
    i
    param.k21 = 1e-6 + (i-1) * k21_step_size;
    k21(i) = param.k21;
    
    for j = 1:1:num_steps + 1;
        j;
        param.k0 = 1e-6 + (j-1) * k0_step_size;
        k0(j) = param.k0;      

        [t,y] = ode15s('Bistability_ode', tspan, init_off, [], param);
        H2O2_switch_on(i,j) = y(end,1);

        [t,y] = ode15s('Bistability_ode', tspan, init_on, [], param);
        H2O2_switch_off(i,j) = y(end,1);
    end    
end

% save Two_par_k0k21_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0k21_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0k21_k21.mat k21
% save Two_par_k0k21_k0.mat k0

% Figure 5B left panel
figure(5021)
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;
Bistability_mag = log10(H2O2_ratio);
surf(k0, k21, Bistability_mag');  
pbaspect([1.1 1 1])
xlim([1e-6,300])
ylim([1e-6,3500])
xlabel('k0')
ylabel('k21')
zlabel('Bistability magnitude')
shading interp 
colorbar
hold on

% Overlay XPP-AUT results to confirm
[num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 2); 
column_num = length(txt);
row_num = length(num);
figure(5021)
z_offset = zeros(row_num,1) + max(max(H2O2_ratio));
scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
xlabel(txt(1))
ylabel(txt(2))
box on

% Figure S2B
figure(5023) 
mesh(k0(1:20:end), k21(1:20:end), H2O2_switch_on(1:20:end,1:20:end));
hold on
mesh(k0(1:20:end), k21(1:20:end), H2O2_switch_off(1:20:end,1:20:end));
xlabel('k0')
ylabel('k21')
zlabel('H2O2 (uM)')
colorbar
%
%% --------------- Generate Fig. 5B right panel - Flux analysis --------------------------------------------------------------- %%
init = default_init;
param = default_param;
param.PRXSO2H_null_curve_switch = 0;

k21_vector = 2.^[-2:1:2] * default_param.k21; 

for i = 1 : 1 : length(k21_vector)   
    param.k21 = k21_vector(i);       
    H2O2_vector = [];
    Flux_k1_vector = [];
    Flux_k3_vector = [];
    Flux_k5_vector = [];
    Flux_total_removal_vector = [];

    init.H2O2 = 0.001;
    while init.H2O2 <= 100
        
        init.PRX = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('Bistability_ode', tspan, y0, [], param);

        Flux_k1 = param.k1 * y(end, 1) * y(end, 2);
        Flux_k1_vector =  [Flux_k1_vector, Flux_k1];

        Flux_k3 = param.k3 * y(end, 1) * y(end, 3);
        Flux_k3_vector =  [Flux_k3_vector, Flux_k3];

        Flux_k5  = param.k5 * y(end, 1);
        Flux_k5_vector =  [Flux_k5_vector, Flux_k5];

        Flux_total_removal = Flux_k1 + Flux_k3 + Flux_k5;
        Flux_total_removal_vector =  [Flux_total_removal_vector, Flux_total_removal];

        init.H2O2 =  init.H2O2 * 1.1;      
    end

    % Figure 5B right panel
    figure(5022) 
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  
    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);   
    xlim([1e-3,5])
    ylim([1e-2,1e3])
    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')
    pbaspect([1.3 1 1])   
end
%

%% --------------- Generate Fig. 5C left panel and S2C - 2-parameter bifurcation k22 vs. k0 -------------------------------------------------------- %%
%Expect several hours to complete
param = default_param;
% H2O2, PRX, PRXSOH, PRXSS, PRXSO2H
init_off = [0, 100, 0, 0, 0]; % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100]; % initial concentrations for ON state;

num_steps = 400; 
k22_step_size = 1.5/num_steps; 
k0_step_size = 150/num_steps;  

for i=1:1:num_steps + 1;
    i
    param.k22 = 1e-6 + (i-1) * k22_step_size;   
    k22(i) = param.k22;                         
    
    for j = 1:1:num_steps + 1;
        j;
        param.k0 = 1e-6 + (j-1) * k0_step_size;
        k0(j) = param.k0;      

        [t,y] = ode15s('Bistability_ode', tspan, init_off, [], param);
        H2O2_switch_on(i,j) = y(end,1);

        [t,y] = ode15s('Bistability_ode', tspan, init_on, [], param);
        H2O2_switch_off(i,j) = y(end,1);
    end    
end

% save Two_par_k0k22_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0k22_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0k22_k22.mat k22
% save Two_par_k0k22_k0.mat k0

% Figure 5C left panel
figure(5031)                                        
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;
Bistability_mag = log10(H2O2_ratio);
surf(k0, k22, Bistability_mag');                         
pbaspect([1.1 1 1])

xlim([1e-6,150])
ylim([1e-6,1.5])

xlabel('k0')
ylabel('k22')                                       
zlabel('Bistability magnitude')
shading interp 
colorbar
hold on

% Overlay XPP-AUT results to confirm
[num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 3);
column_num = length(txt);
row_num = length(num);
figure(5031)                                         
z_offset = zeros(row_num,1) + max(max(H2O2_ratio));
scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
xlabel(txt(1))
ylabel(txt(2))
box on

% Figure S2C
figure(5033) 
mesh(k0(1:20:end), k22(1:20:end), H2O2_switch_on(1:20:end,1:20:end));
hold on
mesh(k0(1:20:end), k22(1:20:end), H2O2_switch_off(1:20:end,1:20:end));
xlabel('k0')
ylabel('k22')                                       
zlabel('H2O2 (uM)')
colorbar
%
%% --------------- Generate Fig. 5C right panel - Flux analysis --------------------------------------------------------------- %%
init = default_init;
param = default_param;
param.PRXSO2H_null_curve_switch = 0;

k22_vector = 2.^[-2:1:2] * default_param.k22; 

for i = 1 : 1 : length(k22_vector)   
    param.k22 = k22_vector(i);      
    H2O2_vector = [];
    Flux_k1_vector = [];
    Flux_k3_vector = [];
    Flux_k5_vector = [];
    Flux_total_removal_vector = [];

    init.H2O2 = 0.001;
    while init.H2O2 <= 100
        
        init.PRX = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('Bistability_ode', tspan, y0, [], param);

        Flux_k1 = param.k1 * y(end, 1) * y(end, 2);
        Flux_k1_vector =  [Flux_k1_vector, Flux_k1];

        Flux_k3 = param.k3 * y(end, 1) * y(end, 3);
        Flux_k3_vector =  [Flux_k3_vector, Flux_k3];

        Flux_k5  = param.k5 * y(end, 1);
        Flux_k5_vector =  [Flux_k5_vector, Flux_k5];

        Flux_total_removal = Flux_k1 + Flux_k3 + Flux_k5;
        Flux_total_removal_vector =  [Flux_total_removal_vector, Flux_total_removal];

        init.H2O2 =  init.H2O2 * 1.1;    
    end
    
    % Figure 5C right panel
    figure(5032) 
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  
    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);    
    xlim([1e-3,5])
    ylim([1e-2,1e3])
    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')
    pbaspect([1.3 1 1])  
end
%

%% --------------- Generate Fig. 5D left panel and S2D - 2-parameter bifurcation k3 vs. k0 --------------------------------------------------------- %%
%Expect several hours to complete
param = default_param;
% H2O2, PRX, PRXSOH, PRXSS, PRXSO2H
init_off = [0, 100, 0, 0, 0]; % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100]; % initial concentrations for ON state;

num_steps = 400; 
k3_step_size = log10(10/1e-5)/num_steps; 
k0_step_size = 220/num_steps;

for i=1:1:num_steps + 1;
    i
    param.k3 = 10^(-5 + (i-1) * k3_step_size);  
    k3(i) = param.k3;                           
    
    for j = 1:1:num_steps + 1;
        j;
        param.k0 = 1e-6 + (j-1) * k0_step_size;
        k0(j) = param.k0;      

        [t,y] = ode15s('Bistability_ode', tspan, init_off, [], param);
        H2O2_switch_on(i,j) = y(end,1);

        [t,y] = ode15s('Bistability_ode', tspan, init_on, [], param);
        H2O2_switch_off(i,j) = y(end,1);
    end    
end

% save Two_par_k0k3_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0k3_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0k3_k3.mat k3
% save Two_par_k0k3_k0.mat k0

% Figure 5D left panel
figure(5041)                                        
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;
Bistability_mag = log10(H2O2_ratio);
surf(k0, k3, Bistability_mag');                          
pbaspect([1.1 1 1])
xlim([1e-6,220])
ylim([1e-6,10])
xlabel('k0')
ylabel('k3')                                       
zlabel('Bistability magnitude')
shading interp 
colorbar
hold on

% Overlay XPP-AUT results to confirm
[num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 4);    
column_num = length(txt);
row_num = length(num);
figure(5041)                                                   
z_offset = zeros(row_num,1) + max(max(Bistability_mag));
scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
xlabel(txt(1))
ylabel(txt(2))
box on

% Figure S2D
figure(5043)                            
mesh(k0(1:20:end), k3(1:20:end), H2O2_switch_on(1:20:end,1:20:end));
hold on
mesh(k0(1:20:end), k3(1:20:end), H2O2_switch_off(1:20:end,1:20:end));
xlabel('k0')
ylabel('k3')                                                    
zlabel('H2O2 (uM)')
colorbar
%
%% --------------- Generate Fig. 5D right panel - Flux analysis ---------------------------------------------------------------- %%
init = default_init;
param = default_param;
param.PRXSO2H_null_curve_switch = 0;

k3_vector = 2.^[-2:1:2] * default_param.k3;                     

for i = 1 : 1 : length(k3_vector)                               
    param.k3 = k3_vector(i);                                    
    H2O2_vector = [];
    Flux_k1_vector = [];
    Flux_k3_vector = [];
    Flux_k5_vector = [];
    Flux_total_removal_vector = [];

    init.H2O2 = 0.001;
    while init.H2O2 <= 100
        
        init.PRX = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('Bistability_ode', tspan, y0, [], param);

        Flux_k1 = param.k1 * y(end, 1) * y(end, 2);
        Flux_k1_vector =  [Flux_k1_vector, Flux_k1];

        Flux_k3 = param.k3 * y(end, 1) * y(end, 3);
        Flux_k3_vector =  [Flux_k3_vector, Flux_k3];

        Flux_k5  = param.k5 * y(end, 1);
        Flux_k5_vector =  [Flux_k5_vector, Flux_k5];

        Flux_total_removal = Flux_k1 + Flux_k3 + Flux_k5;
        Flux_total_removal_vector =  [Flux_total_removal_vector, Flux_total_removal];

        init.H2O2 =  init.H2O2 * 1.1;      
    end
    
    % Figure 5D right panel
    figure(5042)                                                
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  

    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);
       
    xlim([5e-3,2])
    ylim([1,2e2])
    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')
    pbaspect([1.3 1 1])    
end
%

%% --------------- Generate Fig. 5E and left panel and S5E - 2-parameter bifurcation k4f vs. k0 -------------------------------------------------------- %%
%Expect several hours to complete
param = default_param;
% H2O2, PRX, PRXSOH, PRXSS, PRXSO2H
init_off = [0, 100, 0, 0, 0];               % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100];               % initial concentrations for ON state;

num_steps = 400; 
k4f_step_size = 0.005/num_steps;             
k0_step_size = 150/num_steps;  

for i=1:1:num_steps + 1;
    i
    param.k4f = 1e-6 + (i-1) * k4f_step_size;  
    k4f(i) = param.k4f;                         
    
    for j = 1:1:num_steps + 1;
        j;
        param.k0 = 1e-6 + (j-1) * k0_step_size;
        k0(j) = param.k0;      
        
        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t,y] = ode15s('Bistability_ode', tspan, init_off, options, param);
        H2O2_switch_on(i,j) = y(end,1);

%         figure(100)
%         plot(t, y(:,1))
%         xlabel ('Time')
%         ylabel ('H2O2')
%         hold on
        
        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t,y] = ode15s('Bistability_ode', tspan, init_on, options, param);
        H2O2_switch_off(i,j) = y(end,1);

%         figure(100)
%         plot(t, y(:,1))
%         xlabel ('Time')
%         ylabel ('H2O2')
%         hold on
    end    
end

% save Two_par_k0k4f_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0k4f_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0k4f_k4f.mat k4f
% save Two_par_k0k4f_k0.mat k0

% Figure 5E left panel
figure(5051)                                         
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;
Bistability_mag = log10(H2O2_ratio);
surf(k0, k4f, Bistability_mag');   
pbaspect([1.1 1 1])
xlim([1e-6,150])
ylim([1e-6,0.005])
xlabel('k0')
ylabel('k4f')                                       
zlabel('Bistability magnitude')
shading interp 
colorbar
hold on

% Overlay XPP-AUT results to confirm
[num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 5);   
column_num = length(txt);
row_num = length(num);
figure(5051)                                                                                      
z_offset = zeros(row_num,1) + max(max(H2O2_ratio));
scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
xlabel(txt(1))
ylabel(txt(2))
box on

% Figure S2E
figure(5053)                             
mesh(k0(1:20:end), k4f(1:20:end), H2O2_switch_on(1:20:end,1:20:end));
hold on
mesh(k0(1:20:end), k4f(1:20:end), H2O2_switch_off(1:20:end,1:20:end));
xlabel('k0')
ylabel('k4f')                                                   
zlabel('H2O2 (uM)')
colorbar
%
%% --------------- Generate Fig. 5E right panel - Flux analysis --------------------------------------------------------------- %%
init = default_init;
param = default_param;
param.PRXSO2H_null_curve_switch = 0;

k4f_vector = 5.^[-2:1:2] * default_param.k4f;                     

for i = 1 : 1 : length(k4f_vector)                                
    param.k4f = k4f_vector(i);                                    
    H2O2_vector = [];
    Flux_k1_vector = [];
    Flux_k3_vector = [];
    Flux_k5_vector = [];
    Flux_total_removal_vector = [];

    init.H2O2 = 0.001;
    while init.H2O2 <= 100
        
        init.PRX = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('Bistability_ode', tspan, y0, [], param);

        Flux_k1 = param.k1 * y(end, 1) * y(end, 2);
        Flux_k1_vector =  [Flux_k1_vector, Flux_k1];

        Flux_k3 = param.k3 * y(end, 1) * y(end, 3);
        Flux_k3_vector =  [Flux_k3_vector, Flux_k3];

        Flux_k5  = param.k5 * y(end, 1);
        Flux_k5_vector =  [Flux_k5_vector, Flux_k5];

        Flux_total_removal = Flux_k1 + Flux_k3 + Flux_k5;
        Flux_total_removal_vector =  [Flux_total_removal_vector, Flux_total_removal];

        init.H2O2 =  init.H2O2 * 1.1;      
    end
    
    % Figure 5E right panel
    figure(5052)                                               
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  
    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);
      
    xlim([5e-3,2])
    ylim([1,2e2])
    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')
    pbaspect([1.3 1 1])
end
%

%% --------------- Generate Fig. 5F and left panel and S5F - 2-parameter bifurcation k4b vs. k0 -------------------------------------------------------- %%
%Expect several hours to complete
param = default_param;
% H2O2, PRX, PRXSOH, PRXSS, PRXSO2H
init_off = [0, 100, 0, 0, 0]; % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100]; % initial concentrations for ON state;

num_steps = 400; 
k4b_step_size = 0.1/num_steps;
k0_step_size = 150/num_steps;  

for i=1:1:num_steps + 1;
    i
    param.k4b = 1e-6 + (i-1) * k4b_step_size;     
    k4b(i) = param.k4b;                           
    
    for j = 1:1:num_steps + 1;
        j;
        param.k0 = 1e-6 + (j-1) * k0_step_size;
        k0(j) = param.k0;

        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t,y] = ode15s('Bistability_ode', tspan, init_off, options, param);
        H2O2_switch_on(i,j) = y(end,1);

        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t,y] = ode15s('Bistability_ode', tspan, init_on, options, param);
        H2O2_switch_off(i,j) = y(end,1);
    end    
end

% save Two_par_k0k4b_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0k4b_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0k4b_k4b.mat k4b
% save Two_par_k0k4b_k0.mat k0

% Figure 5F left panel
figure(5061)                                        
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;
Bistability_mag = log10(H2O2_ratio);     
surf(k0, k4b, Bistability_mag');       
pbaspect([1.1 1 1])
xlim([1e-6,150])
ylim([1e-6,0.1])
xlabel('k0')
ylabel('k4b')                                        
zlabel('Bistability magnitude')
shading interp 
colorbar
hold on

% Overlay XPP-AUT results to confirm
[num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 6);   
column_num = length(txt);
row_num = length(num);
figure(5061)                                                   
z_offset = zeros(row_num,1) + max(max(H2O2_ratio));
scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
xlabel(txt(1))
ylabel(txt(2))
box on

% Figure S2F
figure(5063)                              
mesh(k0(1:20:end), k4b(1:20:end), H2O2_switch_on(1:20:end,1:20:end));
hold on
mesh(k0(1:20:end), k4b(1:20:end), H2O2_switch_off(1:20:end,1:20:end));
xlabel('k0')
ylabel('k4b')                                                    
zlabel('H2O2 (uM)')
colorbar
%
%% --------------- Generate Fig. 5F right panel - Flux analysis --------------------------------------------------------------- %%
init = default_init;
param = default_param;
param.PRXSO2H_null_curve_switch = 0;

k4b_vector = 8.^[-2:1:2] * default_param.k4b;                    

for i = 1 : 1 : length(k4b_vector)                                
    param.k4b = k4b_vector(i);                                   
    H2O2_vector = [];
    Flux_k1_vector = [];
    Flux_k3_vector = [];
    Flux_k5_vector = [];
    Flux_total_removal_vector = [];

    init.H2O2 = 0.001;
    while init.H2O2 <= 100
        
        init.PRX = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('Bistability_ode', tspan, y0, [], param);

        Flux_k1 = param.k1 * y(end, 1) * y(end, 2);
        Flux_k1_vector =  [Flux_k1_vector, Flux_k1];

        Flux_k3 = param.k3 * y(end, 1) * y(end, 3);
        Flux_k3_vector =  [Flux_k3_vector, Flux_k3];

        Flux_k5  = param.k5 * y(end, 1);
        Flux_k5_vector =  [Flux_k5_vector, Flux_k5];

        Flux_total_removal = Flux_k1 + Flux_k3 + Flux_k5;
        Flux_total_removal_vector =  [Flux_total_removal_vector, Flux_total_removal];

        init.H2O2 =  init.H2O2 * 1.1;      
    end
    
    % Figure 5F right panel
    figure(5063)                                               
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  
    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);
    xlim([5e-3,2])
    ylim([1,2e2])
    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')
    pbaspect([1.3 1 1])  
end
%

%% --------------- Generate Fig. 5G and left panel and S5G - 2-parameter bifurcation k4c vs. k0 -------------------------------------------------------- %%
%Expect several hours to complete
param = default_param;
tspan1 = [0:100:1e7]; % need 1e7 to reach steady-state

% H2O2, PRX, PRXSOH, PRXSS, PRXSO2H
init_off = [0, 100, 0, 0, 0]; % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100]; % initial concentrations for ON state;

num_steps = 400; 
k4c_step_size = log10(10/1e-6)/num_steps;       
k0_step_size = 150/num_steps;  

for i=1:1:num_steps + 1;
    i
    param.k4c = 10^(-6 + (i-1) * k4c_step_size);   
    k4c(i) = param.k4c;                         
    
    for j = 1:1:num_steps + 1;
        j;
        param.k0 = 1e-5 + (j-1) * k0_step_size;
        k0(j) = param.k0;      

        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t,y] = ode15s('Bistability_ode', tspan1, init_off, options, param);
        H2O2_switch_on(i,j) = y(end,1);

        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t,y] = ode15s('Bistability_ode', tspan1, init_on, options, param);
        H2O2_switch_off(i,j) = y(end,1);
    end    
end

% save Two_par_k0k4c_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0k4c_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0k4c_k4c.mat k4c
% save Two_par_k0k4c_k0.mat k0

% Figure 5G left panel
figure(5071)                                        
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;
Bistability_mag = log10(H2O2_ratio);
surf(k0, k4c, Bistability_mag');                        
pbaspect([1.1 1 1])
xlabel('k0')
ylabel('k4c')                                        
zlabel('Bistability magnitude')
shading interp 
colorbar
hold on

% Overlay XPP-AUT results to confirm
[num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 7); 
row_num = length(num);
figure(5071)                                         
z_offset = zeros(row_num,1) + max(max(H2O2_ratio));
scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
xlabel(txt(1))
ylabel(txt(2))
box on
xlim([1e-6,150])
ylim([1e-6,10])

% Figure S2G
figure(5073) 
mesh(k0(1:20:end), k4c(1:20:end), H2O2_switch_on(1:20:end,1:20:end));
hold on
mesh(k0(1:20:end), k4c(1:20:end), H2O2_switch_off(1:20:end,1:20:end));
xlabel('k0')
ylabel('k4c')                                       
zlabel('H2O2 (uM)')
colorbar
%
%% --------------- Generate Fig. 5G right panel - Flux analysis --------------------------------------------------------------- %%
init = default_init;
param = default_param;
param.PRXSO2H_null_curve_switch = 0;

k4c_vector = 2.^[-2:1:2] * default_param.k4c; 

for i = 1 : 1 : length(k4c_vector)   
    param.k4c = k4c_vector(i);       
    H2O2_vector = [];
    Flux_k1_vector = [];
    Flux_k3_vector = [];
    Flux_k5_vector = [];
    Flux_total_removal_vector = [];

    init.H2O2 = 0.001;
    while init.H2O2 <= 100
        
        init.PRX = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('Bistability_ode', tspan, y0, [], param);

        Flux_k1 = param.k1 * y(end, 1) * y(end, 2);
        Flux_k1_vector =  [Flux_k1_vector, Flux_k1];

        Flux_k3 = param.k3 * y(end, 1) * y(end, 3);
        Flux_k3_vector =  [Flux_k3_vector, Flux_k3];

        Flux_k5  = param.k5 * y(end, 1);
        Flux_k5_vector =  [Flux_k5_vector, Flux_k5];

        Flux_total_removal = Flux_k1 + Flux_k3 + Flux_k5;
        Flux_total_removal_vector =  [Flux_total_removal_vector, Flux_total_removal];

        init.H2O2 =  init.H2O2 * 1.1;    
    end
    
    % Figure 5G right panel
    figure(5072) 
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  
    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);
    xlim([5e-3,2])
    ylim([1,2e2])
    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')
    pbaspect([1.3 1 1])   
end
%

%% --------------- Generate Fig. 5H and left panel and S5H - 2-parameter bifurcation SRXtot vs. k0 ----------------------------------------------------- %%
%Expect several hours to complete
param = default_param;
tspan1 = [0:100:1e7]; % need 1e7 to reach steady-state

% H2O2, PRX, PRXSOH, PRXSS, PRXSO2H
init_off = [0, 100, 0, 0, 0]; % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100]; % initial concentrations for ON state;

num_steps = 400; 
SRXtot_step_size = log10(100/1e-3)/num_steps; 
k0_step_size = 150/num_steps;  

for i=1:1:num_steps + 1;
    i
    param.SRXtot = 10^(-3 + (i-1) * SRXtot_step_size);        
    SRXtot(i) = param.SRXtot;                                 
    
    for j = 1:1:num_steps + 1;
        j;
        param.k0 = 1e-6 + (j-1) * k0_step_size; 
        k0(j) = param.k0;      

        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t,y] = ode15s('Bistability_ode', tspan1, init_off, options, param);
        H2O2_switch_on(i,j) = y(end,1);

        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t,y] = ode15s('Bistability_ode', tspan1, init_on, options, param);
        H2O2_switch_off(i,j) = y(end,1);
    end    
end

% save Two_par_k0SRXtot_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0SRXtot_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0SRXtot_SRXtot.mat SRXtot
% save Two_par_k0SRXtot_k0.mat k0

% Figure 5H left panel
figure(5081)                                                
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;   
Bistability_mag = log10(H2O2_ratio);
surf(k0, SRXtot, Bistability_mag');                         
pbaspect([1.1 1 1])
xlim([1e-6, 150])
ylim([1e-6, 100])
xlabel('k0')
ylabel('SRXtot')                                         
zlabel('Bistability magnitude')
shading interp 
colorbar
hold on

% Overlay XPP-AUT results to confirm
[num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 8); 
column_num = length(txt);
row_num = length(num);
figure(5081)                                                 
z_offset = zeros(row_num,1) + max(max(H2O2_ratio));
scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
xlabel(txt(1))
ylabel(txt(2))
box on

% Figure S2H
figure(5083) 
mesh(k0(1:20:end), SRXtot(1:20:end), H2O2_switch_on(1:20:end,1:20:end));
hold on
mesh(k0(1:20:end), SRXtot(1:20:end), H2O2_switch_off(1:20:end,1:20:end));
xlabel('k0')
ylabel('SRXtot')                                       
zlabel('H2O2 (uM)')
colorbar
%
%% --------------- Generate Fig. 5H right panel - Flux analysis ------------------------------------------------------------ %%
init = default_init;
param = default_param;
param.PRXSO2H_null_curve_switch = 0;

SRXtot_vector = 2.^[-2:1:2] * default_param.SRXtot;    

for i = 1 : 1 : length(SRXtot_vector)                    
    param.SRXtot = SRXtot_vector(i);                  
    H2O2_vector = [];
    Flux_k1_vector = [];
    Flux_k3_vector = [];
    Flux_k5_vector = [];
    Flux_total_removal_vector = [];

    init.H2O2 = 0.001;
    while init.H2O2 <= 100
        
        init.PRX = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('Bistability_ode', tspan, y0, [], param);

        Flux_k1 = param.k1 * y(end, 1) * y(end, 2);
        Flux_k1_vector =  [Flux_k1_vector, Flux_k1];

        Flux_k3 = param.k3 * y(end, 1) * y(end, 3);
        Flux_k3_vector =  [Flux_k3_vector, Flux_k3];

        Flux_k5  = param.k5 * y(end, 1);
        Flux_k5_vector =  [Flux_k5_vector, Flux_k5];

        Flux_total_removal = Flux_k1 + Flux_k3 + Flux_k5;
        Flux_total_removal_vector =  [Flux_total_removal_vector, Flux_total_removal];

        init.H2O2 =  init.H2O2 * 1.1;     
    end
    
    % Figure 5H right panel
    figure(5102) 
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  
    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);
       
    xlim([5e-3,2])
    ylim([1,2e2])
    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')
    pbaspect([1.3 1 1])
end
%

%% --------------- Generate Fig. 5I and left panel and S5I - 2-parameter bifurcation k5 vs. k0 --------------------------------------------------------- %%
%Expect several hours to complete
param = default_param;
% H2O2, PRX, PRXSOH, PRXSS, PRXSO2H
init_off = [0, 100, 0, 0, 0]; % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100]; % initial concentrations for ON state;

num_steps = 400; 
k5_step_size = 1000/num_steps;
k0_step_size = 150/num_steps;  

for i=1:1:num_steps + 1;
    i
    param.k5 = 1e-6 + (i-1) * k5_step_size;    
    k5(i) = param.k5;                          
    
    for j = 1:1:num_steps + 1;
        j;
        param.k0 = 1e-6 + (j-1) * k0_step_size;
        k0(j) = param.k0;      

        [t,y] = ode15s('Bistability_ode', tspan, init_off, [], param);
        H2O2_switch_on(i,j) = y(end,1);

        [t,y] = ode15s('Bistability_ode', tspan, init_on, [], param);
        H2O2_switch_off(i,j) = y(end,1);
    end    
end

% save Two_par_k0k5_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0k5_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0k5_k5.mat k5
% save Two_par_k0k5_k0.mat k0

% Figure 5I left panel
figure(5091)                                       
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;
Bistability_mag = log10(H2O2_ratio);
surf(k0, k5, Bistability_mag');  
pbaspect([1.1 1 1])
xlim([1e-6,150])
ylim([1e-6,1000])
xlabel('k0')
ylabel('k5')                                      
zlabel('Bistability magnitude')
shading interp 
colorbar
hold on

% Overlay XPP-AUT results to confirm
[num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 9); 
column_num = length(txt);
row_num = length(num);
figure(5091)                                        
z_offset = zeros(row_num,1) + max(max(H2O2_ratio));
scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
xlabel(txt(1))
ylabel(txt(2))
box on

% Figure S2I
figure(5093) 
mesh(k0(1:10:end), k5(1:1:10), H2O2_switch_on(1:10:end,1:1:10)');
hold on
mesh(k0(1:10:end), k5(1:1:10), H2O2_switch_off(1:10:end,1:1:10)');
xlabel('k0')
ylabel('k5')                                       
zlabel('H2O2 (uM)')
colorbar
%
%% --------------- Generate Fig. 5I right panel - Flux analysis ---------------------------------------------------------------- %%
init = default_init;
param = default_param;
param.PRXSO2H_null_curve_switch = 0;

k5_vector = 2.^[-2:1:2] * default_param.k5; 

for i = 1 : 1 : length(k5_vector)   
    param.k5 = k5_vector(i);       
    H2O2_vector = [];
    Flux_k1_vector = [];
    Flux_k3_vector = [];
    Flux_k5_vector = [];
    Flux_total_removal_vector = [];

    init.H2O2 = 0.001;
    while init.H2O2 <= 100
        
        init.PRX = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('Bistability_ode', tspan, y0, [], param);

        Flux_k1 = param.k1 * y(end, 1) * y(end, 2);
        Flux_k1_vector =  [Flux_k1_vector, Flux_k1];

        Flux_k3 = param.k3 * y(end, 1) * y(end, 3);
        Flux_k3_vector =  [Flux_k3_vector, Flux_k3];

        Flux_k5  = param.k5 * y(end, 1);
        Flux_k5_vector =  [Flux_k5_vector, Flux_k5];

        Flux_total_removal = Flux_k1 + Flux_k3 + Flux_k5;
        Flux_total_removal_vector =  [Flux_total_removal_vector, Flux_total_removal];

        init.H2O2 =  init.H2O2 * 1.1;    
    end
    
    % Figure 5I right panel
    figure(5082) 
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  
    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);
       
    xlim([5e-3,2])
    ylim([1,2e2])
    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')
    pbaspect([1.3 1 1])
end
%

%% --------------- Generate Fig. 5J and left panel and S5J - 2-parameter bifurcation PRXtot vs. k0 ----------------------------------------------------- %%
%Expect several hours to complete
param = default_param;

num_steps = 400; 
PRXtot_step_size = 1000/num_steps;
k0_step_size = 200/num_steps;  

for i=1:1:num_steps + 1;
    i
    param.PRXtot = 1e-6 + (i-1) * PRXtot_step_size; 
    PRXtot(i) = param.PRXtot;                             
    
    % H2O2, PRX, PRXSOH, PRXSS, PRXSO2H
    init_off = [0, PRXtot(i), 0, 0, 0]; % initial concentrations for OFF state;
    init_on = [10, 0, 0, 0, PRXtot(i)]; % initial concentrations for ON state;

    for j = 1:1:num_steps + 1;
        j;
        param.k0 = 1e-6 + (j-1) * k0_step_size;   
        k0(j) = param.k0;      

        [t,y] = ode15s('Bistability_ode', tspan, init_off, [], param);
        H2O2_switch_on(i,j) = y(end,1);

        [t,y] = ode15s('Bistability_ode', tspan, init_on, [], param);
        H2O2_switch_off(i,j) = y(end,1);
    end    
end

% save Two_par_k0PRXtot_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0PRXtot_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0PRXtot_PRXtot.mat PRXtot
% save Two_par_k0PRXtot_k0.mat k0

% Figure 5J left panel
figure(5101)                                              
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;
Bistability_mag = log10(H2O2_ratio);
surf(k0, PRXtot, Bistability_mag');                        
pbaspect([1.1 1 1])
xlim([1e-6,200])
ylim([1e-6,1000])
xlabel('k0')
ylabel('PRXtot')                                      
zlabel('Bistability magnitude')
shading interp 
colorbar
hold on

% Overlay XPP-AUT results to confirm
[num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 10); 
column_num = length(txt);
row_num = length(num);
figure(5101)                                               
z_offset = zeros(row_num,1) + max(max(H2O2_ratio));
scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
xlabel(txt(1))
ylabel(txt(2))
box on
 
% Figure S2J
figure(5103) 
mesh(k0(1:20:end), PRXtot(1:20:end), H2O2_switch_on(1:20:end,1:20:end));
hold on
mesh(k0(1:20:end), PRXtot(1:20:end), H2O2_switch_off(1:20:end,1:20:end));
xlabel('k0')
ylabel('PRXtot')                                         
zlabel('H2O2 (uM)')
colorbar
%
%% --------------- Generate Fig. 5J right panel - Flux analysis ------------------------------------------------------------ %%
init = default_init;
param = default_param;
param.PRXSO2H_null_curve_switch = 0;

PRXtot_vector = 2.^[-2:1:2] * default_param.PRXtot;   

for i = 1 : 1 : length(PRXtot_vector)                   
    param.PRXtot = PRXtot_vector(i);                 
    H2O2_vector = [];
    Flux_k1_vector = [];
    Flux_k3_vector = [];
    Flux_k5_vector = [];
    Flux_total_removal_vector = [];

    init.H2O2 = 0.001;
    while init.H2O2 <= 100
        
        init.PRX = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('Bistability_ode', tspan, y0, [], param);

        Flux_k1 = param.k1 * y(end, 1) * y(end, 2);
        Flux_k1_vector =  [Flux_k1_vector, Flux_k1];

        Flux_k3 = param.k3 * y(end, 1) * y(end, 3);
        Flux_k3_vector =  [Flux_k3_vector, Flux_k3];

        Flux_k5  = param.k5 * y(end, 1);
        Flux_k5_vector =  [Flux_k5_vector, Flux_k5];

        Flux_total_removal = Flux_k1 + Flux_k3 + Flux_k5;
        Flux_total_removal_vector =  [Flux_total_removal_vector, Flux_total_removal];

        init.H2O2 =  init.H2O2 * 1.1;      
    end
    
    % Figure 5J right panel
    figure(5092) 
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  
    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);
       
    xlim([5e-3,2])
    ylim([1,2e2])
    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')
    pbaspect([1.3 1 1])
end
%

%% --------------- Generate Fig. S3 --------------------------------------------------------- %%
% Fig. S3 is reproduced from Fig. 5.
%

%% --------------- Generate Fig. 6 - Time course for varying k0 ---------------------------------------------------------------- %%
param = default_param;
tspan = [0:100:3600*60]; 
% H2O2, PRX, PRXSOH, PRXSS, PRXSO2H
init_off = [0, 100, 0, 0, 0];               % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100];               % initial concentrations for ON state;

% Starting from initial OFF state
param.k0 = 10
[t,y0] = ode15s('Bistability_ode', tspan, init_off, [], param); 
init_off = y0(end,:);  

num_steps = 10;         
k0_step_size = 100/num_steps; 
k0 = [];
for j = 1:1:num_steps;
    j;
    param.k0 = 10 + (j-1) * k0_step_size;   
    k0(j) = param.k0; 

    options = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [t,y] = ode15s('Bistability_ode', tspan, init_off, options, param);
    H2O2_switch_on(j) = y(end,1);

    if  j <= 4  
        line_color =  [0.47,0.67,0.19];
    else
        line_color =  [0.85,0.33,0.10];
    end

    figure(601)
    semilogy(t/3600, y(:,1), 'Color', line_color, 'LineWidth', 2)
    xlabel ('Time (h)')
    ylabel ('H2O2')
    ylim([1e-3,10])
    pbaspect([1.2 1 1])
    hold on  

    figure(602)
    semilogy(t/3600, y(:,2), 'Color', line_color, 'LineWidth', 2)
    xlabel ('Time (h)')
    ylabel ('PRX')
    ylim([1e-3,1e3])
    pbaspect([1.2 1 1])
    hold on

    figure(603)
    semilogy(t/3600, y(:,3) + y(:,4), 'Color', line_color, 'LineWidth', 2)
    xlabel ('Time (h)')
    ylabel ('PRXSOH + PRXSS')
    ylim([1e-1,1e2])
    pbaspect([1.2 1 1])
    hold on

    figure(604)
    semilogy(t/3600, (param.PRXtot-y(:,2)-y(:,3)-y(:,4)), 'Color', line_color, 'LineWidth', 2)
    xlabel ('Time (h)')
    ylabel ('PRXSO2H + PS')
    ylim([1e-2,1e3])
    pbaspect([1.2 1 1])
    hold on
end   

% Starting from initial ON state
param.k0 = 100
[t,y0] = ode15s('Bistability_ode', tspan, init_on, [], param);
init_on = y0(end,:);  

for i = 1:1:num_steps;
    i;
    param.k0 = 5 + (i-1) * k0_step_size;   % 5
    k0(i) = param.k0; 

    options = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [t,y] = ode15s('Bistability_ode', tspan, init_on, options, param);
    H2O2_switch_off(i) = y(end,1);

    if  i <= 2
        line_color =  [0.47,0.67,0.19];
    else
        line_color =  [0.85,0.33,0.10];
    end

    figure(605)
    semilogy(t/3600, y(:,1), 'Color', line_color, 'LineWidth', 2)
    xlabel ('Time (h)')
    ylabel ('H2O2')
    ylim([1e-3,10])
    pbaspect([1.2 1 1])
    hold on  

    figure(606)
    semilogy(t/3600, y(:,2), 'Color', line_color, 'LineWidth', 2)
    xlabel ('Time (h)')
    ylabel ('PRX')
    ylim([1e-3,1e3])
    pbaspect([1.2 1 1])
    hold on

    figure(607)
    semilogy(t/3600, y(:,3) + y(:,4), 'Color', line_color, 'LineWidth', 2)
    xlabel ('Time (h)')
    ylabel ('PRXSOH + PRXSS')
    ylim([1e-1,1e2])
    pbaspect([1.2 1 1])
    hold on

    figure(608)
    semilogy(t/3600, (param.PRXtot-y(:,2)-y(:,3)-y(:,4)), 'Color', line_color, 'LineWidth', 2)
    xlabel ('Time (h)')
    ylabel ('PRXSO2H + PS')
    ylim([1e-2,1e3])
    pbaspect([1.2 1 1])
    hold on
end   
%

%% --------------- Generate Fig. 7A and 7B - Saddle-node bifurcation from XPP-AUT results ----------------------------------------------------------------------------------- %%
% Figure 7A
% k0 = 2.5, 5, 10, 20, 40, 80, 160 (number of sheets = 7)
for i = 1 : 1 : 7 
    [num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_SRXtot_versus_H2O2', i);
    column_num = length(txt);
    
    % Obtain branch index: 
    branch_1_index = find(num(:, column_num) == 1); 
    branch_2_index = find(num(:, column_num) == 2);
    branch_3_index = find(num(:, column_num) == 3);
    
    % Obtain branch data: 
    branch_1 = num(branch_1_index, [1, 2]);
    branch_2 = num(branch_2_index, [1, 2]);
    branch_3 = num(branch_3_index, [1, 2]);
       
    figure(701)
    loglog(branch_1(:, 1), branch_1(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    hold on
    loglog(branch_2(:, 1), branch_2(:, 2), ':', 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_3(:, 1), branch_3(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    pbaspect([1.2 0.9 1])
    xlabel(txt(1))
    ylabel(txt(2))
    box on
end

% Figure 7B
% k0 = 2.5, 5, 10, 20, 40, 80, 160 (number of sheets = 7)
for i = 1 : 1 : 7 
    [num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_SRXtot_versus_PRXSO2H.xlsx', i);
    column_num = length(txt);
    
    % Obtain branch index: 
    branch_1_index = find(num(:, column_num) == 1); 
    branch_2_index = find(num(:, column_num) == 2);
    branch_3_index = find(num(:, column_num) == 3);
    
    % Obtain branch data: 
    branch_1 = num(branch_1_index, [1, 2]);
    branch_2 = num(branch_2_index, [1, 2]);
    branch_3 = num(branch_3_index, [1, 2]);
    
    figure(702)
    loglog(branch_1(:, 1), branch_1(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    hold on
    loglog(branch_2(:, 1), branch_2(:, 2), ':', 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_3(:, 1), branch_3(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    pbaspect([1.2 0.9 1])
    xlabel(txt(1))
    ylabel(txt(2))
    box on
end
%

%% --------------- Generate Fig. 7C and 7D - Time course for varying SRXtot under different k0   ------------------------------- %%
param = default_param;
param.k0 = 10; % Manually change k0 to different values: k0 = 10, 20, 40, 80
tspan = [0:100:3600*400]; 

% H2O2, PRX, PRXSOH, PRXSS, PRXSO2H
init_off = [0, 100, 0, 0, 0];               % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100];               % initial concentrations for ON state;

% obtain initial OFF steady state values
param.SRXtot = 10;
[t,y0] = ode15s('Bistability_ode', tspan, init_off, [], param);
init_off = y0(end,:);  

% obtain initial ON steady state values
param.SRXtot = 0.001;
[t,y0] = ode15s('Bistability_ode', tspan, init_on, [], param);
init_on = y0(end,:);  

SRXtot_max_on = 3;
SRXtot_max_off = 30;
SRXtot_min = 1e-4;
num_steps = 10;     

for j = 0:1:num_steps;
    j;

    SRXtot_step_size_on = SRXtot_max_on/num_steps; 
    if  j == 0
        param.SRXtot = SRXtot_min;
    else
        param.SRXtot = j * SRXtot_step_size_on;
    end
    SRXtot_on(j+1) = param.SRXtot;  

    [t,y] = ode15s('Bistability_ode', tspan, init_off, [], param);
    H2O2_switch_on(j+1) = y(end,1);

    if  y(end,1) > 0.1
        line_color =  [0.85,0.33,0.10];
    else
        line_color =  [0.47,0.67,0.19];
    end

    figure(703)
    semilogy(t/3600, y(:,1), 'Color', line_color, 'LineWidth', 2);
    xlabel ('Time (h)')
    ylabel ('H2O2')
    set(gca, 'xtick', 0:100:400)
    set(gca, 'ytick', [0.001,0.01,0.1,1,10])
    xlim([0,400])
    ylim([1e-3,10])
    set(gca,'fontsize',28);
    set(get(gca,'XLabel'),'FontSize',28);
    set(get(gca,'YLabel'),'FontSize',28);
    pbaspect([1.2 1 1])
    hold on  

%     figure(7031)
%     semilogy(t/3600, y(:,2), 'Color', line_color, 'LineWidth', 2);
%     xlabel ('Time (h)')
%     ylabel ('PRX')
%     ylim([1e-3,1e3])
%     hold on
% 
%     figure(7032)
%     semilogy(t/3600, y(:,3) + y(:,4), 'Color', line_color, 'LineWidth', 2);
%     xlabel ('Time (h)')
%     ylabel ('PRXSOH + PRXSS')
%     ylim([1e-1,1e2])
%     hold on
% 
%     figure(7033)
%     semilogy(t/3600, (param.PRXtot-y(:,2)-y(:,3)-y(:,4)), 'Color', line_color, 'LineWidth', 2);
%     xlabel ('Time (h)')
%     ylabel ('PRXSO2H + PS')
%     ylim([1e-2,1e3])
%     hold on 

    SRXtot_step_size_off = SRXtot_max_off/num_steps;  
    if  j == 0
        param.SRXtot = SRXtot_min;
    else
        param.SRXtot = j * SRXtot_step_size_off;
    end
    SRXtot_off(j+1) = param.SRXtot;    

    [t,y] = ode15s('Bistability_ode', tspan, init_on, [], param);
    H2O2_switch_off(j+1) = y(end,1);

    if  y(end,1) > 0.1
        line_color =  [0.85,0.33,0.10];
    else
        line_color =  [0.47,0.67,0.19];
    end

    figure(704)
    semilogy(t/3600, y(:,1), 'Color', line_color, 'LineWidth', 2);
    xlabel ('Time (h)')
    ylabel ('H2O2')
    set(gca, 'xtick', 0:2:10)
    set(gca, 'ytick', [0.001,0.01,0.1,1,10])
    xlim([0,10])
    ylim([1e-3,10])
    set(gca,'fontsize',28);
    set(get(gca,'XLabel'),'FontSize',28);
    set(get(gca,'YLabel'),'FontSize',28);
    pbaspect([1.2 1 1])
    hold on  

%     figure(7041)
%     semilogy(t/3600, y(:,2), 'Color', line_color, 'LineWidth', 2);
%     xlabel ('Time (h)')
%     ylabel ('PRX')
%     xlim([0,10])
%     ylim([1e-3,1e3])
%     hold on
% 
%     figure(7042)
%     semilogy(t/3600, y(:,3) + y(:,4), 'Color', line_color, 'LineWidth', 2);
%     xlabel ('Time (h)')
%     ylabel ('PRXSOH + PRXSS')
%     xlim([0,10])
%     ylim([1e-1,1e2])
%     hold on
% 
%     figure(7043)
%     semilogy(t/3600, (param.PRXtot-y(:,2)-y(:,3)-y(:,4)), 'Color', line_color, 'LineWidth', 2);
%     xlabel ('Time (h)')
%     ylabel ('PRXSO2H + PS')
%     xlim([0,10])
%     ylim([1e-2,1e3])
%     hold on
end
