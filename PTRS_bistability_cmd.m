%%  -------------- Matlab codes of Bistability Model -------------------------------------------------------------- %%
% Concentration unit: uM, time unit: second
close all
clear all
clc
tspan = [0:100:1000000];

%

%% --------------- PARAMETERS ------------------------------------------------------------------------------------- %%
param.k0 = 10;
param.k1 = 50;
param.k2a = 10;
param.k2b = 1;
param.k2c = 60;         
param.Km2c = 3;    
param.k3 = 0.002;
param.k4f = 0.0014;
param.k4b = 0.001;
param.k4c = 0.006;     
param.k5 = 110;         
param.PRXtot = 100;
param.SRXtot = 0.3;     
param.TRXtot = 30;      
param.TRtot = 5.52;    
param.PRXSO2H_switch = 1;
param.TRXSS_switch = 1; 
param.PRXSS_switch = 1; 
param.PS_switch = 1; 
param.PRXSH_switch = 1;
param.H2O2_switch = 1; 
param.k3_switch = 1;
default_param = param;

%

%% --------------- INITIAL CONDITION ------------------------------------------------------------------------------ %%
init.H2O2 = 0.01;
init.PRXSH = param.PRXtot;
init.PRXSOH = 0;
init.PRXSS = 0;
init.PRXSO2H = 0; 
init.TRXSS = 0;
default_init = init;

%

%% --------------- Generate final_fig. 5A-5I - Saddle-node bifurcation from XPP-AUT results ----------------------- %%

[num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_k0.xlsx'); % Bifurcation_k0.xlsx contains bufucation results obtained from XPP-AUT
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
    xlim([10, 1000])
    xlabel(txt(1))
    ylabel(txt(1+i))
    box on
    pbaspect([1.1 1 1])
    set(gca,'Fontsize', 16);

end

%

%% --------------- Generate final_fig. 6A - H2O2 null curves ------------------------------------------------------ %%
init = default_init;
param = default_param;
param.PRXSO2H_switch = 0;
k0_vector = [61.31, 100, 152.3];
LRC_H2O2 = [];
increment = 1.01;

for i = 1 : 1 : 3
        param.k0 = k0_vector(i);
        PRXSO2H_vector = [];
        H2O2_vector = [];
        init.PRXSO2H = 0.1;
        
                while init.PRXSO2H <= 1000
                    init.PRXSH = param.PRXtot - init.PRXSO2H;
                    PRXSO2H_vector = [PRXSO2H_vector, init.PRXSO2H];
                    y0 = cell2mat(struct2cell(init));
                    [t,y]=ode15s('PTRS_ode', tspan, y0, [], param);
                    H2O2_vector = [H2O2_vector, y(end,1)];
                    init.PRXSO2H = init.PRXSO2H * increment;

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
        loglog(PRXSO2H_vector(1 : index_100), H2O2_vector(1 : index_100),  'Color', [51,102,255]/255, 'LineWidth', 3);
        hold on  
        loglog(PRXSO2H_vector(index_100 : end), H2O2_vector(index_100 : end), ':', 'Color', [51,102,255]/255, 'LineWidth', 3);        
        xlabel ('PRXSO2H')
        ylabel ('H2O2')      
        pbaspect([1.3 1 1])
        set(gca,'Fontsize',22);
                          
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
        lineH2O2 = plot(PRXSO2H_vector(2:(end-1)), LRC_H2O2(i,:), 'LineWidth', 3);
        xlabel('PRXSO2H');
        ylabel('LRC of H2O2');
        hold on
        box on
        pbaspect([1.3 1 1])
        set(gca,'Fontsize',22);
end
%

%% --------------- Generate final_fig. 6A - PRXSO2H null curve ---------------------------------------------------- %%
init = default_init;
param = default_param;
param.H2O2_switch = 0;
PRXSO2H_vector = [];
H2O2_vector = [];
init.H2O2 = 0.001;
LRC_PRXSO2H = [];
increment = 1.01;

while init.H2O2 <= 100    
        init.PRXSH = param.PRXtot - init.PRXSO2H;
        H2O2_vector = [H2O2_vector, init.H2O2];
        y0 = cell2mat(struct2cell(init));
        [t,y]=ode15s('PTRS_ode', tspan, y0, [], param);
        PRXSO2H_vector =  [PRXSO2H_vector, y(end,5)];
        init.H2O2 =  init.H2O2 * increment;

%         #Visually check if steady state is reached
%         figure(3012)
%         plot(t, y(:,5))
%         xlabel ('Time')
%         ylabel ('PRXSO2H')
%         hold on   
end

figure(301)
loglog(PRXSO2H_vector, H2O2_vector, 'Color', [255,50,50]/255, 'LineWidth', 3);
xlabel ('PRXSO2H')
ylabel ('H2O2')
hold on
xlim([0.1, 300])
ylim([0.001, 5])
pbaspect([1.3 1 1])
set(gca,'Fontsize', 22);

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
linePRXSO2H = semilogx(H2O2_vector(2:(end-1)), LRC_PRXSO2H, 'LineWidth', 3);
xlabel('H2O2');
ylabel('LRC of PRXSO2H');
hold on
box on
pbaspect([1.3 1 1])
set(gca,'Fontsize', 22);
%

%% --------------- Generate final_fig. 6B ------------------------------------------------------------------------- %%
init = default_init;
param = default_param;
param.PRXSO2H_switch = 0;
PRXSO2H_vector = [];
H2O2_vector = [];
Total_nonPRXSO2H_vector = [];
Flux_k1_vector = [];
Flux_k3_vector = [];
Flux_k5_vector = [];
init.PRXSO2H = 1;
increment = 1.01;

        while init.PRXSO2H <= 1000
            init.PRXSH = param.PRXtot - init.PRXSO2H;
            PRXSO2H_vector = [PRXSO2H_vector, init.PRXSO2H];
            y0 = cell2mat(struct2cell(init));
      
            [t,y]=ode15s('PTRS_ode', tspan, y0, [], param);
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

%% --------------- Generate final_fig. 6C ------------------------------------------------------------------------- %%
init = default_init;
param = default_param;
param.H2O2_switch = 0;
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
            init.PRXSH = param.PRXtot - init.PRXSO2H;    
            H2O2_vector = [H2O2_vector, init.H2O2];   
            y0 = cell2mat(struct2cell(init));   
            [t,y]=ode15s('PTRS_ode', tspan, y0, [], param);
            PRXSO2H_vector =  [PRXSO2H_vector, y(end,5)];
            init.H2O2 =  init.H2O2 * increment;
    end

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

    figure(303)
    yyaxis left
    set(gca,'YColor','k');
    loglog(H2O2_vector, PRXSO2H_vector, '-', 'Color', [255,50,50]/255, 'LineWidth', 3, 'DisplayName', 'PRXSO2H');
    hold on
    xlabel ('H2O2 (uM)')
    ylabel ('Concentration (uM)')
    xlim([1e-6, 1])
    xticks([1e-3, 1e-2, 1e-1, 1, 10])
    ylim([0.01, 100])
    yticks([0.0001, 0.01, 1, 100])
    pbaspect([1.3 1 1])
    set(gca,'Fontsize', 22);
    legend
    
    yyaxis right
    semilogx(H2O2_vector(1:(end-2)), LRC_PRXSO2H, ':', 'Color', [0,0.45,0.74], 'LineWidth', 2, 'DisplayName', 'LRC PRXSO2H');
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










end
%

%% --------------- Generate final_fig. 6D ------------------------------------------------------------------------- %%
init = default_init;
param = default_param;
param.H2O2_switch = 0; 
k4f_vector = [0.2*default_param.k4f, default_param.k4f, 5*default_param.k4f];
LRC_PRXSO2H = [];
increment = 1.01;

for i = 1 : 1 : 3
    param.k4f = k4f_vector(i);
    PRXSO2H_vector = [];
    H2O2_vector = [];
    init.H2O2 = 1e-6;
    
    while init.H2O2 <= 1        
            init.PRXSH = param.PRXtot - init.PRXSO2H;  
            H2O2_vector = [H2O2_vector, init.H2O2];    
            y0 = cell2mat(struct2cell(init));    
            [t,y]=ode15s('PTRS_ode', tspan, y0, [], param);    
            PRXSO2H_vector =  [PRXSO2H_vector, y(end,5)];   
            init.H2O2 =  init.H2O2 * increment;
    end
 
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


    figure(304)
    yyaxis left
    set(gca,'YColor','k');
    loglog(H2O2_vector, PRXSO2H_vector, '-', 'Color', [255,50,50]/255, 'LineWidth', 3, 'DisplayName', 'PRXSO2H');
    hold on
    xlabel ('H2O2 (uM)')
    ylabel ('Concentration (uM)')
    xlim([1e-3, 0.1])
    xticks([1e-3, 1e-2, 1e-1, 1, 10])
    ylim([0.05, 200])
    yticks([0.0001, 0.01, 1, 100])
    pbaspect([1.3 1 1])
    set(gca,'Fontsize', 22);
    legend
    
    yyaxis right
    semilogx(H2O2_vector(1:(end-2)), LRC_PRXSO2H, ':', 'Color', [0,0.45,0.74], 'LineWidth', 2, 'DisplayName', 'LRC PRXSO2H');
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

end

%

%% --------------- Generate final_fig. S3A and S3B ------------------------------------------------------------------ %%
% Fig. S3A is reproduced from Fig. 6A and Fig. S3B is reproduced from Fig. 5A.

%% --------------- Generate final_fig. S3C ------------------------------------------------------------------------- %%
init = default_init;
param = default_param;
param.H2O2_switch = 0;
H2O2_vector = [];
Flux_k1_vector = [];
Flux_k3_vector = [];
Flux_k5_vector = [];
Flux_total_removal_vector = [];
init.H2O2 = 0.001;
increment = 1.01;

while init.H2O2 <= 100  

        init.PRXSH = param.PRXtot - init.PRXSO2H;
        H2O2_vector = [H2O2_vector, init.H2O2];
        y0 = cell2mat(struct2cell(init));
        [t,y]=ode15s('PTRS_ode', tspan, y0, [], param);

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

figure(403)

loglog(H2O2_vector, Flux_k1_vector, '-r', 'LineWidth', 2);
hold on
loglog(H2O2_vector, Flux_k3_vector, '-y', 'LineWidth', 2);
loglog(H2O2_vector, Flux_k5_vector, '-g', 'LineWidth', 2);
loglog(H2O2_vector, Flux_total_removal_vector, '-b', 'LineWidth', 2);

xlabel('H2O2 (uM)')
ylabel('Flux (uM/S)')

pbaspect([1.05 1 1])

% yline(152.3, '--')
% yline(100, '--')
% yline(61.3, '--')

set(gca,'Fontsize', 22);

xlim([3e-3, 5])
ylim([2, 300])

%


%% --------------- Generate final_fig. 7A left panel - 2-parameter bifurcation k1 vs. k0 ------------------ %%
% Expect several hours to complete
param = default_param;
% H2O2, PRXSH, PRXSOH, PRXSS, PRXSO2H, TRXSS
init_off = [0, 100, 0, 0, 0, 0]; % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100, 0]; % initial concentrations for ON state;
num_steps = 300; 
k1_step_size = 300/num_steps;
k0_step_size = 300/num_steps;  

for i=1:1:num_steps + 1
    i
    param.k1 = 1e-6 + (i-1) * k1_step_size;
    k1(i) = param.k1;
    for j = 1:1:num_steps + 1
        j;
        param.k0 = 1e-6 + (j-1) * k0_step_size;
        k0(j) = param.k0;      

        [t,y] = ode15s('PTRS_ode', tspan, init_off, [], param);
        H2O2_switch_on(i,j) = y(end,1);

        [t,y] = ode15s('PTRS_ode', tspan, init_on, [], param);
        H2O2_switch_off(i,j) = y(end,1);
    end    
end

% save Two_par_k0k1_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0k1_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0k1_k1.mat k1
% save Two_par_k0k1_k0.mat k0


figure(5011)
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;
Bistability_mag = log10(H2O2_ratio);
surf(k0, k1, Bistability_mag);  

xlim([0, 300])
ylim([0, 300])

xlabel('k0')
ylabel('k1')
zlabel('Bistability magnitude')
shading interp 

set(gca,'Fontsize', 26);
pbaspect([1.26 1 1])
cb = colorbar; 
set(cb,'position',[.86 0.236 .04 .69])

hold on

% % Overlay XPP-AUT results to confirm
% [num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 1);    
% column_num = length(txt);
% row_num = length(num);
% figure(5011)
% z_offset = zeros(row_num,1) + max(max(H2O2_ratio));
% scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
% xlabel(txt(1))
% ylabel(txt(2))
% box on

%

%% --------------- Generate final_fig. 7A right panel - Flux analysis --------------------------------------------- %%
init = default_init;
param = default_param;
param.H2O2_switch = 0;
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

        init.PRXSH = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('PTRS_ode', tspan, y0, [], param);

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

    figure(5012)
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  
    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);
    

    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')
    
    set(gca,'Fontsize', 26);
    pbaspect([1.26 1 1])

    xlim([5e-3,2])
    xticks([0.01, 0.1, 1])

    ylim([1,1e3])
    yticks([1, 10, 100, 1000])

end

%

%% --------------- Generate final_fig. 7B left panel - 2-parameter bifurcation k2a vs. k0  ---------------- %%
%Expect several hours to complete
param = default_param;
% H2O2, PRXSH, PRXSOH, PRXSS, PRXSO2H, TRXSS
init_off = [0, 100, 0, 0, 0, 0]; % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100, 0]; % initial concentrations for ON state;

num_steps = 600; 
k2a_step_size = 800/num_steps;
k0_step_size = 600/num_steps;  

for i=1:1:num_steps + 1
    i
    param.k2a = 1e-6 + (i-1) * k2a_step_size;
    k2a(i) = param.k2a;
    
    for j = 1:1:num_steps + 1
        j;
        param.k0 = 1e-6 + (j-1) * k0_step_size;
        k0(j) = param.k0;      

        [t,y] = ode15s('PTRS_ode', tspan, init_off, [], param);
        H2O2_switch_on(i,j) = y(end,1);

        [t,y] = ode15s('PTRS_ode', tspan, init_on, [], param);
        H2O2_switch_off(i,j) = y(end,1);
    end    
end

% save Two_par_k0k2a_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0k2a_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0k2a_k2a.mat k2a
% save Two_par_k0k2a_k0.mat k0


figure(5021)
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;
Bistability_mag = log10(H2O2_ratio);
surf(k0, k2a, Bistability_mag);  
pbaspect([1.1 1 1])
xlim([1e-6, 600])
ylim([1e-6, 800])
xlabel('k0')
ylabel('k2a')
zlabel('Bistability magnitude')
shading interp 
colorbar
hold on

% % Overlay XPP-AUT results to confirm
% [num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 2); 
% column_num = length(txt);
% row_num = length(num);
% figure(5021)
% z_offset = zeros(row_num,1) + max(max(H2O2_ratio));
% scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
% xlabel(txt(1))
% ylabel(txt(2))
% box on

%

%% --------------- Generate final_fig. 7B right panel - Flux analysis --------------------------------------------- %%
init = default_init;
param = default_param;
param.H2O2_switch = 0;

k2a_vector = 2.^[-2:1:2] * default_param.k2a; 

for i = 1 : 1 : length(k2a_vector)   
    param.k2a = k2a_vector(i);       
    H2O2_vector = [];
    Flux_k1_vector = [];
    Flux_k3_vector = [];
    Flux_k5_vector = [];
    Flux_total_removal_vector = [];

    init.H2O2 = 0.001;
    while init.H2O2 <= 100
        
        init.PRXSH = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('PTRS_ode', tspan, y0, [], param);

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


    figure(5022) 
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  
    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);   

    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')

    set(gca,'Fontsize', 26);
    pbaspect([1.26 1 1])

    xlim([5e-3,2])
    xticks([0.01, 0.1, 1])

    ylim([1,1e3])
    yticks([1, 10, 100, 1000])
end

%

%% --------------- Generate final_fig. 7C left panel - 2-parameter bifurcation k2b vs. k0 ----------------- %%
%Expect several hours to complete
param = default_param;
% H2O2, PRXSH, PRXSOH, PRXSS, PRXSO2H, TRXSS
init_off = [0, 100, 0, 0, 0, 0]; % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100, 0]; % initial concentrations for ON state;

num_steps = 300; 
k2b_step_size = 2/num_steps; 
k0_step_size = 200/num_steps;  

for i=1:1:num_steps + 1;
    i
    param.k2b = 1e-6 + (i-1) * k2b_step_size;   
    k2b(i) = param.k2b;                         
    
    for j = 1:1:num_steps + 1;
        j;
        param.k0 = 1e-6 + (j-1) * k0_step_size;
        k0(j) = param.k0;      

        [t,y] = ode15s('PTRS_ode', tspan, init_off, [], param);
        H2O2_switch_on(i,j) = y(end,1);

        [t,y] = ode15s('PTRS_ode', tspan, init_on, [], param);
        H2O2_switch_off(i,j) = y(end,1);
    end    
end

% save Two_par_k0k2b_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0k2b_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0k2b_k2b.mat k2b
% save Two_par_k0k2b_k0.mat k0


figure(5031)                                        
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;
Bistability_mag = log10(H2O2_ratio);
surf(k0, k2b, Bistability_mag);                         
pbaspect([1.1 1 1])

xlim([1e-6,200])
ylim([1e-6,2])

xlabel('k0')
ylabel('k2b')                                       
zlabel('Bistability magnitude')
shading interp 
colorbar
hold on
yline(default_param.k2b,'--r')
% 
% % Overlay XPP-AUT results to confirm
% [num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 3);
% column_num = length(txt);
% row_num = length(num);
% figure(5031)                                         
% z_offset = zeros(row_num,1) + max(max(H2O2_ratio));
% scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
% xlabel(txt(1))
% ylabel(txt(2))
% box on


%
%% --------------- Generate final_fig. 7C right panel - Flux analysis --------------------------------------------- %%
init = default_init;
param = default_param;
param.H2O2_switch = 0;

k2b_vector = 2.^[-2:1:2] * default_param.k2b; 

for i = 1 : 1 : length(k2b_vector)   
    param.k2b = k2b_vector(i);      
    H2O2_vector = [];
    Flux_k1_vector = [];
    Flux_k3_vector = [];
    Flux_k5_vector = [];
    Flux_total_removal_vector = [];

    init.H2O2 = 0.001;
    while init.H2O2 <= 100
        
        init.PRXSH = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('PTRS_ode', tspan, y0, [], param);

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
    

    figure(5032) 
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  
    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);    

    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')

    set(gca,'Fontsize', 26);
    pbaspect([1.26 1 1])

    xlim([5e-3,2])
    xticks([0.01, 0.1, 1])

    ylim([1,1e3])
    yticks([1, 10, 100, 1000])
end
%

%% --------------- Generate final_fig. 7D left panel - 2-parameter bifurcation k3 vs. k0 ------------------ %%
%Expect several hours to complete
param = default_param;
% H2O2, PRXSH, PRXSOH, PRXSS, PRXSO2H, TRXSS
init_off = [0, 100, 0, 0, 0, 0]; % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100, 0]; % initial concentrations for ON state;

num_steps = 300; 
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

        [t,y] = ode15s('PTRS_ode', tspan, init_off, [], param);
        H2O2_switch_on(i,j) = y(end,1);

        [t,y] = ode15s('PTRS_ode', tspan, init_on, [], param);
        H2O2_switch_off(i,j) = y(end,1);
    end    
end

% save Two_par_k0k3_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0k3_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0k3_k3.mat k3
% save Two_par_k0k3_k0.mat k0


figure(5041)                                        
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;
Bistability_mag = log10(H2O2_ratio);
surf(k0, k3, Bistability_mag);                          
pbaspect([1.1 1 1])
xlim([1e-6,220])
ylim([1e-6,6])
xlabel('k0')
ylabel('k3')                                       
zlabel('Bistability magnitude')
shading interp 
colorbar
hold on
yline(default_param.k3,'--r')
% % Overlay XPP-AUT results to confirm
% [num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 4);    
% column_num = length(txt);
% row_num = length(num);
% figure(5041)                                                   
% z_offset = zeros(row_num,1) + max(max(Bistability_mag));
% scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
% xlabel(txt(1))
% ylabel(txt(2))
% box on
%
%% --------------- Generate final_fig. 7D right panel - Flux analysis --------------------------------------------- %%
init = default_init;
param = default_param;
param.H2O2_switch = 0;

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
        
        init.PRXSH = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('PTRS_ode', tspan, y0, [], param);

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
    
    figure(5042)                                                
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  

    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);
       
    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')

    set(gca,'Fontsize', 26);
    pbaspect([1.26 1 1])

    xlim([5e-3,2])
    xticks([0.01, 0.1, 1])

    ylim([1,1e3])
    yticks([1, 10, 100, 1000])
end
%

%% --------------- Generate final_fig. 7E left panel - 2-parameter bifurcation k4f vs. k0 ------------- %%
%Expect several hours to complete
param = default_param;
% H2O2, PRXSH, PRXSOH, PRXSS, PRXSO2H, TRXSS
init_off = [0, 100, 0, 0, 0, 0];               % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100, 0];               % initial concentrations for ON state;

num_steps = 300; 
k4f_step_size = 0.005/num_steps;             
k0_step_size = 200/num_steps;  

for i=1:1:num_steps + 1;
    i
    param.k4f = 1e-6 + (i-1) * k4f_step_size;  
    k4f(i) = param.k4f;                         
    
    for j = 1:1:num_steps + 1;
        j;
        param.k0 = 1e-6 + (j-1) * k0_step_size;
        k0(j) = param.k0;      
        
        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t,y] = ode15s('PTRS_ode', tspan, init_off, options, param);
        H2O2_switch_on(i,j) = y(end,1);

%         figure(100)
%         plot(t, y(:,1))
%         xlabel ('Time')
%         ylabel ('H2O2')
%         hold on
        
        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t,y] = ode15s('PTRS_ode', tspan, init_on, options, param);
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


figure(5051)                                         
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;
Bistability_mag = log10(H2O2_ratio);
surf(k0, k4f, Bistability_mag);   
pbaspect([1.1 1 1])
xlim([1e-6, 200])
ylim([1e-6, 0.005])
xlabel('k0')
ylabel('k4f')                                       
zlabel('Bistability magnitude')
shading interp 
colorbar
hold on

yline(default_param.k4f,'--r')

% % Overlay XPP-AUT results to confirm
% [num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 5);   
% column_num = length(txt);
% row_num = length(num);
% figure(5051)                                                                                      
% z_offset = zeros(row_num,1) + max(max(H2O2_ratio));
% scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
% xlabel(txt(1))
% ylabel(txt(2))
% box on
%
%% --------------- Generate final_fig. 7E right panel - Flux analysis --------------------------------------------- %%
init = default_init;
param = default_param;
param.H2O2_switch = 0;

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
        
        init.PRXSH = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('PTRS_ode', tspan, y0, [], param);

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
    

    figure(5052)                                               
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  
    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);
      
    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')

    set(gca,'Fontsize', 26);
    pbaspect([1.26 1 1])

    xlim([5e-3,2])
    xticks([0.01, 0.1, 1])

    ylim([1,1e3])
    yticks([1, 10, 100, 1000])    
end
%

%% --------------- Generate final_fig. 7F left panel - 2-parameter bifurcation k4b vs. k0 ------------- %%
%Expect several hours to complete
param = default_param;
% H2O2, PRXSH, PRXSOH, PRXSS, PRXSO2H, TRXSS
init_off = [0, 100, 0, 0, 0, 0]; % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100, 0]; % initial concentrations for ON state;

num_steps = 300; 
k4b_step_size = 0.1/num_steps;
k0_step_size = 200/num_steps;  

for i=1:1:num_steps + 1;
    i
    param.k4b = 1e-6 + (i-1) * k4b_step_size;     
    k4b(i) = param.k4b;                           
    
    for j = 1:1:num_steps + 1;
        j;
        param.k0 = 1e-6 + (j-1) * k0_step_size;
        k0(j) = param.k0;

        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t,y] = ode15s('PTRS_ode', tspan, init_off, options, param);
        H2O2_switch_on(i,j) = y(end,1);

        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t,y] = ode15s('PTRS_ode', tspan, init_on, options, param);
        H2O2_switch_off(i,j) = y(end,1);
    end    
end

% save Two_par_k0k4b_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0k4b_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0k4b_k4b.mat k4b
% save Two_par_k0k4b_k0.mat k0


figure(5061)                                        
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;
Bistability_mag = log10(H2O2_ratio);     
surf(k0, k4b, Bistability_mag);       
pbaspect([1.1 1 1])
xlim([1e-6, 200])
ylim([1e-6, 0.1])
xlabel('k0')
ylabel('k4b')                                        
zlabel('Bistability magnitude')
shading interp 
colorbar
hold on

yline(default_param.k4b,'--r')


% % Overlay XPP-AUT results to confirm
% [num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 6);   
% column_num = length(txt);
% row_num = length(num);
% figure(5061)                                                   
% z_offset = zeros(row_num,1) + max(max(H2O2_ratio));
% scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
% xlabel(txt(1))
% ylabel(txt(2))
% box on
% 

%% --------------- Generate final_fig. 7F right panel - Flux analysis --------------------------------------------- %%
init = default_init;
param = default_param;
param.H2O2_switch = 0;

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
        
        init.PRXSH = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('PTRS_ode', tspan, y0, [], param);

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
    
    figure(5063)                                               
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  
    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);

    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')

    set(gca,'Fontsize', 26);
    pbaspect([1.26 1 1])

    xlim([5e-3,2])
    xticks([0.01, 0.1, 1])

    ylim([1,1e3])
    yticks([1, 10, 100, 1000])  
end
%

%% --------------- Generate final_fig. 7G left panel - 2-parameter bifurcation k4c vs. k0 ------------- %%
%Expect several hours to complete
param = default_param;
tspan1 = [0:100:1e7]; % need 1e7 to reach steady-state

% H2O2, PRXSH, PRXSOH, PRXSS, PRXSO2H, TRXSS
init_off = [0, 100, 0, 0, 0, 0]; % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100, 0]; % initial concentrations for ON state;

num_steps = 300; 
k4c_step_size = log10(10/1e-6)/num_steps;       
k0_step_size = 350/num_steps;  

for i=1:1:num_steps + 1;
    i
    param.k4c = 10^(-6 + (i-1) * k4c_step_size);   
    k4c(i) = param.k4c;                         
    
    for j = 1:1:num_steps + 1;
        j;
        param.k0 = 1e-5 + (j-1) * k0_step_size;
        k0(j) = param.k0;      

        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t,y] = ode15s('PTRS_ode', tspan1, init_off, options, param);
        H2O2_switch_on(i,j) = y(end,1);

        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t,y] = ode15s('PTRS_ode', tspan1, init_on, options, param);
        H2O2_switch_off(i,j) = y(end,1);
    end    
end

% save Two_par_k0k4c_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0k4c_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0k4c_k4c.mat k4c
% save Two_par_k0k4c_k0.mat k0


figure(5071)                                        
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;
Bistability_mag = log10(H2O2_ratio);
surf(k0, k4c, Bistability_mag);    
pbaspect([1.1 1 1])
xlabel('k0')
ylabel('k4c')                                        
zlabel('Bistability magnitude')
shading interp 
colorbar
hold on

yline(default_param.k4c,'--r')

% % Overlay XPP-AUT results to confirm
% [num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 7); 
% row_num = length(num);
% figure(5071)                                         
% z_offset = zeros(row_num,1) + max(max(H2O2_ratio));
% scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
% xlabel(txt(1))
% ylabel(txt(2))
% box on
% xlim([1e-6,150])
% ylim([1e-6,10])
% 

%% --------------- Generate final_fig. 7G right panel - Flux analysis --------------------------------------------- %%
init = default_init;
param = default_param;
param.H2O2_switch = 0;

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
        
        init.PRXSH = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('PTRS_ode', tspan, y0, [], param);

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
    

    figure(5072) 
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  
    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);

    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')

    set(gca,'Fontsize', 26);
    pbaspect([1.26 1 1])

    xlim([5e-3,2])
    xticks([0.01, 0.1, 1])

    ylim([1,1e3])
    yticks([1, 10, 100, 1000])  
end
%

%% --------------- Generate final_fig. 7H left panel - 2-parameter bifurcation SRXtot vs. k0 ---------- %%
%Expect several hours to complete
param = default_param;
tspan1 = [0:100:1e7]; % need 1e7 to reach steady-state

% H2O2, PRXSH, PRXSOH, PRXSS, PRXSO2H, TRXSS
init_off = [0, 100, 0, 0, 0, 0]; % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100, 0]; % initial concentrations for ON state;

num_steps = 300; 
SRXtot_step_size = log10(100/1e-3)/num_steps; 
k0_step_size = 500/num_steps;  

for i=1:1:num_steps + 1;
    i
    param.SRXtot = 10^(-3 + (i-1) * SRXtot_step_size);        
    SRXtot(i) = param.SRXtot;                                 
    
    for j = 1:1:num_steps + 1;
        j;
        param.k0 = 1e-6 + (j-1) * k0_step_size; 
        k0(j) = param.k0;      

        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t,y] = ode15s('PTRS_ode', tspan1, init_off, options, param);
        H2O2_switch_on(i,j) = y(end,1);

        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t,y] = ode15s('PTRS_ode', tspan1, init_on, options, param);
        H2O2_switch_off(i,j) = y(end,1);
    end    
end

% save Two_par_k0SRXtot_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0SRXtot_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0SRXtot_SRXtot.mat SRXtot
% save Two_par_k0SRXtot_k0.mat k0


figure(5081)                                                
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;   
Bistability_mag = log10(H2O2_ratio);
surf(k0, SRXtot, Bistability_mag);                         
pbaspect([1.1 1 1])
xlim([1e-6, 500])
ylim([1e-6, 20])
xlabel('k0')
ylabel('SRXtot')                                         
zlabel('Bistability magnitude')
shading interp 
colorbar
hold on

yline(default_param.SRXtot,'--r')

% % Overlay XPP-AUT results to confirm
% [num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 8); 
% column_num = length(txt);
% row_num = length(num);
% figure(5081)                                                 
% z_offset = zeros(row_num,1) + max(max(H2O2_ratio));
% scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
% xlabel(txt(1))
% ylabel(txt(2))
% box on
% 

%% --------------- Generate final_fig. 7H right panel - Flux analysis --------------------------------------------- %%
init = default_init;
param = default_param;
param.H2O2_switch = 0;

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
        
        init.PRXSH = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('PTRS_ode', tspan, y0, [], param);

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
    
    figure(5082) 
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  
    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);
       
    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')

    set(gca,'Fontsize', 26);
    pbaspect([1.26 1 1])

    xlim([5e-3,2])
    xticks([0.01, 0.1, 1])

    ylim([1,1e3])
    yticks([1, 10, 100, 1000]) 
end
%

%% --------------- Generate final_fig. 7I left panel - 2-parameter bifurcation k5 vs. k0 -------------- %%
%Expect several hours to complete
param = default_param;
% H2O2, PRXSH, PRXSOH, PRXSS, PRXSO2H, TRXSS
init_off = [0, 100, 0, 0, 0, 0]; % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100, 0]; % initial concentrations for ON state;

num_steps = 600; 
k5_step_size = 1600/num_steps;
k0_step_size = 250/num_steps;  

for i=1:1:num_steps + 1
    i
    param.k5 = 1e-6 + (i-1) * k5_step_size;    
    k5(i) = param.k5;                          
    
    for j = 1:1:num_steps + 1
        j;
        param.k0 = 1e-6 + (j-1) * k0_step_size;
        k0(j) = param.k0;      

        [t,y] = ode15s('PTRS_ode', tspan, init_off, [], param);
        H2O2_switch_on(i,j) = y(end,1);

        [t,y] = ode15s('PTRS_ode', tspan, init_on, [], param);
        H2O2_switch_off(i,j) = y(end,1);
    end    
end

% save Two_par_k0k5_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0k5_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0k5_k5.mat k5
% save Two_par_k0k5_k0.mat k0


figure(5091)                                       
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;
Bistability_mag = log10(H2O2_ratio);
surf(k0, k5, Bistability_mag);  
pbaspect([1.1 1 1])
xlim([1e-6, 250])
ylim([1e-6, 2000])
xlabel('k0')
ylabel('k5')                                      
zlabel('Bistability magnitude')
shading interp 
colorbar
hold on

yline(default_param.k5,'--r')

% % Overlay XPP-AUT results to confirm
% [num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 9); 
% column_num = length(txt);
% row_num = length(num);
% figure(5091)                                        
% z_offset = zeros(row_num,1) + max(max(H2O2_ratio));
% scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
% xlabel(txt(1))
% ylabel(txt(2))
% box on
% 

%% --------------- Generate final_fig. 7I right panel - Flux analysis --------------------------------------------- %%
init = default_init;
param = default_param;
param.H2O2_switch = 0;

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
        
        init.PRXSH = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('PTRS_ode', tspan, y0, [], param);

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
    

    figure(5092) 
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  
    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);
 
    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')

    set(gca,'Fontsize', 26);
    pbaspect([1.26 1 1])

    xlim([5e-3,2])
    xticks([0.01, 0.1, 1])

    ylim([1,1e3])
    yticks([1, 10, 100, 1000]) 
end
%

%% --------------- Generate final_fig. 7J left panel - 2-parameter bifurcation PRXtot vs. k0 ---------- %%
%Expect several hours to complete
param = default_param;

num_steps = 300; 
PRXtot_step_size = 1000/num_steps;
k0_step_size = 350/num_steps;  

for i=1:1:num_steps + 1;
    i
    param.PRXtot = 1e-6 + (i-1) * PRXtot_step_size; 
    PRXtot(i) = param.PRXtot;                             
    
    % H2O2, PRXSH, PRXSOH, PRXSS, PRXSO2H, TRXSS
    init_off = [0, PRXtot(i), 0, 0, 0, 0]; % initial concentrations for OFF state;
    init_on = [10, 0, 0, 0, PRXtot(i), 0]; % initial concentrations for ON state;

    for j = 1:1:num_steps + 1;
        j;
        param.k0 = 1e-6 + (j-1) * k0_step_size;   
        k0(j) = param.k0;      

        [t,y] = ode15s('PTRS_ode', tspan, init_off, [], param);
        H2O2_switch_on(i,j) = y(end,1);

        [t,y] = ode15s('PTRS_ode', tspan, init_on, [], param);
        H2O2_switch_off(i,j) = y(end,1);
    end    
end

% save Two_par_k0PRXtot_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0PRXtot_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0PRXtot_PRXtot.mat PRXtot
% save Two_par_k0PRXtot_k0.mat k0

figure(5101)                                              
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;
Bistability_mag = log10(H2O2_ratio);
surf(k0, PRXtot, Bistability_mag);                        
pbaspect([1.1 1 1])
xlim([1e-6, 350])
ylim([1e-6, 1000])
xlabel('k0')
ylabel('PRXtot')                                      
zlabel('Bistability magnitude')
shading interp 
colorbar
hold on

yline(default_param.PRXtot,'--r')

% % Overlay XPP-AUT results to confirm
% [num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 10); 
% column_num = length(txt);
% row_num = length(num);
% figure(5101)                                               
% z_offset = zeros(row_num,1) + max(max(H2O2_ratio));
% scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
% xlabel(txt(1))
% ylabel(txt(2))
% box on
%  

%% --------------- Generate final_fig. 7J right panel - Flux analysis --------------------------------------------- %%
init = default_init;
param = default_param;
param.H2O2_switch = 0;

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
        
        init.PRXSH = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('PTRS_ode', tspan, y0, [], param);

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
    

    figure(5102) 
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  
    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);
       
    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')

    set(gca,'Fontsize', 26);
    pbaspect([1.26 1 1])

    xlim([5e-3,2])
    xticks([0.01, 0.1, 1])

    ylim([1,1e3])
    yticks([1, 10, 100, 1000]) 
end
%
%% --------------- Generate final_fig. S4A left panel - 2-parameter bifurcation k2c vs. k0 ------------------------ %%
%Expect several hours to complete
param = default_param;
% H2O2, PRXSH, PRXSOH, PRXSS, PRXSO2H, TRXSS
init_off = [0, 100, 0, 0, 0, 0]; % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100, 0]; % initial concentrations for ON state;

num_steps = 300; 
k2c_step_size = 500/num_steps; 
k0_step_size = 200/num_steps;  

for i=1:1:num_steps + 1;
    i
    param.k2c = 1e-6 + (i-1) * k2c_step_size;   
    k2c(i) = param.k2c;                         
    
    for j = 1:1:num_steps + 1;
        j;
        param.k0 = 1e-6 + (j-1) * k0_step_size;
        k0(j) = param.k0;      

        [t,y] = ode15s('PTRS_ode', tspan, init_off, [], param);
        H2O2_switch_on(i,j) = y(end,1);

        [t,y] = ode15s('PTRS_ode', tspan, init_on, [], param);
        H2O2_switch_off(i,j) = y(end,1);
    end    
end

% save Two_par_k0k2c_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0k2c_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0k2c_k2c.mat k2c
% save Two_par_k0k2c_k0.mat k0


figure(5111)                                        
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;
Bistability_mag = log10(H2O2_ratio);
surf(k0, k2c, Bistability_mag);                         
pbaspect([1.1 1 1])

xlim([1e-6,200])
ylim([1e-6,2])

xlabel('k0')
ylabel('k2c')                                       
zlabel('Bistability magnitude')
shading interp 
colorbar
hold on

yline(default_param.k2c,'--r')


% 
% % Overlay XPP-AUT results to confirm
% [num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 3);
% column_num = length(txt);
% row_num = length(num);
% figure(5031)                                         
% z_offset = zeros(row_num,1) + max(max(H2O2_ratio));
% scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
% xlabel(txt(1))
% ylabel(txt(2))
% box on
%
%% --------------- Generate final_fig. S4A right panel - Flux analysis -------------------------------------------- %%
init = default_init;
param = default_param;
param.H2O2_switch = 0;

k2c_vector = 2.^[-2:1:2] * default_param.k2c;   

for i = 1 : 1 : length(k2c_vector)                   
    param.k2c = k2c_vector(i);                 
    H2O2_vector = [];
    Flux_k1_vector = [];
    Flux_k3_vector = [];
    Flux_k5_vector = [];
    Flux_total_removal_vector = [];

    init.H2O2 = 0.001;
    while init.H2O2 <= 100
        
        init.PRXSH = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('PTRS_ode', tspan, y0, [], param);

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
    

    figure(5112) 
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  
    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);
       
    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')

    set(gca,'Fontsize', 26);
    pbaspect([1.26 1 1])

    xlim([5e-3,2])
    xticks([0.01, 0.1, 1])

    ylim([1,1e3])
    yticks([1, 10, 100, 1000]) 
end
%
%% --------------- Generate final_fig. S4B left panel - 2-parameter bifurcation Km2c vs. k0 ----------------------- %%
%Expect several hours to complete
param = default_param;
% H2O2, PRXSH, PRXSOH, PRXSS, PRXSO2H, TRXSS
init_off = [0, 100, 0, 0, 0, 0]; % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100, 0]; % initial concentrations for ON state;

num_steps = 300; 
Km2c_step_size = 300/num_steps; 
k0_step_size = 200/num_steps;  

for i=1:1:num_steps + 1;
    i
    param.Km2c = 1e-6 + (i-1) * Km2c_step_size;   
    Km2c(i) = param.Km2c;                         
    
    for j = 1:1:num_steps + 1;
        j;
        param.k0 = 1e-6 + (j-1) * k0_step_size;
        k0(j) = param.k0;      

        [t,y] = ode15s('PTRS_ode', tspan, init_off, [], param);
        H2O2_switch_on(i,j) = y(end,1);

        [t,y] = ode15s('PTRS_ode', tspan, init_on, [], param);
        H2O2_switch_off(i,j) = y(end,1);
    end    
end

% save Two_par_k0k2c_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0k2c_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0k2c_k2c.mat k2c
% save Two_par_k0k2c_k0.mat k0


figure(5121)                                        
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;
Bistability_mag = log10(H2O2_ratio);
surf(k0, Km2c, Bistability_mag);                         
pbaspect([1.1 1 1])

xlim([1e-6,200])
ylim([1e-6,2])

xlabel('k0')
ylabel('Km2c')                                       
zlabel('Bistability magnitude')
shading interp 
colorbar
hold on

yline(default_param.Km2c,'--r')

% 
% % Overlay XPP-AUT results to confirm
% [num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 3);
% column_num = length(txt);
% row_num = length(num);
% figure(5031)                                         
% z_offset = zeros(row_num,1) + max(max(H2O2_ratio));
% scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
% xlabel(txt(1))
% ylabel(txt(2))
% box on
%
%% --------------- Generate final_fig. S4B right panel - Flux analysis -------------------------------------------- %%
init = default_init;
param = default_param;
param.H2O2_switch = 0;
 
Km2c_vector = [0.1 0.2 1 5 10] * default_param.Km2c;   

for i = 1 : 1 : length(Km2c_vector)                   
    param.Km2c = Km2c_vector(i);                 
    H2O2_vector = [];
    Flux_k1_vector = [];
    Flux_k3_vector = [];
    Flux_k5_vector = [];
    Flux_total_removal_vector = [];

    init.H2O2 = 0.001;
    while init.H2O2 <= 100
        
        init.PRXSH = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('PTRS_ode', tspan, y0, [], param);

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
    
    figure(5122) 
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  
    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);
       
    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')

    set(gca,'Fontsize', 26);
    pbaspect([1.26 1 1])

    xlim([5e-3,2])
    xticks([0.01, 0.1, 1])

    ylim([1,1e3])
    yticks([1, 10, 100, 1000]) 
end
%


%% --------------- Generate final_fig. S4C left panel - 2-parameter bifurcation TRXtot vs. k0 ------------- %%
% Expect several hours to complete
param = default_param;
% H2O2, PRXSH, PRXSOH, PRXSS, PRXSO2H, TRXSS
init_off = [0, 100, 0, 0, 0, 0]; % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100, 0]; % initial concentrations for ON state;

num_steps = 300; 
TRXtot_step_size = 200/num_steps;
k0_step_size = 200/num_steps;  

for i=1:1:num_steps + 1
    i
    param.TRXtot = 1e-6 + (i-1) * TRXtot_step_size;
    TRXtot(i) = param.TRXtot;
    for j = 1:1:num_steps + 1
        j;
        param.k0 = 1e-6 + (j-1) * k0_step_size;
        k0(j) = param.k0;      

        [t,y] = ode15s('PTRS_ode', tspan, init_off, [], param);
        H2O2_switch_on(i,j) = y(end,1);

        [t,y] = ode15s('PTRS_ode', tspan, init_on, [], param);
        H2O2_switch_off(i,j) = y(end,1);
    end    
end

% save Two_par_k0TRXtot_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0TRXtot_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0TRXtot_TRXtot.mat TRXtot
% save Two_par_k0TRXtot_k0.mat k0


figure(5131)
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;
Bistability_mag = log10(H2O2_ratio);
surf(k0, TRXtot, Bistability_mag);  

xlim([0, 300])
ylim([0, 100])

xlabel('k0')
ylabel('TRXtot')
zlabel('Bistability magnitude')
shading interp 

set(gca,'Fontsize', 26);
pbaspect([1.26 1 1])
cb = colorbar; 
set(cb,'position',[.86 0.236 .04 .69])

hold on

% % Overlay XPP-AUT results to confirm
% [num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 1);    
% column_num = length(txt);
% row_num = length(num);
% figure(5011)
% z_offset = zeros(row_num,1) + max(max(H2O2_ratio));
% scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
% xlabel(txt(1))
% ylabel(txt(2))
% box on
%
%% --------------- Generate final_fig. S4C right panel - Flux analysis -------------------------------------------- %%
init = default_init;
param = default_param;
param.H2O2_switch = 0;

TRXtot_vector = [0.1 0.2 1 5 10] * default_param.TRXtot;   

for i = 1 : 1 : length(TRXtot_vector)                   
    param.TRXtot = TRXtot_vector(i);                 
    H2O2_vector = [];
    Flux_k1_vector = [];
    Flux_k3_vector = [];
    Flux_k5_vector = [];
    Flux_total_removal_vector = [];

    init.H2O2 = 0.001;
    while init.H2O2 <= 100
        
        init.PRXSH = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('PTRS_ode', tspan, y0, [], param);

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
    

    figure(5132) 
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  
    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);
       
    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')

    set(gca,'Fontsize', 26);
    pbaspect([1.26 1 1])

    xlim([5e-3,2])
    xticks([0.01, 0.1, 1])

    ylim([1,1e3])
    yticks([1, 10, 100, 1000]) 
end
%
%% --------------- Generate final_fig. S4D left panel - 2-parameter bifurcation TRtot vs. k0 -------------- %%
% Expect several hours to complete
param = default_param;
% H2O2, PRXSH, PRXSOH, PRXSS, PRXSO2H, TRXSS
init_off = [0, 100, 0, 0, 0, 0]; % initial concentrations for OFF state;
init_on = [10, 0, 0, 0, 100, 0]; % initial concentrations for ON state;

num_steps = 300; 
TRtot_step_size = 10/num_steps;
k0_step_size = 200/num_steps;  

for i=1:1:num_steps + 1
    i
    param.TRtot = 1e-6 + (i-1) * TRtot_step_size;
    TRtot(i) = param.TRtot;
    for j = 1:1:num_steps + 1
        j;
        param.k0 = 1e-6 + (j-1) * k0_step_size;
        k0(j) = param.k0;      

        [t,y] = ode15s('PTRS_ode', tspan, init_off, [], param);
        H2O2_switch_on(i,j) = y(end,1);

        [t,y] = ode15s('PTRS_ode', tspan, init_on, [], param);
        H2O2_switch_off(i,j) = y(end,1);
    end    
end

% save Two_par_k0TRtot_H2O2_switch_on.mat H2O2_switch_on
% save Two_par_k0TRtot_H2O2_switch_off.mat H2O2_switch_off
% save Two_par_k0TRtot_TRtot.mat TRtot
% save Two_par_k0TRtot_k0.mat k0

figure(5141)
H2O2_ratio = H2O2_switch_off./H2O2_switch_on;
Bistability_mag = log10(H2O2_ratio);
surf(k0, TRtot, Bistability_mag);  

xlim([0, 300])
ylim([0, 300])

xlabel('k0')
ylabel('TRtot')
zlabel('Bistability magnitude')
shading interp 

set(gca,'Fontsize', 26);
pbaspect([1.26 1 1])
cb = colorbar; 
set(cb,'position',[.86 0.236 .04 .69])

hold on

% % Overlay XPP-AUT results to confirm
% [num, txt] = xlsread('XPP-AUT\Bistability_bifurcation_two_parameters.xlsx', 1);    
% column_num = length(txt);
% row_num = length(num);
% figure(5011)
% z_offset = zeros(row_num,1) + max(max(H2O2_ratio));
% scatter3(num(:, 1), num(:, 2), z_offset,'red','SizeData', 2);
% xlabel(txt(1))
% ylabel(txt(2))
% box on
%
%% --------------- Generate final_fig. S4D right panel - Flux analysis -------------------------------------------- %%
init = default_init;
param = default_param;
param.H2O2_switch = 0;

TRtot_vector = 2.^[-2:1:2] * default_param.TRtot;   

for i = 1 : 1 : length(TRtot_vector)                   
    param.TRtot = TRtot_vector(i);                 
    H2O2_vector = [];
    Flux_k1_vector = [];
    Flux_k3_vector = [];
    Flux_k5_vector = [];
    Flux_total_removal_vector = [];

    init.H2O2 = 0.001;
    while init.H2O2 <= 100
        
        init.PRXSH = param.PRXtot - init.PRXSO2H;

        H2O2_vector = [H2O2_vector, init.H2O2];

        y0 = cell2mat(struct2cell(init));

        [t,y]=ode15s('PTRS_ode', tspan, y0, [], param);

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
    
    figure(5142) 
    lh = loglog(H2O2_vector, Flux_total_removal_vector, 'Color', [0.00,0.45,0.74], 'LineWidth', 5); 
    lh.Color = [lh.Color i*0.2];
    hold on  

    lh = loglog(H2O2_vector, Flux_k1_vector, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];
  
    lh = loglog(H2O2_vector, Flux_k3_vector, 'Color', [0.93,0.69,0.13], 'LineWidth', 3);
    lh.Color = [lh.Color i*0.2];

    loglog(H2O2_vector, Flux_k5_vector, 'Color', [0.47,0.67,0.19], 'LineWidth', 3);
       
    xlabel('H2O2 (uM)')
    ylabel('Flux (uM/S)')

    set(gca,'Fontsize', 26);
    pbaspect([1.26 1 1])

    xlim([5e-3,2])
    xticks([0.01, 0.1, 1])

    ylim([1,1e3])
    yticks([1, 10, 100, 1000]) 

end

%



%% --------------- Generate final_fig. S5A-S5F - Time course for varying k0 --------------------------------------- %%
y0 = cell2mat(struct2cell(default_init)); %
param = default_param;

tspan0 = [0 : 100 : 3600*10]; % (0) run to steady state; 3600*10 seconds;

tspan1 = [0 : 100 : 3600*100];     % (1) run for 2e5 seconds with k0 = 10;
tspan2 = [3600*100 : 100 : 3600*200];   % (2) run with different k0 values(20:20:200) between 2e5-7e5 seconds;
tspan3 = [3600*200 : 100 : 3600*300];   % (3) run for 7e5-1e6 seconds with k0 = 10;

% (0) run to steady state; 3600*10 seconds;
param.k0 = 10;
[t,y] = ode15s('PTRS_ode', tspan0, y0, [], param); 
y1_0 = y(end,:);  


% run with different k0 values(80:20:200);
num_steps = 7; 
k0_vector = [];
for j = 1:1:num_steps 
    j
    % (1) run for 2e5 seconds with k0 = 10;
    param.k0 = 10;
    [t,y1] = ode15s('PTRS_ode', tspan1, y1_0, [], param); 
    y2_0 = y1(end,:);  

    % (2) run with different k0 values(80:20:200) between 2e5-7e5 seconds;
    param.k0 = 60 + 20*j;
    k0_vector = [k0_vector, param.k0];
    [t,y2] = ode15s('PTRS_ode', tspan2, y2_0, [], param); 
    y3_0 = y2(end,:);  

    % (3) run for 7e5-1e6 seconds with k0 = 10;
    param.k0 = 10;
    [t,y3] = ode15s('PTRS_ode', tspan3, y3_0, [], param); 
    
    % Plotting time course;
    t123 = [tspan1, tspan2, tspan3];
    y123 = [y1; y2; y3];

    pb_ratio = 1.8;
    font_size = 22;
    line_width = 3;

    % H2O2, PRXSH, PRXSOH, PRXSS, PRXSO2H, TRXSS
    figure(601)
    semilogy(t123/3600, y123(:,1), 'LineWidth', line_width)
    xlabel ('Time (h)')
    ylabel ('H2O2')
    pbaspect([pb_ratio 1 1])
    set(gca,'fontsize',font_size);
    hold on  

    figure(602)
    semilogy(t123/3600, y123(:,2), 'LineWidth', line_width)
    xlabel ('Time (h)')
    ylabel ('PRXSH')
    pbaspect([pb_ratio 1 1])
    set(gca,'fontsize',font_size);
    hold on

    figure(603)
    semilogy(t123/3600, y123(:,3), 'LineWidth', line_width)
    xlabel ('Time (h)')
    ylabel ('PRXSOH')
    pbaspect([pb_ratio 1 1])
    set(gca,'fontsize',font_size);
    hold on

    figure(604)
    semilogy(t123/3600, y123(:,4), 'LineWidth', line_width)
    xlabel ('Time (h)')
    ylabel ('PRXSS')
    pbaspect([pb_ratio 1 1])
    set(gca,'fontsize',font_size);
    hold on

    figure(605)
    semilogy(t123/3600, y123(:,5), 'LineWidth', line_width)
    xlabel ('Time (h)')
    ylabel ('PRXSO2H')
    pbaspect([pb_ratio 1 1])
    set(gca,'fontsize',font_size);
    hold on

    figure(606)
    semilogy(t123/3600, y123(:,6), 'LineWidth', line_width)
    xlabel ('Time (h)')
    ylabel ('TRXSS')
    pbaspect([pb_ratio 1 1])
    set(gca,'fontsize',font_size);
    hold on

end

%

%% --------------- Generate figure 9C and 9D -------------------------------------------------------------------- %%
y0 = cell2mat(struct2cell(default_init));
param = default_param;
param.k0 = 40; % run with different param.k0 values (40, 80, 160);
tspan0 = [0 : 1 : 3600*100]; 

tspan1 = [0 : 1 : 3600*100]; 
tspan2 = [3600*100 : 1 : 3600*600]; 
tspan3 = [3600*600 : 1 : 3600*700]; 

% (0) run to steady state;
param.SRXtot = 5;
[t,y] = ode15s('PTRS_ode', tspan0, y0, [], param); 
y1_0 = y(end,:);  

% run with different param.SRXtot values;
SRXtot_vector = [5, 0.5, 0.25, 0.1, 0.05, 0.01, 0.005, 0];
num_steps = length(SRXtot_vector); 

for j = 1:1:num_steps 
    j
    % (1) run for 100 seconds with param.SRXtot = 0.6;
    param.SRXtot = 5;
    [t,y1] = ode15s('PTRS_ode', tspan1, y1_0, [], param); 
    y2_0 = y1(end,:);  

    % (2) run with different param.SRXtot values between 100-300 seconds;
    param.SRXtot = SRXtot_vector(j);
    [t,y2] = ode15s('PTRS_ode', tspan2, y2_0, [], param); 
    y3_0 = y2(end,:);  

    % (3) run for 300-400 seconds with k0 = 10;
    param.SRXtot = 5;
    [t,y3] = ode15s('PTRS_ode', tspan3, y3_0, [], param); 
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
    set(gca,'fontsize',font_size);
    hold on  

    % PRXSO2H
    figure(913)
    plot(t123/3600, y123(:,5), 'LineWidth', line_width)
    xlabel ('Time (h)')
    ylabel ('PRXSO2H')
    pbaspect([pb_ratio 1 1])
    set(gca,'fontsize',font_size);
    hold on  
    
end

