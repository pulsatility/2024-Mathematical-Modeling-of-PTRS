%% --------------- Generate Fig. 8A and 8B ------------------------- %%

for i = 1 : 1 : 4
    [num, txt] = xlsread('Ultrasensitivity_k0_H2O2_PRX1 & PRX2 & PRX3', i);
    column_num = length(txt);
    
    % Obtain branch index: 
    branch_1_index = find(num(:, column_num) == 1); 
    branch_2_index = find(num(:, column_num) == 2);
    branch_3_index = find(num(:, column_num) == 3);
    branch_4_index = find(num(:, column_num) == 4);
    branch_5_index = find(num(:, column_num) == 5);

    % Obtain branch data: 
    branch_1 = num(branch_1_index, [1, 2]);
    branch_2 = num(branch_2_index, [1, 2]);
    branch_3 = num(branch_3_index, [1, 2]);
    branch_4 = num(branch_4_index, [1, 2]);
    branch_5 = num(branch_5_index, [1, 2]);

    figure(801)
    loglog(branch_1(:, 1), branch_1(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    hold on
    loglog(branch_2(:, 1), branch_2(:, 2), ':', 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_3(:, 1), branch_3(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_4(:, 1), branch_4(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_5(:, 1), branch_5(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);

    pbaspect([1.2 0.9 1])
    xlabel(txt(1))
    ylabel(txt(2))
    box on
    set(gca,'Fontsize', 16);

end

for i = 1 : 1 : 4
    [num, txt] = xlsread('Ultrasensitivity_k0_PRXSO2H_PRX1 & PRX2 & PRX3', i);
    column_num = length(txt);
    
    % Obtain branch index: 
    branch_1_index = find(num(:, column_num) == 1); 
    branch_2_index = find(num(:, column_num) == 2);
    branch_3_index = find(num(:, column_num) == 3);
    branch_4_index = find(num(:, column_num) == 4);
    branch_5_index = find(num(:, column_num) == 5);

    % Obtain branch data: 
    branch_1 = num(branch_1_index, [1, 2]);
    branch_2 = num(branch_2_index, [1, 2]);
    branch_3 = num(branch_3_index, [1, 2]);
    branch_4 = num(branch_4_index, [1, 2]);
    branch_5 = num(branch_5_index, [1, 2]);

    figure(802)
    loglog(branch_1(:, 1), branch_1(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    hold on
    loglog(branch_2(:, 1), branch_2(:, 2), ':', 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_3(:, 1), branch_3(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_4(:, 1), branch_4(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_5(:, 1), branch_5(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);

    pbaspect([1.2 0.9 1])
    xlabel(txt(1))
    ylabel(txt(2))
    box on
    set(gca,'Fontsize', 16);

end

%

%% --------------- Generate Fig. 8C and 8D ------------------------- %%

for i = 1 : 1 : 4
    [num, txt] = xlsread('Bistability_k0_H2O2_PRX1 & PRX2 & PRX3', i);
    column_num = length(txt);
    
    % Obtain branch index: 
    branch_1_index = find(num(:, column_num) == 1); 
    branch_2_index = find(num(:, column_num) == 2);
    branch_3_index = find(num(:, column_num) == 3);
    branch_4_index = find(num(:, column_num) == 4);
    branch_5_index = find(num(:, column_num) == 5);

    % Obtain branch data: 
    branch_1 = num(branch_1_index, [1, 2]);
    branch_2 = num(branch_2_index, [1, 2]);
    branch_3 = num(branch_3_index, [1, 2]);
    branch_4 = num(branch_4_index, [1, 2]);
    branch_5 = num(branch_5_index, [1, 2]);

    figure(803)
    loglog(branch_1(:, 1), branch_1(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    hold on
    loglog(branch_2(:, 1), branch_2(:, 2), ':', 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_3(:, 1), branch_3(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_4(:, 1), branch_4(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_5(:, 1), branch_5(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);

    pbaspect([1.2 0.9 1])
    xlabel(txt(1))
    ylabel(txt(2))
    box on
    set(gca,'Fontsize', 16);

end

for i = 1 : 1 : 4
    [num, txt] = xlsread('Bistability_k0_PRXSO2H_PRX1 & PRX2 & PRX3', i);
    column_num = length(txt);
    
    % Obtain branch index: 
    branch_1_index = find(num(:, column_num) == 1); 
    branch_2_index = find(num(:, column_num) == 2);
    branch_3_index = find(num(:, column_num) == 3);
    branch_4_index = find(num(:, column_num) == 4);
    branch_5_index = find(num(:, column_num) == 5);

    % Obtain branch data: 
    branch_1 = num(branch_1_index, [1, 2]);
    branch_2 = num(branch_2_index, [1, 2]);
    branch_3 = num(branch_3_index, [1, 2]);
    branch_4 = num(branch_4_index, [1, 2]);
    branch_5 = num(branch_5_index, [1, 2]);

    figure(804)
    loglog(branch_1(:, 1), branch_1(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    hold on
    loglog(branch_2(:, 1), branch_2(:, 2), ':', 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_3(:, 1), branch_3(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_4(:, 1), branch_4(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_5(:, 1), branch_5(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);

    pbaspect([1.2 0.9 1])
    xlabel(txt(1))
    ylabel(txt(2))
    box on
    set(gca,'Fontsize', 16);

end

%

%% --------------- Generate Fig. 8E and 8F ------------------------- %%

for i = [5, 3, 6]
    [num, txt] = xlsread('Ultrasensitivity_k0_H2O2_PRX1 & PRX2 & PRX3', i);
    column_num = length(txt);
    
    % Obtain branch index: 
    branch_1_index = find(num(:, column_num) == 1); 
    branch_2_index = find(num(:, column_num) == 2);
    branch_3_index = find(num(:, column_num) == 3);
    branch_4_index = find(num(:, column_num) == 4);
    branch_5_index = find(num(:, column_num) == 5);

    % Obtain branch data: 
    branch_1 = num(branch_1_index, [1, 2]);
    branch_2 = num(branch_2_index, [1, 2]);
    branch_3 = num(branch_3_index, [1, 2]);
    branch_4 = num(branch_4_index, [1, 2]);
    branch_5 = num(branch_5_index, [1, 2]);

    figure(805)
    loglog(branch_1(:, 1), branch_1(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    hold on
    loglog(branch_2(:, 1), branch_2(:, 2), ':', 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_3(:, 1), branch_3(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_4(:, 1), branch_4(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_5(:, 1), branch_5(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);

    pbaspect([1.2 0.9 1])
    xlabel(txt(1))
    ylabel(txt(2))
    box on
    set(gca,'Fontsize', 16);

end

for i = [5, 3, 6]
    [num, txt] = xlsread('Ultrasensitivity_k0_PRXSO2H_PRX1 & PRX2 & PRX3', i);
    column_num = length(txt);
    
    % Obtain branch index: 
    branch_1_index = find(num(:, column_num) == 1); 
    branch_2_index = find(num(:, column_num) == 2);
    branch_3_index = find(num(:, column_num) == 3);
    branch_4_index = find(num(:, column_num) == 4);
    branch_5_index = find(num(:, column_num) == 5);

    % Obtain branch data: 
    branch_1 = num(branch_1_index, [1, 2]);
    branch_2 = num(branch_2_index, [1, 2]);
    branch_3 = num(branch_3_index, [1, 2]);
    branch_4 = num(branch_4_index, [1, 2]);
    branch_5 = num(branch_5_index, [1, 2]);

    figure(806)
    loglog(branch_1(:, 1), branch_1(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    hold on
    loglog(branch_2(:, 1), branch_2(:, 2), ':', 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_3(:, 1), branch_3(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_4(:, 1), branch_4(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_5(:, 1), branch_5(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);

    pbaspect([1.2 0.9 1])
    xlabel(txt(1))
    ylabel(txt(2))
    box on
    set(gca,'Fontsize', 16);

end

%

%% --------------- Generate Fig. 8G and 8H ------------------------- %%

for i = [5, 3, 6]
    [num, txt] = xlsread('Bistability_k0_H2O2_PRX1 & PRX2 & PRX3', i);
    column_num = length(txt);
    
    % Obtain branch index: 
    branch_1_index = find(num(:, column_num) == 1); 
    branch_2_index = find(num(:, column_num) == 2);
    branch_3_index = find(num(:, column_num) == 3);
    branch_4_index = find(num(:, column_num) == 4);
    branch_5_index = find(num(:, column_num) == 5);

    % Obtain branch data: 
    branch_1 = num(branch_1_index, [1, 2]);
    branch_2 = num(branch_2_index, [1, 2]);
    branch_3 = num(branch_3_index, [1, 2]);
    branch_4 = num(branch_4_index, [1, 2]);
    branch_5 = num(branch_5_index, [1, 2]);

    figure(807)
    loglog(branch_1(:, 1), branch_1(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    hold on
    loglog(branch_2(:, 1), branch_2(:, 2), ':', 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_3(:, 1), branch_3(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_4(:, 1), branch_4(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_5(:, 1), branch_5(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);

    pbaspect([1.2 0.9 1])
    xlabel(txt(1))
    ylabel(txt(2))
    box on
    set(gca,'Fontsize', 16);

end

for i = [5, 3, 6]
    [num, txt] = xlsread('Bistability_k0_PRXSO2H_PRX1 & PRX2 & PRX3', i);
    column_num = length(txt);
    
    % Obtain branch index: 
    branch_1_index = find(num(:, column_num) == 1); 
    branch_2_index = find(num(:, column_num) == 2);
    branch_3_index = find(num(:, column_num) == 3);
    branch_4_index = find(num(:, column_num) == 4);
    branch_5_index = find(num(:, column_num) == 5);

    % Obtain branch data: 
    branch_1 = num(branch_1_index, [1, 2]);
    branch_2 = num(branch_2_index, [1, 2]);
    branch_3 = num(branch_3_index, [1, 2]);
    branch_4 = num(branch_4_index, [1, 2]);
    branch_5 = num(branch_5_index, [1, 2]);

    figure(808)
    loglog(branch_1(:, 1), branch_1(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    hold on
    loglog(branch_2(:, 1), branch_2(:, 2), ':', 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_3(:, 1), branch_3(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_4(:, 1), branch_4(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);
    loglog(branch_5(:, 1), branch_5(:, 2), 'Color', [51,102,255]/255, 'LineWidth', 3);

    pbaspect([1.2 0.9 1])
    xlabel(txt(1))
    ylabel(txt(2))
    box on
    set(gca,'Fontsize', 16);

end

%
