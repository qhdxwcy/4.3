phi0 = 0.8;
R_all = 2.5;
R_54 = 2.5^(1/3);
R_43 = 2.5^(1/3);
R_32 = 2.5^(1/3);
phi0_values = [0.5,0.8,1.1,1.4,1.7,2.0];

C_T = phi0*(1-1/R_all);

%% alpha3
%alpha2 = -1*pi/4:5*pi/180:pi/4;
%tan(alpha3)*R_32/cos(alpha3)-tan(alpha2)/cos(alpha2)=tan(alpha4)*R_43*R_32/cos(alpha4);
[alpha2_grid, alpha4_grid] = meshgrid(linspace(-pi/3, pi/3, 200));
% 初始化 alpha3 的结果矩阵
alpha3_grid = zeros(size(alpha2_grid));

for i = 1:size(alpha2_grid, 1)
    for j = 1:size(alpha2_grid, 2)
        a2 = alpha2_grid(i, j);
        a4 = alpha4_grid(i, j);

        fun = @(a3) tan(a3) * R_32 / cos(a3) - tan(a2) / cos(a2) - tan(a4) * R_43 * R_32 / cos(a4);
        a3_guess = 0;
        % 调用 fzero 求解
        try
            alpha3_grid(i, j) = fzero(fun, a3_guess);
        catch
            alpha3_grid(i, j) = NaN;  % 无解时标记为 NaN
        end
    end
end
figure;
contourf(rad2deg(alpha2_grid), rad2deg(alpha4_grid), rad2deg(alpha3_grid), 50,'LineStyle', ':');
xlabel('\alpha_2 (deg)');
ylabel('\alpha_4 (deg)');
title('\alpha_3 云图');
colorbar;

%% beta2
alpha2 = linspace(-pi/3,pi/3,200);
for k = 1:length(phi0_values)
    beta2 = atan(R_54*R_43*R_32/phi0_values(k)+tan(alpha2));
    alpha2_c = alpha2*180/pi;
    beta2_c = beta2*180/pi;
    subplot(2,3,k);
    plot(alpha2_c,beta2_c,'LineWidth',2);
    xlabel('\alpha2','FontSize', 12);
    ylabel('\beta2','FontSize', 12);
    title(['\beta2 可能取值曲线 (\phi_0 = ', num2str(phi0_values(k)), ')'],'FontSize', 14);
    grid on;
end

%% beta3
beta3_grid = zeros(size(alpha2_grid));
for k = 1:length(phi0_values)
    for i = 1:size(alpha2_grid, 1)
        for j = 1:size(alpha2_grid, 2)
            beta3_grid(i,j) = atan(R_54*R_43/phi0_values(k)+tan(alpha3_grid(i,j)));
        end
    end
    subplot(2,3,k);
    contourf(rad2deg(alpha2_grid), rad2deg(alpha4_grid), rad2deg(beta3_grid), 50,'LineStyle', ':');
    xlabel('\alpha_2 (deg)');
    ylabel('\alpha_4 (deg)');
    title(['\beta3 云图 (\phi_0 = ', num2str(phi0_values(k)), ')'],'FontSize', 14);
    colorbar;
end

%% beta4
alpha4 = linspace(-pi/3,pi/3,200);
for k = 1:length(phi0_values)
    beta4 = atan(R_54/phi0_values(k)+tan(alpha4));
    subplot(2,3,k);
    plot(rad2deg(alpha4),rad2deg(beta4),'LineWidth',2);
    xlabel('\alpha4(deg)');
    ylabel('\beta4(deg)');
    title(['\beta4 可能取值曲线 (\phi_0 = ', num2str(phi0_values(k)), ')'],'FontSize', 14);
    grid on;
end

%% beta5
beta5 = zeros(1,length(phi0_values));
for k = 1:length(phi0_values)
    beta5(k) = atan(1/phi0_values(k));
end
scatter(phi0_values,rad2deg(beta5),'filled','MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black');
title('\beta5在不同\phi_0下的取值','FontSize', 14);
xlabel('\phi_0');
ylabel('\beta5(deg)');


%% W3/W2
W3W2_grid = zeros(size(alpha2_grid));
beta2_phi = zeros(length(phi0_values),size(alpha2_grid,1));
beta3_phi = zeros(length(phi0_values),size(alpha2_grid,1),size(alpha2_grid,1));
for k = 1:length(phi0_values)
    beta2_phi(k,:) = atan(R_54*R_43*R_32/phi0_values(k)+tan(alpha2));
end
for k = 1:length(phi0_values)
    beta3_phi(k,:,:) = atan(R_54*R_43/phi0_values(k)+tan(alpha3_grid));
end

for k = 1:length(phi0_values)
    for i = 1:size(alpha2_grid, 1)
        for j = 1:size(alpha2_grid, 2)
            W3W2_grid(i,j) = cos(beta2_phi(k,i))*R_32/cos(beta3_phi(k,i,j));
        end
    end
    mask1 = (W3W2_grid <= 1);
    
    subplot(2,3,k);
    contourf(rad2deg(alpha2_grid), rad2deg(alpha4_grid), W3W2_grid, 50,'LineStyle', ':');
    hold on;
    contourf(rad2deg(alpha2_grid), rad2deg(alpha4_grid), mask1, [0.5, 0.5], ...
        'LineColor', 'none', 'FaceColor', 'w', 'FaceAlpha', 0.6);
    [C, h] = contour(rad2deg(alpha2_grid), rad2deg(alpha4_grid), W3W2_grid, [1,1], ...
        'LineColor', 'r', 'LineWidth', 2, 'ShowText', 'on');
    xlabel('\alpha_2 (deg)','FontSize', 12);
    ylabel('\alpha_4 (deg)','FontSize', 12);
    title(['W3/W2云图 (\phi_0 = ', num2str(phi0_values(k)), ')'],'FontSize', 14);
    colorbar;
    hold off;
end


%% V4/V3
V4V3_grid = zeros(size(alpha2_grid));
for k = 1:length(phi0_values)
    for i = 1:size(alpha2_grid, 1)
        for j = 1:size(alpha2_grid, 2)
            V4V3_grid(i,j) = cos(alpha3_grid(i,j))*R_43/cos(alpha4_grid(i,j));
        end
    end
    mask2 = (V4V3_grid <= 1);
    
    subplot(2,3,k);
    contourf(rad2deg(alpha2_grid), rad2deg(alpha4_grid), V4V3_grid, 50,'LineStyle', ':');
    hold on;
    contourf(rad2deg(alpha2_grid), rad2deg(alpha4_grid), mask2, [0.5, 0.5], ...
        'LineColor', 'none', 'FaceColor', 'w', 'FaceAlpha', 0.6);
    [C, h] = contour(rad2deg(alpha2_grid), rad2deg(alpha4_grid), V4V3_grid, [1,1], ...
        'LineColor', 'r', 'LineWidth', 2, 'ShowText', 'on');
    clabel(C, h, 'FontSize', 12, 'Color', 'r');
    xlabel('\alpha_2 (deg)','FontSize', 12);
    ylabel('\alpha_4 (deg)','FontSize', 12);
    title(['V4/V3云图 (\phi_0 = ', num2str(phi0_values(k)), ')'],'FontSize', 14);
    colorbar;
    hold off;
end

%% W5/W4
W5W4 = zeros(1,size(alpha2_grid,2));
beta4_phi = zeros(length(phi0_values),size(alpha2_grid,2));
for k = 1:length(phi0_values)
    beta4_phi(k,:) = atan(R_54/phi0_values(k)+tan(alpha4));
end
for k = 1:length(phi0_values)
    for i = 1:size(alpha2_grid, 2)
        W5W4(i) = R_54*cos(beta4_phi(k,i))/cos(beta5(k));
    end
    subplot(2, 3, k);
    plot(rad2deg(alpha4), W5W4, 'LineWidth', 2);
    hold on;
    idx = find(W5W4 < 1, 1, 'first'); % 返回第一个满足条件的索引
    alpha4_critical = rad2deg(alpha4(idx));   % 对应的 α4 值
    % 绘制临界点的竖线
    xline(alpha4_critical, 'r--', 'LineWidth', 1.5, ...
    'Label', sprintf('α4 = %.2f°', alpha4_critical)); 
    xlabel('\alpha4','FontSize', 12);
    ylabel('W5/W4','FontSize', 12);
    title(['W5/W4可能取值曲线 (\phi_0 = ', num2str(phi0_values(k)), ')'],'FontSize', 14);
    hold off;
end

%% alpha3-alpha4
s1_grid = zeros(size(alpha2_grid));
for i = 1:size(alpha2_grid, 1)
    for j = 1:size(alpha2_grid, 2)
        s1_grid(i,j) = alpha3_grid(i,j)-alpha4(1,j);
    end
end
mask3 = (rad2deg(s1_grid) <= -60 | rad2deg(s1_grid) >= -10 & rad2deg(s1_grid) <= 10 | rad2deg(s1_grid) >= 60);

figure;
contourf(rad2deg(alpha2_grid), rad2deg(alpha4_grid), rad2deg(s1_grid), 50,'LineStyle', ':');
hold on;
contourf(rad2deg(alpha2_grid), rad2deg(alpha4_grid), mask3, [0.5, 0.5], ...
    'LineColor', 'none', 'FaceColor', 'w', 'FaceAlpha', 0.6);
[C, h] = contour(rad2deg(alpha2_grid), rad2deg(alpha4_grid), rad2deg(s1_grid), [-60,-10], ...
    'LineColor', 'r', 'LineWidth', 2, 'ShowText', 'on');
clabel(C, h, 'FontSize', 12, 'Color', 'r');
[D, m] = contour(rad2deg(alpha2_grid), rad2deg(alpha4_grid), rad2deg(s1_grid), [10,60], ...
    'LineColor', 'r', 'LineWidth', 2, 'ShowText', 'on');
clabel(D, m, 'FontSize', 12, 'Color', 'r');
xlabel('\alpha_2 (deg)');
ylabel('\alpha_4 (deg)');
title('\alpha3-\alpha4 云图');
colorbar;
hold off;

%% beta2-beta3
s2_grid = zeros(size(alpha2_grid));
for k = 1:length(phi0_values)
    for i = 1:size(alpha2_grid, 1)
        for j = 1:size(alpha2_grid, 2)
            s2_grid(i,j) = beta2_phi(k,i)-beta3_phi(k,i,j);
        end
    end
    mask4 = (rad2deg(s2_grid) <= -10 | rad2deg(s2_grid) >= 10);
    
    subplot(2,3,k);
    contourf(rad2deg(alpha2_grid), rad2deg(alpha4_grid), rad2deg(s2_grid), 50,'LineStyle', ':');
    hold on;
    contourf(rad2deg(alpha2_grid), rad2deg(alpha4_grid), mask4, [0.5, 0.5], ...
        'LineColor', 'none', 'FaceColor', 'w', 'FaceAlpha', 0.6);
    [C, h] = contour(rad2deg(alpha2_grid), rad2deg(alpha4_grid), rad2deg(s2_grid), [-10,10], ...
        'LineColor', 'r', 'LineWidth', 2, 'ShowText', 'on');
    xlabel('\alpha_2 (deg)','FontSize', 12);
    ylabel('\alpha_4 (deg)','FontSize', 12);
    title(['\beta2-\beta3云图 (\phi_0 = ', num2str(phi0_values(k)), ')'],'FontSize', 14);
    colorbar;
    hold off;
end

%% beta4-beta5
s3 = zeros(1,size(alpha2_grid,2));
for k = 1:length(phi0_values)
    for i = 1:size(alpha2_grid, 2)
        s3(i) = beta4_phi(k,i)-beta5(k);
    end
    subplot(2,3,k);
    plot(rad2deg(alpha4),rad2deg(s3),'LineWidth',2);
    hold on;
    % 找出 s3 < -10 的第一个临界点（确保只取第一个）
    idx_low = find(rad2deg(s3) < -10, 1, 'last');
    if ~isempty(idx_low)
        alpha4_low = rad2deg(alpha4(idx_low));
        xline(alpha4_low, 'r--', 'LineWidth', 1.5, ...
              'Label', sprintf('α4=%.2f°', alpha4_low));
    end
    
    % 找出 s3 > 10 的第一个临界点（确保只取第一个）
    idx_high = find(rad2deg(s3) > 10, 1, 'first');
    if ~isempty(idx_high)
        alpha4_high = rad2deg(alpha4(idx_high));
        xline(alpha4_high, 'r--', 'LineWidth', 1.5, ...
              'Label', sprintf('α4=%.2f°', alpha4_high));
    end
    xlabel('\alpha4','FontSize', 12);
    ylabel('\beta4-\beta5','FontSize', 12);
    title(['\beta4-\beta5可能取值曲线 (\phi_0 = ', num2str(phi0_values(k)), ')'],'FontSize', 14);
end

%% psi
psi = zeros(1,size(alpha2_grid,2));
for k = 1:length(phi0_values)
    for i = 1:size(alpha2_grid, 2)
        psi(i) = phi0_values(k)*tan(alpha4(i))/(R_54*cos(alpha4(i)));
    end
    subplot(2,3,k);
    plot(rad2deg(alpha4),psi,'LineWidth',2);
    xlabel('\alpha4','FontSize', 12);
    ylabel('\psi','FontSize', 12);
    title(['\psi可能取值曲线 (\phi_0 = ', num2str(phi0_values(k)), ')'],'FontSize', 14);
end

%% eta
eta = zeros(1,size(alpha2_grid,2));
for i = 1:size(alpha2_grid, 2)
    eta(i) = (R_all-1)*R_54*cos(alpha4(i))/(2*tan(alpha4(i))*R_all);
end
exclude_range = [-0.1, 0.1];
% 找出不在排除范围内的索引
idx1 = alpha4 < exclude_range(1);  % alpha4 < -0.1
idx2 = alpha4 > exclude_range(2);  % alpha4 > 0.1
plot(rad2deg(alpha4(idx1)),eta(idx1),'b-','LineWidth',2);
hold on;
plot(rad2deg(alpha4(idx2)),eta(idx2),'b-','LineWidth',2);
hold off;
xlabel('\alpha4');
ylabel('\eta');
title('\eta可能取值曲线');


%% 总约束
ss = zeros(size(alpha2_grid));
for k = 1:length(phi0_values)
    for i = 1:size(alpha2_grid, 1)
        for j = 1:size(alpha2_grid, 2)
            W3W2_grid(i,j) = cos(beta2_phi(k,i))*R_32/cos(beta3_phi(k,i,j));
            V4V3_grid(i,j) = cos(alpha3_grid(i,j))*R_43/cos(alpha4_grid(i,j));
            s1_grid(i,j) = alpha3_grid(i,j)-alpha4(1,j);
            s2_grid(i,j) = beta2_phi(k,i)-beta3_phi(k,i,j);
            ss(i,j) = i+j;
        end
    end
    W5W4 = R_54*cos(beta4_phi(k,:))./cos(beta5(k));
    alpha4_critical = rad2deg(alpha4(find(W5W4 < 1, 1, 'first')));
    mask5 = (rad2deg(alpha4_grid) >= alpha4_critical); % W5/W4<1的约束区域
    s3 = beta4_phi(k,:) - beta5(k);
    s3_deg = rad2deg(s3);
    idx_low = find(s3_deg < -10, 1, 'last');
    idx_high = find(s3_deg > 10, 1, 'first');
    alpha4_low = rad2deg(alpha4(idx_low));
    alpha4_high = rad2deg(alpha4(idx_high));
    mask6 = (rad2deg(alpha4_grid) <= alpha4_low) | (rad2deg(alpha4_grid) >= alpha4_high);

    mask1 = (W3W2_grid <= 1);
    mask2 = (V4V3_grid <= 1);
    mask3 = (rad2deg(s1_grid) <= -60 | rad2deg(s1_grid) >= -10 & rad2deg(s1_grid) <= 10 | rad2deg(s1_grid) >= 60);
    mask4 = (rad2deg(s2_grid) <= -10 | rad2deg(s2_grid) >= 10);
    mask_total = mask1 | mask2 | mask3 | mask4 | mask5 | mask6;
    
    subplot(2,3,k);
    contourf(rad2deg(alpha2_grid), rad2deg(alpha4_grid), rad2deg(s2_grid), 50,'LineStyle', 'none');
    hold on;
    contourf(rad2deg(alpha2_grid), rad2deg(alpha4_grid), mask_total, [0.5, 0.5], ...
        'LineColor', 'none', 'FaceColor', 'w', 'FaceAlpha', 0.6);
    [C_total, h_total] = contour(rad2deg(alpha2_grid), rad2deg(alpha4_grid), double(mask_total), [0.5, 0.5], ...
        'LineColor', 'r', 'LineWidth', 2);
    xlabel('\alpha_2 (deg)','FontSize', 12);
    ylabel('\alpha_4 (deg)','FontSize', 12);
    title(['综合约束(\phi_0 = ', num2str(phi0_values(k)), ')'],'FontSize', 14);
    colorbar;
    hold off;
end