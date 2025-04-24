clc
clear
close all

rng(100, 'twister');
alpha = 0.5;
epsilon = 0.05;
num_N = 5000;
dt = 0.0001;
T = dt * num_N;
c_alpha = alpha * gamma((1+alpha)/2) / (2^(1-alpha) * pi^(0.5) * gamma(1-alpha/2));
Lt1 = generate_Lt(alpha,epsilon,num_N,dt,c_alpha);
Lt2 = generate_Lt(alpha,epsilon,num_N,dt,c_alpha);

Xzero = 5;
X_simul = zeros(num_N + 1,1);
Yzero = -5;
Y_simul = zeros(num_N + 1,1);
X_simul(1) = Xzero;
Y_simul(1) = Yzero;
for i = 2:num_N + 1
        X_simul(i) = X_simul(i-1) + ( 5*X_simul(i-1) - Y_simul(i-1)^2) * dt + sqrt(1) * (Lt1(i) - Lt1(i-1));
        Y_simul(i) = Y_simul(i-1) + ( 5*X_simul(i-1) + Y_simul(i-1)) * dt + sqrt(1) * (Lt2(i) - Lt2(i-1));
end

any(isinf(X_simul))   % 检查是否有Inf



%%
X = X_simul(1:end-1);
Y = Y_simul(1:end-1);
FX_X =(X_simul(2:end)-X_simul(1:end-1))./dt; 
FX_Y =(Y_simul(2:end)-Y_simul(1:end-1))./dt; 

order = 4;  % 最大幂次
[deg_funcs, deg_strs] = deg_func_xy(order);
num_funcs = length(deg_funcs);
Phi = zeros(length(X), num_funcs); 
for k = 1:num_funcs
    func = deg_funcs{k}; % 获取第 k 个函数句柄
    Phi(:, k) = func(X , Y); % 计算函数值并存储
end
Phi_old = Phi;

[Drift1_X,Drift_var1_X] = MSBL(Phi, FX_X, 1, 1);
[Drift1_Y,Drift_var1_Y] = MSBL(Phi, FX_Y, 1, 1);
p = 1;
[~,Drift2_X,Drift_var2_X] = StepwiseRSBL(Phi, FX_X, 1, 1, p);
[~,Drift2_Y,Drift_var2_Y] = StepwiseRSBL(Phi, FX_Y, 1, 1, p);

var_alpha = 2 * dt * c_alpha /(2 - alpha);
XG =  X_simul(1:end-1);
YG = Y_simul(1:end-1);
GX1_X = ((X_simul(2:end,:)- X_simul(1:end-1,:)-Phi_old*Drift1_X*dt)).^2/ var_alpha;
GY1_Y = ((Y_simul(2:end,:) - Y_simul(1:end-1,:) - Phi_old * Drift1_Y * dt)).^2/ var_alpha;

GX2_X = ((X_simul(2:end,:)- X_simul(1:end-1,:)-Phi_old*Drift2_X*dt)).^2/ var_alpha;
GY2_Y = ((Y_simul(2:end,:) - Y_simul(1:end-1,:) - Phi_old*Drift2_Y*dt)).^2/ var_alpha;

indices_GX_X2 = find(abs(GX2_X) >300);
indices_GY_Y2 = find(abs(GY2_Y) > 1000);
total_indices2 = unique([indices_GX_X2; indices_GY_Y2]);

GX1_X(total_indices2) = [];
GY1_Y(total_indices2) = [];

GX2_X(total_indices2) = [];
GY2_Y(total_indices2) = [];

XG(total_indices2) = [];
YG(total_indices2) = [];

Phi_old(total_indices2, :) = []; 

[Diff1_X,Diff_var1_X] = MSBL(Phi_old, GX1_X*100, 1, 1);
[Diff1_Y,Diff_var1_Y] = MSBL(Phi_old, GY1_Y*100, 1, 1);
p2 = 0.01;
[~,Diff2_X,Diff_var2_X] = StepwiseRSBL(Phi_old, GX2_X*100, 1, 1, p2);
[~,Diff2_Y,Diff_var2_Y] = StepwiseRSBL(Phi_old, GY2_Y*100, 1, 1, p2);
%%
Drift1_true = Drift2_X;
Drift2_true = Drift2_Y;

%%
MSE1 = mse(Drift1_true - Drift1_X);
MSE2 = mse(Drift1_true - Drift2_X);
MSE3 = mse(Drift2_true - Drift1_Y);
MSE4 = mse(Drift2_true - Drift2_Y);

disp(MSE1);
disp(MSE2);
disp(MSE3);
disp(MSE4);

MAE1 = abs(Drift1_true - Drift1_X);
MAE2 = abs(Drift1_true - Drift2_X);
MAE3 = abs(Drift2_true - Drift1_Y);
MAE4 = abs(Drift2_true - Drift2_Y);
%%
Diff1_true = Diff2_X;
Diff2_true = Diff2_Y;
%%
MSE5 = mse(Diff1_true - Diff1_X);
MSE6 = mse(Diff1_true - Diff2_X);
MSE7 = mse(Diff2_true - Diff1_Y);
MSE8 = mse(Diff2_true - Diff2_Y);

disp(MSE5);
disp(MSE6);
disp(MSE7);
disp(MSE8);

MAE5 = abs(Diff1_true - Diff1_X);
MAE6 = abs(Diff1_true - Diff2_X);
MAE7 = abs(Diff2_true - Diff1_Y);
MAE8 = abs(Diff2_true - Diff2_Y);

%% 画样本路径图
map = colormap(nclCM(232)); % color包里选颜色
map = flipud(map);
close all;

figureUnits = 'centimeters';
figureWidth = 28; % 宽度增加以适应两个子图
figureHeight = 10;

% 窗口设置
figureHandle = figure;
set(gcf, 'Units', figureUnits, 'Position', [2 5 figureWidth figureHeight]);
T = 0:dt:dt * num_N;
subplotGap = 0.08; 
subplotWidth = (1 - 2.8* subplotGap) / 2; 
subplotHeight = 0.7;  % 降低子图高度
subplotBottom = 0.25;  % 提升子图起始位置



% ====================== 第一个子图 ======================
subplot('Position', [subplotGap, subplotBottom, subplotWidth, subplotHeight]);
hold on;
scatter(T, X_simul, 5, X_simul, 'filled');
hXLabel1 = xlabel('$T$', 'Interpreter', 'latex');
hYLabel1 = ylabel('$X_t$', 'Interpreter', 'latex');

colormap(map);
colorbar;
set(gca, 'Box', 'off', ...
         'LineWidth', 1, ...
         'XGrid', 'off', 'YGrid', 'off', ...
         'TickDir', 'out', 'TickLength', [.005 .005], ...
         'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', 0:0.2:0.5, ...
         'XLim', [0 0.51], ...
         'YTick', -20:10:15, ...
         'YLim', [-20 15]);

set(gca, 'FontName', 'Times New Roman', 'FontSize', 15,'FontWeight' , 'normal')
set(hXLabel1, 'FontSize', 17, 'FontName', 'Times New Roman','FontWeight' , 'normal')
set(hYLabel1, 'FontSize', 17, 'FontName', 'Times New Roman','FontWeight' , 'normal')

% ====================== 第二个子图 ======================
subplot('Position', [2*subplotGap + subplotWidth, subplotBottom, subplotWidth, subplotHeight]);
hold on;

scatter(T, Y_simul, 5, Y_simul, 'filled');
hXLabel2 = xlabel('$T$', 'Interpreter', 'latex');
hYLabel2= ylabel('$Y_t$', 'Interpreter', 'latex');

colormap(map);
colorbar;
set(gca, 'Box', 'off', ...
         'LineWidth', 1, ...
         'XGrid', 'off', 'YGrid', 'off', ...
         'TickDir', 'out', 'TickLength', [.005 .005], ...
         'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', 0:0.2:0.5, ...
         'XLim', [0 0.51], ...
         'YTick', -7:5:18, ...
         'YLim', [-7 18]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15,'FontWeight' , 'normal')
set(hXLabel2, 'FontSize', 17, 'FontName', 'Times New Roman','FontWeight' , 'normal')
set(hYLabel2, 'FontSize', 17, 'FontName', 'Times New Roman','FontWeight' , 'normal')
% 背景颜色
set(gcf, 'Color', [1 1 1]);
%% 图片输出
figW = figureWidth;%将图形的宽度和高度赋值给变量 figW 和 figH
figH = figureHeight;
set(figureHandle,'PaperUnits',figureUnits);%设置图形的打印单位为 figureUnits
set(figureHandle,'PaperPosition',[0 0 figW figH]);
fileout = '样本路径图';
print(figureHandle,[fileout,'.png'],'-r600','-dpng');

%% 画漂移项和耗散项的分布图
C1 = addcolorplus(252); 
C2 = addcolorplus(119); 
close all
data = struct();
data(1).mu1 = Drift1_X(2);  
data(1).sigma1 = sqrt(Drift_var1_X(2)); 
data(1).mu2 = Drift2_X(2); 
data(1).sigma2 = sqrt(Drift_var2_X(2));

data(2).mu1 = Drift1_X(7);  
data(2).sigma1 = sqrt(Drift_var1_X(7)); 
data(2).mu2 = Drift2_X(7); 
data(2).sigma2 = sqrt(Drift_var2_X(7));

data(3).mu1 = Diff1_X(1);  
data(3).sigma1 = sqrt(Diff_var1_X(1)); 
data(3).mu2 = Diff2_X(1); 
data(3).sigma2 = sqrt(Diff_var2_X(1));

data(4).mu1 = Drift1_Y(2);  
data(4).sigma1 = sqrt(Drift_var1_Y(2)); 
data(4).mu2 = Drift2_Y(2); 
data(4).sigma2 = sqrt(Drift_var2_Y(2));

data(5).mu1 = Drift1_Y(6);  
data(5).sigma1 = sqrt(Drift_var1_Y(6)); 
data(5).mu2 = Drift2_Y(6); 
data(5).sigma2 = sqrt(Drift_var2_Y(6));

data(6).mu1 = Diff1_Y(1);  
data(6).sigma1 = sqrt(Diff_var1_Y(1)); 
data(6).mu2 = Diff2_Y(1); 
data(6).sigma2 = sqrt(Diff_var2_Y(1));

% 图片尺寸设置
figureUnits = 'centimeters';
figureWidth = 28; % 宽度增加以适应 3 列
figureHeight = 18; % 高度增加以适应 2 行

% 创建画布
figureHandle = figure;
set(gcf, 'Units', figureUnits, 'Position', [2 5 figureWidth figureHeight]);

topMargin = 0.04;     % 上边距减小
bottomMargin = 0.15;  % 下边距保持较大（用于标题）
leftMargin = 0.08;
rightMargin = 0.04;

% 间隙参数
subplotGapX = 0.07;  
subplotGapY = 0.15;

% 计算子图尺寸
subplotWidth = (1 - leftMargin - rightMargin - 2*subplotGapX)/3;
subplotHeight = (1 - topMargin - bottomMargin - subplotGapY)/2;
% 循环绘制 6 张图
for i = 1:6
    row = ceil(i/3);          % 行号（1或2）
    col = mod(i-1,3)+1;       % 列号（1-3）
    
    left = leftMargin + (col-1)*(subplotWidth + subplotGapX);
    bottom = 1 - topMargin - row*subplotHeight - (row-1)*subplotGapY;
    % 创建子图（使用绝对定位）
    subplot('Position', [left, bottom, subplotWidth, subplotHeight]); 
    hold on;

    x1 = linspace(data(i).mu1 - 4 * data(i).sigma1, data(i).mu1 + 4 * data(i).sigma1, 1000); 
    pdf1 = normpdf(x1, data(i).mu1, data(i).sigma1); 
    x2 = linspace(data(i).mu2 - 4 * data(i).sigma2, data(i).mu2 + 4 * data(i).sigma2, 1000); 
    pdf2 = normpdf(x2, data(i).mu2, data(i).sigma2); 
    
    % 绘制分布
    area(x1, pdf1, 'LineWidth', 1.5, 'FaceColor', C1, 'EdgeColor', C1, ...
        'FaceAlpha', .3, 'EdgeAlpha', 1);
    area(x2, pdf2, 'LineWidth', 1.5, 'FaceColor', C2, 'EdgeColor', C2, ...
        'FaceAlpha', .3, 'EdgeAlpha', 1);
    
    % 绘制均值虚线
    line([data(i).mu1 data(i).mu1], [0 normpdf(data(i).mu1, data(i).mu1, data(i).sigma1)], ...
        'LineStyle', '--', 'Color', C1, 'LineWidth', 1.5);
    line([data(i).mu2 data(i).mu2], [0 normpdf(data(i).mu2, data(i).mu2, data(i).sigma2)], ...
        'LineStyle', '--', 'Color', C2, 'LineWidth', 1.5);

    set(gca, 'Box', 'off', 'XGrid', 'off', 'YGrid', 'off', ...
             'TickDir', 'out', 'TickLength', [.01 .01], ...
             'XMinorTick', 'off', 'YMinorTick', 'off', ...
             'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
             'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'normal');
    
    xlabel(sprintf('$\\theta_{%d}$', i), ...
           'FontSize', 17, ...
           'FontName', 'Times New Roman', ...
           'FontWeight', 'normal', ...
           'Interpreter', 'latex'); 
    if i == 1
        ylabel('PDF', 'FontSize', 15, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
    end
    if i == 4
        ylabel('PDF', 'FontSize', 15, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
    end

    if i ==1
        legend({'SBL', 'AT-RSBL'}, 'Location', 'northwest', 'FontSize', 12, 'Box', 'off', 'FontWeight', 'normal');
    end

end

% 背景颜色
set(gcf, 'Color', [1 1 1]);
%% 图片输出
figW = figureWidth;
figH = figureHeight;
set(figureHandle,'PaperUnits',figureUnits);
set(figureHandle,'PaperPosition',[0 0 figW figH]);
fileout = '漂移和耗散项的分布';
print(figureHandle,[fileout,'.png'],'-r600','-dpng');

%% 画漂移函数的三维图像
close all
% 定义函数
fun = @(a,b) 5*a - b.^2;
fun1 = @(a,b) Drift2_X(2)*a + Drift2_X(7)*b.^2;
gun = @(a,b) 5*a + b;
gun1 = @(a,b) Drift2_Y(2)*a + Drift2_Y(6)*b;

map = colormap(nclCM(232)); 
map = flipud(map);

% 图形参数设置
figureUnits = 'centimeters';
figureWidth = 30; 
figureHeight = 18;
figureHandle = figure;
set(gcf, 'Units', figureUnits, 'Position', [2 5 figureWidth figureHeight]);

% 子图布局参数
topMargin = 0.04;      % 上边距
bottomMargin = 0.15;   % 下边距（为标题留空间）
leftMargin = 0.06;
rightMargin = 0.04;
subplotGapX = 0.08;    % 横向间距
subplotGapY = 0.15;     % 纵向间距

% 计算子图尺寸
subplotWidth = (1 - leftMargin - rightMargin - 2*subplotGapX)/3;
subplotHeight = (1 - topMargin - bottomMargin - subplotGapY)/2;

% 循环创建子图
for i = 1:6
    % 计算子图位置
    row = ceil(i/3);
    col = mod(i-1,3)+1;
    
    left = leftMargin + (col-1)*(subplotWidth + subplotGapX);
    bottom = 1 - topMargin - row*subplotHeight - (row-1)*subplotGapY;
    
    % 创建子图
    subplot('Position', [left, bottom, subplotWidth, subplotHeight]);
    hold on;
    
    % 根据子图编号选择绘图函数
    switch i
        case 1
            fmesh(fun, 'LineWidth', 1, 'MeshDensity', 25);
            zlabel('$f_1(x,y)$', 'Interpreter', 'latex');
        case 2
            fmesh(fun1, 'LineWidth', 1, 'MeshDensity', 25);
            zlabel('$\hat{f}_1(x,y)$', 'Interpreter', 'latex');
        case 3
            fun_error = @(x,y) abs(fun(x,y) - fun1(x,y));
            fmesh(fun_error, 'LineWidth', 1, 'MeshDensity', 25);
            zlabel('$|f_1 -\hat{f}_1|$', 'Interpreter', 'latex');
            set(gca, 'ZTick', 0:0.02:0.04, 'ZLim', [0 0.04]);
        case 4
            fmesh(gun, 'LineWidth', 1, 'MeshDensity', 25);
            zlabel('$f_2(x,y)$', 'Interpreter', 'latex');
        case 5
            fmesh(gun1, 'LineWidth', 1, 'MeshDensity', 25);
            zlabel('$\hat{f}_2(x,y)$', 'Interpreter', 'latex');
        case 6
            gun_error = @(x,y) abs(gun(x,y) - gun1(x,y));
            fmesh(gun_error, 'LineWidth', 1, 'MeshDensity', 25);
            zlabel('$|f_2 -\hat{f}_2|$', 'Interpreter', 'latex');
            set(gca, 'ZTick', 0:0.04:0.08, 'ZLim', [0 0.08]);
    end
    
    % 公共图形设置
    view(-40, 12);
    colormap(map);
    colorbar;
    
    % 坐标轴标签
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$y$', 'Interpreter', 'latex');
    
    % 图形样式设置
    set(gca, 'Box', 'off', ...
             'LineWidth', 1, 'GridLineStyle', '-', ...
             'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', ...
             'TickDir', 'out', 'TickLength', [.015 .015], ...
             'XMinorTick', 'off', 'YMinorTick', 'off', 'ZMinorTick', 'off', ...
             'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], 'ZColor', [.1 .1 .1]);
    
    % 字体设置
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
    set(findobj(gca, 'Type', 'text'), 'FontName', 'Times New Roman', 'FontSize', 17);
   
end

% 全局设置
set(gcf, 'Color', [1 1 1]);

%% 图片输出
figW = figureWidth;
figH = figureHeight;
set(figureHandle,'PaperUnits',figureUnits);
set(figureHandle,'PaperPosition',[0 0 figW figH]);
fileout = '漂移项三维误差图像';
print(figureHandle,[fileout,'.png'],'-r900','-dpng');
%%
function Lt = generate_Lt(alpha,epsilon,num_N,dt,c_alpha)
    sigma = sqrt(2 * c_alpha * epsilon^(2-alpha) / (2-alpha));
    FUN1 = @(y) (1 - (-y) .^ (-alpha)) / (2 * (1 - epsilon^(-alpha)));
    FUN2 = @(y) (1 + y .^ (-alpha) - 2 * epsilon^(-alpha)) / (2 - 2 * epsilon^(-alpha));
    F1_inv = @(u) - (1 - 2 * u * (1 - epsilon^(-alpha)))^(-1 / alpha);
    F2_inv = @(u) ((2 * u - 1) + epsilon^(-alpha) * (2 - 2 * u))^(-1 / alpha);
    
    dB = sqrt(dt) * randn(1, num_N);  
    Wt = [0 cumsum(dB)];   
    
    lambda = 2 * c_alpha * (epsilon^(-alpha) - 1) / alpha;
    pois = poissrnd(lambda * dt,1, num_N);
    Y_value = compute_Z_at_time(FUN1,FUN2,F1_inv,F2_inv,epsilon, pois);
    N_epsilon = [0 cumsum(Y_value)];
   
    Lt = sigma * Wt + N_epsilon;
end

% 逆变换抽样
function y_samples = inverse_transform_sampling(FUN1,FUN2,F1_inv,F2_inv,epsilon, N)
    y_samples = zeros(1, N);
    u = rand(1, N);
    for i = 1:N
        if u(i) < FUN1(-epsilon)
            y_samples(i) = F1_inv(u(i));
        elseif u(i) >= FUN2(epsilon)
            y_samples(i) = F2_inv(u(i));
        end
    end
end

% 计算Z值
function Z_t = compute_Z_at_time(FUN1,FUN2,F1_inv,F2_inv,epsilon, N_poisson_rand)
    Z_t = zeros(size(N_poisson_rand));
    for i = 1:length(N_poisson_rand)
        Y_samples = inverse_transform_sampling(FUN1,FUN2,F1_inv,F2_inv,epsilon, N_poisson_rand(i));
        Z_t(i) = sum(Y_samples);
    end
end


function err = errBurgers(x_targ,phi_index ,ksi)% regular
        err = (x_targ - phi_index *ksi) / norm(x_targ);
end


function LLH = likelihood_ksi(x_targ,phi_index ,ksi,sigma_ksi2,mu_ksi,mu_e,sigma_e2)
        err = errBurgers(x_targ,phi_index,ksi);
        Jksi = 1/2*sum(1./sigma_ksi2.*(ksi-mu_ksi).^2);
        Jerr = 1/2*sum((err-mu_e).^2/sigma_e2);
        LLH = -Jksi-Jerr;
        LLH = double(LLH);
end

function ksiRec = bayesian_update_model(Drif, P,Dri_var,FX)%两个系数变量参数
    % 贝叶斯更新模型函数
    % 输入:
    %   Drift2: 漂移系数向量
    %   Phi1: 基函数矩阵
    % 输出:
    %   ksiRec: 采样结果，尺寸为 [参数数量, numSamps * numRuns]

    % 初始化参数
    mu_ksi = Drif(Drif ~= 0); % 非零漂移系数
    Drift_index = find(Drif ~= 0); % 非零漂移系数的索引
    phi_index = P(:, Drift_index); % 对应的基函数
    sigma_ksi2 = Dri_var(Drift_index); % 漂移系数的方差

    sigma_mu_e2 = (1/3)^2; % 误差均值的先验方差
    alpha_e = 1; % 误差方差的先验参数
    beta_e = 2; % 误差方差的先验参数
    Ne = size(phi_index, 1) * size(phi_index, 2); % 样本数量
    sigmaM1 = 0.02; % ksi1的提议分布标准差
    sigmaM2 = 0.005; % ksi2的提议分布标准差
    numSamps = 100; % 每次运行的采样次数
    numRuns = 10; % 运行次数

    % 预分配总的结果存储
    all_ksiRec = zeros(length(mu_ksi), numSamps * numRuns);
    all_mu_eRec = zeros(1, numSamps * numRuns);
    all_sigma_e2Rec = zeros(1, numSamps * numRuns);
    rej_counts1 = zeros(numRuns, 1);
    rej_counts2 = zeros(numRuns, 1);

    tic;
    for run = 1:numRuns

        % 初始化参数（每次运行独立）
        ksi = mu_ksi;
        err = errBurgers(FX, phi_index, ksi);
        mu_e = mean(err);
        sigma_e2 = var(err);

        rej_count1 = 0;
        rej_count2 = 0;

        % 预分配当前运行的结果
        ksiRec = zeros(length(ksi), numSamps);
        mu_eRec = zeros(1, numSamps);
        sigma_e2Rec = zeros(1, numSamps);

        % MH抽样过程
        for n = 1:numSamps
            ksiRec(:, n) = ksi;
            mu_eRec(:, n) = mu_e;
            sigma_e2Rec(:, n) = sigma_e2;

            % 更新ksi1
            ksiNew1 = ksi(1) + randn * sigmaM1;
            while ksiNew1 < 0
                ksiNew1 = ksi(1) + randn * sigmaM1;
            end
            ksiNew = [ksiNew1; ksi(2)];
            LLH1_Old = likelihood_ksi(FX, phi_index, ksi, sigma_ksi2, mu_ksi, mu_e, sigma_e2);
            LLH1_New = likelihood_ksi(FX, phi_index, ksiNew, sigma_ksi2, mu_ksi, mu_e, sigma_e2);
            if rand < exp(LLH1_New - LLH1_Old)
                ksi = ksiNew;
            else
                rej_count1 = rej_count1 + 1;
            end

            % 更新ksi2
            ksiNew2 = ksi(2) + randn * sigmaM2;
            ksiNew = [ksi(1); ksiNew2];
            LLH1_Old = likelihood_ksi(FX, phi_index, ksi, sigma_ksi2, mu_ksi, mu_e, sigma_e2);
            LLH1_New = likelihood_ksi(FX, phi_index, ksiNew, sigma_ksi2, mu_ksi, mu_e, sigma_e2);
            if rand < exp(LLH1_New - LLH1_Old)
                ksi = ksiNew;
            else
                rej_count2 = rej_count2 + 1;
            end

            % 更新mu_e和sigma_e2
            err = errBurgers(FX, phi_index, ksi);
            C1 = Ne / sigma_e2 + 1 / sigma_mu_e2;
            mu_e = sum(err) / sigma_e2 / C1;
            sigma_e2 = (beta_e + 0.5 * sum((err - mu_e).^2)) / (Ne / 2 + alpha_e + 1);
        end

        % 保存当前运行结果
        startIdx = (run - 1) * numSamps + 1;
        endIdx = run * numSamps;
        all_ksiRec(:, startIdx:endIdx) = ksiRec;
        all_mu_eRec(:, startIdx:endIdx) = mu_eRec;
        all_sigma_e2Rec(:, startIdx:endIdx) = sigma_e2Rec;
        rej_counts1(run) = rej_count1;
        rej_counts2(run) = rej_count2;
    end

    % 返回结果
    ksiRec = all_ksiRec; % 尺寸为 [参数数量, 10*numSamps]
end


function ksiRec = bayesian_update_model_1(Drif, P, Dri_var, tx)
    % 贝叶斯更新模型函数（一维ksi）
    % 输入:
    %   Drif: 漂移系数向量
    %   P: 基函数矩阵
    %   Dri_var: 漂移系数的方差
    %   tx: 时间-空间数据
    % 输出:
    %   ksiRec: 采样结果，尺寸为 [1, numSamps * numRuns]

    % 初始化参数
    mu_ksi = Drif(Drif ~= 0); % 非零漂移系数
    Drift_index = find(Drif ~= 0); % 非零漂移系数的索引
    phi_index = P(:, Drift_index); % 对应的基函数
    sigma_ksi2 = Dri_var(Drift_index); % 漂移系数的方差

    sigma_mu_e2 = (1/3)^2; % 误差均值的先验方差
    alpha_e = 1; % 误差方差的先验参数
    beta_e = 2; % 误差方差的先验参数
    Ne = size(phi_index, 1) * size(phi_index, 2); % 样本数量
    sigmaM = 0.02; % ksi的提议分布标准差
    numSamps = 100; % 每次运行的采样次数
    numRuns = 10; % 运行次数

    % 预分配总的结果存储
    all_ksiRec = zeros(1, numSamps * numRuns); % 一维ksi
    all_mu_eRec = zeros(1, numSamps * numRuns);
    all_sigma_e2Rec = zeros(1, numSamps * numRuns);
    rej_counts = zeros(numRuns, 1); % 拒绝计数

    tic;
    for run = 1:numRuns
        fprintf('Run %d/%d\n', run, numRuns);

        % 初始化参数（每次运行独立）
        ksi = mu_ksi; % 一维ksi
        err = errBurgers(tx, phi_index, ksi);
        mu_e = mean(err);
        sigma_e2 = var(err);

        rej_count = 0; % 拒绝计数

        % 预分配当前运行的结果
        ksiRec = zeros(1, numSamps); % 一维ksi
        mu_eRec = zeros(1, numSamps);
        sigma_e2Rec = zeros(1, numSamps);

        % MH抽样过程
        for n = 1:numSamps
            ksiRec(:, n) = ksi;
            mu_eRec(:, n) = mu_e;
            sigma_e2Rec(:, n) = sigma_e2;

            % 更新ksi
            ksiNew = ksi + randn * sigmaM; % 生成新的ksi
            while ksiNew < 0 % 确保ksi非负
                ksiNew = ksi + randn * sigmaM;
            end

            % 计算似然函数
            LLH_Old = likelihood_ksi(tx, phi_index, ksi, sigma_ksi2, mu_ksi, mu_e, sigma_e2);
            LLH_New = likelihood_ksi(tx, phi_index, ksiNew, sigma_ksi2, mu_ksi, mu_e, sigma_e2);

            % Metropolis-Hastings接受/拒绝
            if rand < exp(LLH_New - LLH_Old)
                ksi = ksiNew; % 接受新值
            else
                rej_count = rej_count + 1; % 拒绝新值
            end

            % 更新mu_e和sigma_e2
            err = errBurgers(tx, phi_index, ksi);
            C1 = Ne / sigma_e2 + 1 / sigma_mu_e2;
            mu_e = sum(err) / sigma_e2 / C1;
            sigma_e2 = (beta_e + 0.5 * sum((err - mu_e).^2)) / (Ne / 2 + alpha_e + 1);
        end

        % 保存当前运行结果
        startIdx = (run - 1) * numSamps + 1;
        endIdx = run * numSamps;
        all_ksiRec(:, startIdx:endIdx) = ksiRec;
        all_mu_eRec(:, startIdx:endIdx) = mu_eRec;
        all_sigma_e2Rec(:, startIdx:endIdx) = sigma_e2Rec;
        rej_counts(run) = rej_count;
    end

    elapsed_time = toc;
    fprintf('Total time: %.2f seconds\n', elapsed_time);

    % 返回结果
    ksiRec = all_ksiRec; % 尺寸为 [1, numSamps * numRuns]
end