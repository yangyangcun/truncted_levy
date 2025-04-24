clc
clear
close all


rng(100, 'twister');
alpha = 0.5;
epsilon = 0.05;
num_N = 25000;
dt = 0.04;
T = dt * num_N;
c_alpha = alpha * gamma((1+alpha)/2) / (2^(1-alpha) * pi^(0.5) * gamma(1-alpha/2));
Lt = generate_Lt(alpha,epsilon,num_N,dt,c_alpha);

Xzero = 2;
X_simul = zeros(1, num_N + 1);
X_simul(1) = Xzero;
for i = 2:num_N + 1
        X_simul(i) = X_simul(i-1) + (X_simul(i-1) - X_simul(i-1)^3) * dt + sqrt(4 ) * (Lt(i) - Lt(i-1));
end
%
X = X_simul';
u=X(1:end-1);
order = 3;  % 最大幂次
[deg_funcs, deg_strs] = deg_func(order);
FX =(X(2:end)-X(1:end-1))./dt; 
Phi = zeros(length(u), length(deg_funcs));
for i = 1:length(deg_funcs)
    Phi(:,i) = deg_funcs{i}(u);
end
%
[Drift1,Drift_var1] = MSBL(Phi, FX, 1, 1);
p = 1;
[~,Drift2,Drift_var2] = StepwiseRSBL(Phi, FX, 1, 1, p);

var_alpha = 2 * dt * c_alpha /(2 - alpha);
GX1 = (X(2:end,:)-X(1:end-1,:)-Phi*Drift1*dt).^2 / var_alpha;
GX2 = (X(2:end,:)-X(1:end-1,:)-Phi*Drift2*dt).^2 / var_alpha;

[Diff1,Diff_var1] = MSBL(Phi, GX1, 1, 1);
p2 = 0.01;
[~,Diff2,Diff_var2] = StepwiseRSBL(Phi, GX2, 1, 1,p2);

%%
Drift_true = Drift2;
%%
MSE1 = mse(Drift_true - Drift1);
MSE2 = mse(Drift_true - Drift2);
disp(MSE1);
disp(MSE2);
MAE1 = abs(Drift_true - Drift1);
MAE2 = abs(Drift_true - Drift2);

%%
Diff_true = Diff2;
%%
MSE3 = mse(Diff_true - Diff1);
MSE4 = mse(Diff_true - Diff2);
disp(MSE3);
disp(MSE4);
MAE3 = abs(Diff_true - Diff1);
MAE4 = abs(Diff_true - Diff2);

%% 画样本路径的图
map = colormap(nclCM(232));%color包里选颜色
map = flipud(map);
close all;
figureUnits = 'centimeters';
figureWidth = 14;
figureHeight = 8;
figureHandle = figure;
set(gcf, 'Units', figureUnits, 'Position', [0 0 figureWidth figureHeight]);
hold on
scatter(0:dt:T, X_simul, 10,  X_simul, 'filled')
hXLabel = xlabel('$T$', 'Interpreter', 'latex');
hYLabel = ylabel('$X_t$', 'Interpreter', 'latex');
colormap(map)
colorbar
set(gca, 'Box', 'off', ...                                        % 边框
         'LineWidth',1,...                                        % 线宽
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.005 .005], ...         % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...            % 小刻度
         'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1],...          % 坐标轴颜色
         'XTick', 0:200:1000,...                                    % 坐标区刻度、范围
         'XLim', [-4 1000],...
         'YTick', -4:3:5,...
         'YLim', [-4 5])
% 字体和字号
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15,'FontWeight' , 'normal')
set([hXLabel, hYLabel], 'FontSize', 17, 'FontName', 'SimSun','FontWeight' , 'normal')
set(gcf,'Color',[1 1 1])
xc = get(gca,'XColor');
yc = get(gca,'YColor');
unit = get(gca,'units');
ax = axes( 'Units', unit,...
           'Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor',xc,...
           'YColor',yc);
set(ax, 'linewidth',1,...
        'XTick', [],...
        'YTick', []);
%% 图片输出
figW = figureWidth;%将图形的宽度和高度赋值给变量 figW 和 figH
figH = figureHeight;
set(figureHandle,'PaperUnits',figureUnits);%设置图形的打印单位为 figureUnits
set(figureHandle,'PaperPosition',[0 0 figW figH]);
fileout = '样本路径图';
print(figureHandle,[fileout,'.png'],'-r600','-dpng');
%% 画漂移项的函数库情况（sbl和rsbl）
data = Drift1;
labels = deg_strs;
C1 = addcolorplus(252); 
C2 = addcolorplus(119); 
data1 = Drift2;

close all
figureUnits = 'centimeters';% 图片尺寸设置（单位：厘米）
figureWidth = 28;
figureHeight = 8;
figureHandle = figure;% 窗口设置
set(gcf, 'Units', figureUnits, 'Position', [0 0 figureWidth figureHeight]);

st = stem(data,...
    'MarkerEdgeColor', C1, ... % 符号轮廓颜色
    'MarkerFaceColor', C1, ... % 符号填充颜色
    'Marker', 'd', ... % 符号类型
    'MarkerSize', 4, ... % 符号尺寸
    'LineWidth', 2, ... % 线宽
    'LineStyle', '-', ... % 线型
    'Color', C1); % 线的颜色
hold on; % 保持当前图形
for i = 1:length(data1)
    if abs(data1(i)) >= 0.01 % 只绘制绝对值大于等于 0.001 的数据点
        stem(i, data1(i), ...
            'MarkerEdgeColor', C2, ... % 符号轮廓颜色
            'MarkerFaceColor', C2, ... % 符号填充颜色
            'Marker', 'd', ... % 符号类型
            'MarkerSize', 4, ... % 符号尺寸
            'LineWidth', 2, ... % 线宽
            'LineStyle', '--', ... % 线型
            'Color', C2); % 线的颜色
    end
end

for i = 1:length(data)
    [maxVal, maxIdx] = max([abs(data(i)), abs(data1(i))]);
    if maxIdx == 1
        labelPos = data(i);
    else
        labelPos = data1(i);
    end
    if maxVal >= 0.01 
        if labelPos > 0
            text(i, labelPos+ 0.1, labels{i}, ... % y 方向增加 0.1 的偏移量
                'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 18, ...
                'Color', 'k', ...
                'Interpreter', 'tex', 'FontName', 'Times New Roman'); % 使用 LaTeX 语法支持上标
        
        else
            text(i, labelPos-0.1, labels{i}, ... % y 方向减少 0.15 的偏移量
                'VerticalAlignment', 'top', ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 18, ...
                'Color', 'k', ...
                'Interpreter', 'tex', 'FontName', 'Times New Roman'); % 使用 LaTeX 语法支持上标
        end
    end
end
legend({'SBL', 'AT-RSBL'}, ...
    'Location', 'northeast', ... % 图例位置在右上角
    'FontSize', 15, ...
    'Box', 'off','FontWeight' , 'normal','FontName', 'Times New Roman'); % 去掉图例边框
hXLabel = xlabel('Library functions');
hYLabel = ylabel('Component coefficients');

set(gca,'looseInset',[0 0 0 0],'Box', 'off', ... % 边框
    'LineWidth', 1, ... % 线宽
    'XGrid', 'off', 'YGrid', 'off', ... % 网格
    'TickDir', 'out', 'TickLength', [.01 .01], ... % 刻度
    'XMinorTick', 'off', 'YMinorTick', 'off', ... % 小刻度
    'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1],'YTick', -1.7:1.7:1.7, ...
         'YLim', [-1.7 1.7], 'FontName', 'Times New Roman', 'FontSize', 15,'FontWeight' , 'normal') % 坐标轴颜色
xc = get(gca,'XColor');
yc = get(gca,'YColor');
unit = get(gca,'units');
ax = axes( 'Units', unit,...
           'Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor',xc,...
           'YColor',yc);
set(ax, 'linewidth',1,...
        'XTick', [],...
        'YTick', []);
set([hXLabel, hYLabel], 'FontSize', 15, 'FontName', 'Times New Roman','FontWeight' , 'normal')
set(gcf, 'Color', [1 1 1]);% 背景颜色

%% 图片输出
figW = figureWidth;
figH = figureHeight;
set(figureHandle,'PaperUnits',figureUnits);
set(figureHandle,'PaperPosition',[0 0 figW figH]);
fileout = '漂移项的函数项情况';
print(figureHandle,[fileout,'.png'],'-r600','-dpng');

%% 漂移、耗散项的分布图
C1 = addcolorplus(252); 
C2 = addcolorplus(119); 
close all
mu1 = Drift2(2); 
sigma1 = sqrt(Drift_var2(2));
mu2 = Drift2(4); 
sigma2 = sqrt(Drift_var2(4));
mu3 = Diff2(1); 
sigma3 = sqrt(Diff_var2(1)); 

x1 = linspace(0, 2, 2000); 
x2 = linspace(-2, 0, 2000); 
x3 = linspace(3, 5, 2000); 
pdf1 = normpdf(x1, mu1, sigma1); 
pdf2 = normpdf(x2, mu2, sigma2); 
pdf3 = normpdf(x3, mu3, sigma3); 

figureUnits = 'centimeters';
figureWidth = 28;
figureHeight = 10;
figureHandle = figure;
set(gcf, 'Units', figureUnits, 'Position', [2 5 figureWidth figureHeight]);

% 调整 subplot 的位置和间距
subplotGap = 0.06; 
subplotWidth = (1 - 3.5 * subplotGap) / 3; 
subplotHeight = 0.7;  % 降低子图高度
subplotBottom = 0.25;  % 提升子图起始位置

% ========================= 第一个子图 =========================
subplot('Position', [subplotGap, subplotBottom, subplotWidth, subplotHeight]);
hold on
area(x1, pdf1, 'LineWidth', 2, 'FaceColor', C2, 'EdgeColor', C2, ...
    'FaceAlpha', .3, 'EdgeAlpha', 1);
line([mu1 mu1], [0 normpdf(mu1, mu1, sigma1)], ...
    'LineStyle', '--', 'Color', C2, 'LineWidth', 1.5);
text(mu1, normpdf(mu1, mu1, sigma1) + 3, sprintf('N(%.4f, %.4f)', mu1, sigma1^2), ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 15, ...
    'Color', 'black', 'FontName', 'Times New Roman');

set(gca, 'Box', 'off', ...
         'XGrid', 'off', 'YGrid', 'off', ...
         'TickDir', 'out', 'TickLength', [.01 .01], ...
         'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1],...
         'XTick', 0.7:0.3:1.3,...
         'XLim', [0.7 1.3],...
         'YLim', [0 35]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'normal');
hXLabel1 = xlabel('$\theta_1$', 'Interpreter', 'latex');
hYLabel1 = ylabel('PDF');
set(hXLabel1, 'FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal')
set(hYLabel1, 'FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal')

% 第二个 subplot
subplot('Position', [2*subplotGap+subplotWidth, subplotBottom, subplotWidth, subplotHeight]);
hold on
area(x2, pdf2, 'LineWidth', 2, 'FaceColor', C2, 'EdgeColor', C2, ...
    'FaceAlpha', .3, 'EdgeAlpha', 1);
line([mu2 mu2], [0 normpdf(mu2, mu2, sigma2)], ...
    'LineStyle', '--', 'Color', C2, 'LineWidth', 1.5);
text(mu2, normpdf(mu2, mu2, sigma2) + 5, sprintf('N(%.4f, %.5f)', mu2, sigma2^2), ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 15, ...
    'Color', 'black', 'FontName', 'Times New Roman');

set(gca, 'Box', 'off', ...
         'XGrid', 'off', 'YGrid', 'off', ...
         'TickDir', 'out', 'TickLength', [.01 .01], ...
         'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1],...
         'XTick', -1.2:0.2:-0.8,...
         'XLim', [-1.2 -0.8],...
         'YLim', [0 80]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'normal');
hXLabel2 = xlabel('$\theta_2$', 'Interpreter', 'latex');
set(hXLabel2, 'FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal')

% 第三个 subplot
subplot('Position', [3*subplotGap+2*subplotWidth, subplotBottom, subplotWidth, subplotHeight]);
hold on
area(x3, pdf3, 'LineWidth', 2, 'FaceColor', C2, 'EdgeColor', C2, ...
    'FaceAlpha', .3, 'EdgeAlpha', 1);
line([mu3 mu3], [0 normpdf(mu3, mu3, sigma3)], ...
    'LineStyle', '--', 'Color', C2, 'LineWidth', 1.5);
text(mu3, normpdf(mu3, mu3, sigma3) + 0.4, sprintf('N(%.4f, %.4f)', mu3, sigma3^2), ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 15, ...
    'Color', 'black', 'FontName', 'Times New Roman');

set(gca, 'Box', 'off', ...
         'XGrid', 'off', 'YGrid', 'off', ...
         'TickDir', 'out', 'TickLength', [.01 .01], ...
         'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1],...
         'XTick', 3:0.5:5,...
         'XLim', [3 5],...
         'YLim', [0 3.5]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'normal');
hXLabel3 = xlabel('$\theta_3$', 'Interpreter', 'latex');
set(hXLabel3, 'FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal')

% 设置背景颜色为白色
set(gcf, 'Color', [1 1 1]);
set(gcf, 'PaperPositionMode', 'auto'); 
set(gcf, 'InvertHardcopy', 'off');
%% 图片输出
figW = figureWidth;
figH = figureHeight;
set(figureHandle,'PaperUnits',figureUnits);
set(figureHandle,'PaperPosition',[0 0 figW figH]);
fileout = '漂移和耗散项的分布';
print(figureHandle,[fileout,'.png'],'-r500','-dpng');

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

function ksiRec = bayesian_update_model(Drif, P,Dri_var,tx)%两个系数变量参数
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
        err = errBurgers(tx(2,:)', phi_index, ksi);
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
            LLH1_Old = likelihood_ksi(tx(2,:)', phi_index, ksi, sigma_ksi2, mu_ksi, mu_e, sigma_e2);
            LLH1_New = likelihood_ksi(tx(2,:)', phi_index, ksiNew, sigma_ksi2, mu_ksi, mu_e, sigma_e2);
            if rand < exp(LLH1_New - LLH1_Old)
                ksi = ksiNew;
            else
                rej_count1 = rej_count1 + 1;
            end

            % 更新ksi2
            ksiNew2 = ksi(2) + randn * sigmaM2;
            ksiNew = [ksi(1); ksiNew2];
            LLH1_Old = likelihood_ksi(tx(2,:)', phi_index, ksi, sigma_ksi2, mu_ksi, mu_e, sigma_e2);
            LLH1_New = likelihood_ksi(tx(2,:)', phi_index, ksiNew, sigma_ksi2, mu_ksi, mu_e, sigma_e2);
            if rand < exp(LLH1_New - LLH1_Old)
                ksi = ksiNew;
            else
                rej_count2 = rej_count2 + 1;
            end

            % 更新mu_e和sigma_e2
            err = errBurgers(tx(2,:)', phi_index, ksi);
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
        err = errBurgers(tx(2,:)', phi_index, ksi);
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
            LLH_Old = likelihood_ksi(tx(2,:)', phi_index, ksi, sigma_ksi2, mu_ksi, mu_e, sigma_e2);
            LLH_New = likelihood_ksi(tx(2,:)', phi_index, ksiNew, sigma_ksi2, mu_ksi, mu_e, sigma_e2);

            % Metropolis-Hastings接受/拒绝
            if rand < exp(LLH_New - LLH_Old)
                ksi = ksiNew; % 接受新值
            else
                rej_count = rej_count + 1; % 拒绝新值
            end

            % 更新mu_e和sigma_e2
            err = errBurgers(tx(2,:)', phi_index, ksi);
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