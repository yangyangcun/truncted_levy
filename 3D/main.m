clc
clear
close all

rng(100, 'twister');
alpha = 1;
epsilon = 0.05;
num_N = 10000;
dt = 0.01;
T = dt * num_N;
c_alpha = alpha * gamma((1+alpha)/2) / (2^(1-alpha) * pi^(0.5) * gamma(1-alpha/2));
Lt1 = generate_Lt(alpha,epsilon,num_N,dt,c_alpha);
Lt2 = generate_Lt(alpha,epsilon,num_N,dt,c_alpha);
Lt3 = generate_Lt(alpha,epsilon,num_N,dt,c_alpha);

Xzero = 0;
Yzero = 0;
Zzero = 0;
X_simul = zeros(num_N + 1,1);
Y_simul = zeros(num_N + 1,1);
Z_simul = zeros(num_N + 1,1);
X_simul(1) = Xzero;
Y_simul(1) = Yzero;
Z_simul(1) = Zzero;
for i = 2:num_N + 1
        X_simul(i) = X_simul(i-1) + ( -X_simul(i-1) + Y_simul(i-1)) * dt + sqrt(1) * (Lt1(i) - Lt1(i-1));
        Y_simul(i) = Y_simul(i-1) + ( 0.5*X_simul(i-1) - Y_simul(i-1) -X_simul(i-1)*Z_simul(i-1)) * dt + sqrt(1) * (Lt2(i) - Lt2(i-1));
        Z_simul(i) = Z_simul(i-1) + ( -(8/3)*Z_simul(i-1) + X_simul(i-1)*Y_simul(i-1)) * dt + sqrt(1) * (Lt3(i) - Lt3(i-1));
end

any(isinf(X_simul))   % 检查是否有Inf

figure;
plot(0:dt:dt * num_N,X_simul, 'Color', 'b');
hold on;
plot(0:dt:dt * num_N,Y_simul, 'Color', 'r');
hold on;
plot(0:dt:dt * num_N,Z_simul, 'Color', 'g');
grid on;

X = X_simul(1:end-1);
Y = Y_simul(1:end-1);
Z = Z_simul(1:end-1);
FX_X =(X_simul(2:end)-X_simul(1:end-1))./dt; 
FX_Y =(Y_simul(2:end)-Y_simul(1:end-1))./dt; 
FX_Z =(Z_simul(2:end)-Z_simul(1:end-1))./dt; 


order = 2;  % 最大幂次
[deg_funcs, deg_strs] = deg_func_xy(order);
num_funcs = length(deg_funcs);
Phi = zeros(length(X), num_funcs); 
for k = 1:num_funcs
    func = deg_funcs{k}; % 获取第 k 个函数句柄
    Phi(:, k) = func(X , Y ,Z); % 计算函数值并存储
end
Phi_old = Phi;

[~,Drift1_X,Drift_var1_X] = MSBL(Phi, FX_X, 1, 1);
[~,Drift1_Y,Drift_var1_Y] = MSBL(Phi, FX_Y, 1, 1);
[~,Drift1_Z,Drift_var1_Z] = MSBL(Phi, FX_Z, 1, 1);
p = 1;
[~,Drift2_X,Drift_var2_X] = StepwiseRSBL(Phi, FX_X, 1, 1, p);
[~,Drift2_Y,Drift_var2_Y] = StepwiseRSBL(Phi, FX_Y, 1, 1, p);
[~,Drift2_Z,Drift_var2_Z] = StepwiseRSBL(Phi, FX_Z, 1, 1, p);

var_alpha = 2 * dt * c_alpha /(2 - alpha);
XG =  X_simul(1:end-1);
YG = Y_simul(1:end-1);
ZG = Z_simul(1:end-1);
GX1_X = ((X_simul(2:end,:)- X_simul(1:end-1,:)-Phi*Drift1_X*dt)).^2/ var_alpha;
GY1_Y = ((Y_simul(2:end,:) - Y_simul(1:end-1,:) - Phi*Drift1_Y*dt)).^2/ var_alpha;
GZ1_Z = ((Z_simul(2:end,:) - Z_simul(1:end-1,:) - Phi*Drift1_Z*dt)).^2/ var_alpha;

GX2_X = ((X_simul(2:end,:)- X_simul(1:end-1,:)-Phi*Drift2_X*dt)).^2/ var_alpha;
GY2_Y = ((Y_simul(2:end,:) - Y_simul(1:end-1,:) - Phi*Drift2_Y*dt)).^2/ var_alpha;
GZ2_Z = ((Z_simul(2:end,:) - Z_simul(1:end-1,:) - Phi*Drift2_Z*dt)).^2/ var_alpha;

[deg_funcs1, deg_strs1] = deg_func_xy(1);
num_funcs1 = length(deg_funcs1);
Phi1 = zeros(length(X), num_funcs1); 
for k = 1:num_funcs1
    func = deg_funcs1{k}; % 获取第 k 个函数句柄
    Phi1(:, k) = func(X , Y ,Z); % 计算函数值并存储
end
[~,Diff1_X,Diff_var1_X] = MSBL(Phi1, GX1_X, 1, 1);
[~,Diff1_Y,Diff_var1_Y] = MSBL(Phi1, GY1_Y, 1, 1);
[~,Diff1_Z,Diff_var1_Z] = MSBL(Phi1, GZ1_Z, 1, 1);
p2 = 0.02;
[~,Diff2_X,Diff_var2_X] = StepwiseRSBL(Phi1, GX2_X, 1, 1, p2);
[~,Diff2_Y,Diff_var2_Y] = StepwiseRSBL(Phi1, GY2_Y, 1, 1, p2);
[~,Diff2_Z,Diff_var2_Z] = StepwiseRSBL(Phi1, GZ2_Z, 1, 1, p2);
%%
Drift1_true = Drift2_X;
Drift2_true = Drift2_Y;
Drift3_true = Drift2_Z;
%%
MSE1 = mse(Drift1_true - Drift1_X);
MSE2 = mse(Drift1_true - Drift2_X);
MSE3 = mse(Drift2_true - Drift1_Y);
MSE4 = mse(Drift2_true - Drift2_Y);
MSE5 = mse(Drift3_true - Drift1_Z);
MSE6 = mse(Drift3_true - Drift2_Z);
disp(MSE1);
disp(MSE2);
disp(MSE3);
disp(MSE4);
disp(MSE5);
disp(MSE6);
MAE1 = abs(Drift1_true - Drift1_X);
MAE2 = abs(Drift1_true - Drift2_X);
MAE3 = abs(Drift2_true - Drift1_Y);
MAE4 = abs(Drift2_true - Drift2_Y);
MAE5 = abs(Drift3_true - Drift1_Z);
MAE6 = abs(Drift3_true - Drift2_Z);
%%
Diff1_true = Diff2_X;
Diff2_true = Diff2_Y;
Diff3_true = Diff2_Z;
%%
MSE7 = mse(Diff1_true - Diff1_X);
MSE8 = mse(Diff1_true - Diff2_X);
MSE9 = mse(Diff2_true - Diff1_Y);
MSE10 = mse(Diff2_true - Diff2_Y);
MSE11 = mse(Diff3_true - Diff1_Z);
MSE12 = mse(Diff3_true - Diff2_Z);

disp(MSE7);
disp(MSE8);
disp(MSE9);
disp(MSE10);
disp(MSE11);
disp(MSE12);

MAE7 = abs(Diff1_true - Diff1_X);
MAE8 = abs(Diff1_true - Diff2_X);
MAE9 = abs(Diff2_true - Diff1_Y);
MAE10 = abs(Diff2_true - Diff2_Y);
MAE11 = abs(Diff3_true - Diff1_Z);
MAE12 = abs(Diff3_true - Diff2_Z);
%% 画样本路径图
map = colormap(nclCM(232)); % color包里选颜色
map = flipud(map);
close all;

% 设置图形尺寸和边距
figureUnits = 'centimeters';
figureWidth = 28; 
figureHeight = 10;
leftMargin = 0.08;    % 左边距
rightMargin = 0.04;   % 右边距
bottomMargin = 0.25;  % 增加底边距给标签留空间
topMargin = 0.01;      % 顶边距
hSpace = 0.07;         % 子图间水平间距

% 计算子图宽度
subplotWidth = (1 - leftMargin - rightMargin - 2*hSpace)/3;

% 创建图形窗口
figureHandle = figure;
set(gcf, 'Units', figureUnits, 'Position', [2 5 figureWidth figureHeight], 'Color', [1 1 1]);
T = 0:dt:dt * num_N;

% ========== 第一个子图 ==========
subplotPos = [leftMargin, bottomMargin, subplotWidth, 0.9-bottomMargin];
subplot('Position', subplotPos);
hold on;
scatter(T, X_simul, 5, X_simul, 'filled');
hXLabel1 = xlabel('$T$', 'Interpreter', 'latex');
hYLabel1 = ylabel('$X_t$', 'Interpreter', 'latex');
colormap(map);
colorbar;
set(gca, 'Box', 'off', 'LineWidth', 1, ...
         'XGrid', 'off', 'YGrid', 'off', ...
         'TickDir', 'out', 'TickLength', [.005 .005], ...
         'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', 0:50:100, 'XLim', [-1 101], ...
         'YTick', -2.2:2:2.5, 'YLim', [-2.2 2.5]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'normal')
set([hXLabel1, hYLabel1], 'FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal')


% ========== 第二个子图 ==========
subplotPos = [leftMargin + subplotWidth + hSpace, bottomMargin, subplotWidth, 0.9-bottomMargin];
subplot('Position', subplotPos);
hold on;
scatter(T, Y_simul, 5, Y_simul, 'filled');
hXLabel2 = xlabel('$T$', 'Interpreter', 'latex');
hYLabel2 = ylabel('$Y_t$', 'Interpreter', 'latex');
colormap(map);
colorbar;
set(gca, 'Box', 'off', 'LineWidth', 1, ...
         'XGrid', 'off', 'YGrid', 'off', ...
         'TickDir', 'out', 'TickLength', [.005 .005], ...
         'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', 0:50:100, 'XLim', [-1 101], ...
         'YTick', -2:2:2, 'YLim', [-2 2]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'normal')
set([hXLabel2, hYLabel2], 'FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal')

% ========== 第三个子图 ==========
subplotPos = [leftMargin + 2*(subplotWidth + hSpace), bottomMargin, subplotWidth, 0.9-bottomMargin];
subplot('Position', subplotPos);
hold on;
scatter(T, Z_simul, 5, Z_simul, 'filled');
hXLabel3 = xlabel('$T$', 'Interpreter', 'latex');
hYLabel3 = ylabel('$Z_t$', 'Interpreter', 'latex');
colormap(map);
colorbar;
set(gca, 'Box', 'off', 'LineWidth', 1, ...
         'XGrid', 'off', 'YGrid', 'off', ...
         'TickDir', 'out', 'TickLength', [.005 .005], ...
         'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', 0:50:100, 'XLim', [-1 101], ...
         'YTick', -2:2:2, 'YLim', [-2 2]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'normal')
set([hXLabel3, hYLabel3], 'FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal')

%% 图片输出
figW = figureWidth;%将图形的宽度和高度赋值给变量 figW 和 figH
figH = figureHeight;
set(figureHandle,'PaperUnits',figureUnits);%设置图形的打印单位为 figureUnits
set(figureHandle,'PaperPosition',[0 0 figW figH]);
fileout = '样本路径图';
print(figureHandle,[fileout,'.png'],'-r600','-dpng');

%% 画漂移项的函数库情况（sbl和rsbl）

labels = deg_strs;
C1 = addcolorplus(252); 
C2 = addcolorplus(119); 
data = Drift1_X;
data1 = Drift2_X;
data2 = Drift1_Y;
data3 = Drift2_Y;
data4 = Drift1_Z;
data5 = Drift2_Z;

close all
figureUnits = 'centimeters';% 图片尺寸设置（单位：厘米）
figureWidth = 28;
figureHeight = 24;
figureHandle = figure;% 窗口设置
set(gcf, 'Units', figureUnits, 'Position', [0 0 figureWidth figureHeight]);

subplotGap = 0.11; % 子图间垂直间距
bottomMargin = 0.13; % 底部边距
topMargin = 0.05; % 顶部边距
subplotHeight = (1 - topMargin - bottomMargin - 2*subplotGap)/3 ;

% ========== 第一个子图 ==========
subplot('Position', [0.1, bottomMargin+2*(subplotHeight+subplotGap), 0.8, subplotHeight]);
hold on
st = stem(data,...
    'MarkerEdgeColor', C1,'MarkerFaceColor', C1,'Marker', 'd','MarkerSize', 4, ... % 符号尺寸
    'LineWidth', 2,'LineStyle', '-','Color', C1); % 线的颜色
hold on; 
for i = 1:length(data1)
    if abs(data1(i)) >= 0.001 % 只绘制绝对值大于等于 0.001 的数据点
        stem(i, data1(i), ...
            'MarkerEdgeColor', C2,'MarkerFaceColor', C2,'Marker', 'd','MarkerSize', 4,'LineWidth', 2, ... % 线宽
            'LineStyle', '--','Color', C2); % 线的颜色
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
            text(i, labelPos + 0.3, labels{i}, ... % y 方向增加 0.1 的偏移量
                'VerticalAlignment', 'bottom','HorizontalAlignment', 'center','FontSize', 18, ...
                'Color', 'k','Interpreter', 'tex', 'FontName', 'Times New Roman'); % 使用 LaTeX 语法支持上标
        else
            text(i, labelPos - 0.1, labels{i}, ... % y 方向减少 0.15 的偏移量
                'VerticalAlignment', 'top','HorizontalAlignment', 'center','FontSize', 18, ...
                'Color', 'k','Interpreter', 'tex', 'FontName', 'Times New Roman'); % 使用 LaTeX 语法支持上标
        end
    end
end

legend({'SBL', 'AT-RSBL'}, ...
    'Location', 'northwest', ... % 图例位置在右上角
    'FontSize', 15, ...
    'Box', 'off','FontWeight' , 'normal', 'FontName', 'Times New Roman'); % 去掉图例边框

hXLabel = xlabel('Library functions');
hYLabel = ylabel('Coefficients (1st component)');
set(gca, 'Box', 'off', ... % 边框
    'LineWidth', 1, ... % 线宽
    'XGrid', 'off', 'YGrid', 'off', ... % 网格
    'TickDir', 'out', 'TickLength', [.01 .01], ... % 刻度
    'XMinorTick', 'off', 'YMinorTick', 'off', ... % 小刻度
    'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], 'FontName', 'Times New Roman', 'FontSize', 15,'FontWeight' , 'bold', ...
    'YTick', -2:2:2,'YLim', [-2 2]) % 坐标轴颜色
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
set( hYLabel, 'FontSize', 15, 'FontName', 'Times New Roman','FontWeight' , 'normal')
set( hXLabel, 'FontSize', 15, 'FontName', 'Times New Roman','FontWeight' , 'normal')

subplot('Position', [0.1, bottomMargin+(subplotHeight+subplotGap), 0.8, subplotHeight]);
hold on
st = stem(data2,...
    'MarkerEdgeColor', C1,'MarkerFaceColor', C1,'Marker', 'd','MarkerSize', 4, ... % 符号尺寸
    'LineWidth', 2,'LineStyle', '-','Color', C1); % 线的颜色
hold on; 
for i = 1:length(data3)
    if abs(data3(i)) >= 0.001 % 只绘制绝对值大于等于 0.001 的数据点
        stem(i, data3(i), ...
            'MarkerEdgeColor', C2,'MarkerFaceColor', C2,'Marker', 'd','MarkerSize', 4,'LineWidth', 2, ... % 线宽
            'LineStyle', '--','Color', C2); % 线的颜色
    end
end

for i = 1:length(data2)
    [maxVal, maxIdx] = max([abs(data2(i)), abs(data3(i))]);
    if maxIdx == 1
        labelPos = data2(i);
    else
        labelPos = data3(i);
    end
    if maxVal >= 0.01 
        if i == 14 || i == 16 % 特殊处理 i = 14 或 i = 16
            text(i, -0.2, labels{i}, ... % 将标签放在横轴下方，y 坐标为 -0.2
                'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', 18, ...
                'Color', 'k', 'Interpreter', 'tex', 'FontName', 'Times New Roman'); 
        else
            if labelPos > 0
                text(i, labelPos + 0.1, labels{i}, ... % y 方向增加 0.1 的偏移量
                    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 18, ...
                    'Color', 'k', 'Interpreter', 'tex', 'FontName', 'Times New Roman'); 
            else
                text(i, labelPos - 0.1, labels{i}, ... 
                    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', 18, ...
                    'Color', 'k', 'Interpreter', 'tex', 'FontName', 'Times New Roman'); 
            end
        end
    end
end

hXLabel = xlabel('Library functions');
hYLabel = ylabel('Coefficients (2nd component)');
set(gca, 'Box', 'off', ... % 边框
    'LineWidth', 1, ... % 线宽
    'XGrid', 'off', 'YGrid', 'off', ... % 网格
    'TickDir', 'out', 'TickLength', [.01 .01], ... % 刻度
    'XMinorTick', 'off', 'YMinorTick', 'off', ... % 小刻度
    'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], 'FontName', 'Times New Roman', 'FontSize', 15,'FontWeight' , 'bold', ...
    'YTick', -2:2:2,'YLim', [-2 2]) % 坐标轴颜色
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


subplot('Position', [0.1, bottomMargin, 0.8, subplotHeight]);
hold on
st = stem(data4,...
    'MarkerEdgeColor', C1,'MarkerFaceColor', C1,'Marker', 'd','MarkerSize', 4, ... % 符号尺寸
    'LineWidth', 2,'LineStyle', '-','Color', C1); % 线的颜色
hold on; 
for i = 1:length(data5)
    if abs(data5(i)) >= 0.001 % 只绘制绝对值大于等于 0.001 的数据点
        stem(i, data5(i), ...
            'MarkerEdgeColor', C2,'MarkerFaceColor', C2,'Marker', 'd','MarkerSize', 4,'LineWidth', 2, ... % 线宽
            'LineStyle', '--','Color', C2); % 线的颜色
    end
end

for i = 1:length(data4)
    [maxVal, maxIdx] = max([abs(data4(i)), abs(data5(i))]);
    if maxIdx == 1
        labelPos = data4(i);
    else
        labelPos = data5(i);
    end
    if maxVal >= 0.01 
        if labelPos > 0
            text(i, labelPos + 0.3, labels{i}, ... % y 方向增加 0.1 的偏移量
                'VerticalAlignment', 'bottom','HorizontalAlignment', 'center','FontSize', 18, ...
                'Color', 'k','Interpreter', 'tex', 'FontName', 'Times New Roman'); % 使用 LaTeX 语法支持上标
        else
            text(i, labelPos - 0.1, labels{i}, ... % y 方向减少 0.15 的偏移量
                'VerticalAlignment', 'top','HorizontalAlignment', 'center','FontSize', 18, ...
                'Color', 'k','Interpreter', 'tex', 'FontName', 'Times New Roman'); % 使用 LaTeX 语法支持上标
        end
    end
end

hXLabel = xlabel('Library functions');
hYLabel = ylabel('Coefficients (3rd component)');
set(gca, 'Box', 'off', ... % 边框
    'LineWidth', 1, ... % 线宽
    'XGrid', 'off', 'YGrid', 'off', ... % 网格
    'TickDir', 'out', 'TickLength', [.01 .01], ... % 刻度
    'XMinorTick', 'off', 'YMinorTick', 'off', ... % 小刻度
    'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], 'FontName', 'Times New Roman', 'FontSize', 15,'FontWeight' , 'bold', ...
    'YTick', -4:4:4,'YLim', [-4 4]) % 坐标轴颜色
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

%% 画耗散项的分布图

C2 = addcolorplus(119); 
close all;
mu1 = Diff2_X(1); 
sigma1 = sqrt(Diff_var2_X(1));
mu2 = Diff2_Y(1); 
sigma2 = sqrt(Diff_var2_Y(1));
mu3 = Diff2_Z(1); 
sigma3 = sqrt(Diff_var2_Z(1));

x = linspace(0, 2, 2000); 
pdf1 = normpdf(x, mu1, sigma1); 
pdf2 = normpdf(x, mu2, sigma2); 
pdf3 = normpdf(x, mu3, sigma3); 

% 图形设置
figureUnits = 'centimeters';
figureWidth = 28;
figureHeight = 8;
figureHandle = figure;
set(gcf, 'Units', figureUnits, 'Position', [2 5 figureWidth figureHeight], 'Color', [1 1 1]);

% 子图布局参数
subplotGap = 0.06; % 子图间水平间距
subplotWidth = (1 - 3.5 * subplotGap) / 3; % 子图宽度
subplotHeight = 0.65; % 子图高度
labelOffset = 0.1; % 标签垂直偏移量

% ========== 第一个子图 ==========
subplotPos = [subplotGap, 0.3, subplotWidth, subplotHeight];
subplot('Position', subplotPos); 
hold on;
area(x, pdf1, 'LineWidth', 2, 'FaceColor', C2, 'EdgeColor', C2, ...
    'FaceAlpha', .3, 'EdgeAlpha', 1);
line([mu1 mu1], [0 normpdf(mu1, mu1, sigma1)], ...
    'LineStyle', '--', 'Color', C2, 'LineWidth', 1.5);
text(mu1, normpdf(mu1, mu1, sigma1) + 1, sprintf('N(%.4f, %.4f)', mu1, sigma1^2), ...
    'HorizontalAlignment', 'center', 'FontSize', 15, ...
    'Color', 'black', 'FontName', 'Times New Roman');

set(gca, 'Box', 'off', 'XGrid', 'off', 'YGrid', 'off', ...
    'TickDir', 'out', 'TickLength', [.01 .01], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', ...
    'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
    'XTick', 0.5:0.5:1.5, 'XLim', [0.5 1.5], ...
    'YTick', 0:4:10, 'YLim', [0 10]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold');

hXLabel1 = xlabel('$\theta_3$', 'Interpreter', 'latex');
hYLabel1 = ylabel('PDF');
set([hXLabel1, hYLabel1], 'FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');

% ========== 第二个子图 ==========
subplotPos = [2*subplotGap + subplotWidth, 0.3, subplotWidth, subplotHeight];
subplot('Position', subplotPos); 
hold on;
area(x, pdf2, 'LineWidth', 2, 'FaceColor', C2, 'EdgeColor', C2, ...
    'FaceAlpha', .3, 'EdgeAlpha', 1);
line([mu2 mu2], [0 normpdf(mu2, mu2, sigma2)], ...
    'LineStyle', '--', 'Color', C2, 'LineWidth', 1.5);
text(mu2, normpdf(mu2, mu2, sigma2) + 1, sprintf('N(%.4f, %.5f)', mu2, sigma2^2), ...
    'HorizontalAlignment', 'center', 'FontSize', 15, ...
    'Color', 'black', 'FontName', 'Times New Roman');

set(gca, 'Box', 'off', 'XGrid', 'off', 'YGrid', 'off', ...
    'TickDir', 'out', 'TickLength', [.01 .01], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', ...
    'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
    'XTick', 0.5:0.5:1.5, 'XLim', [0.5 1.5], ...
    'YTick', 0:4:10, 'YLim', [0 10]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'normal');

hXLabel2 = xlabel('$\theta_7$', 'Interpreter', 'latex');
set(hXLabel2, 'FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');

% ========== 第三个子图 ==========
subplotPos = [3*subplotGap + 2*subplotWidth, 0.3, subplotWidth, subplotHeight];
subplot('Position', subplotPos); 
hold on;
area(x, pdf3, 'LineWidth', 2, 'FaceColor', C2, 'EdgeColor', C2, ...
    'FaceAlpha', .3, 'EdgeAlpha', 1);
line([mu3 mu3], [0 normpdf(mu3, mu3, sigma3)], ...
    'LineStyle', '--', 'Color', C2, 'LineWidth', 1.5);
text(mu3, normpdf(mu3, mu3, sigma3) + 1, sprintf('N(%.4f, %.4f)', mu3, sigma3^2), ...
    'HorizontalAlignment', 'center', 'FontSize', 15, ...
    'Color', 'black', 'FontName', 'Times New Roman');

set(gca, 'Box', 'off', 'XGrid', 'off', 'YGrid', 'off', ...
    'TickDir', 'out', 'TickLength', [.01 .01], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', ...
    'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
    'XTick', 0.5:0.5:1.5, 'XLim', [0.5 1.5], ...
    'YTick', 0:4:8, 'YLim', [0 8]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'normal');

hXLabel3 = xlabel('$\theta_{10}$', 'Interpreter', 'latex');
set(hXLabel3, 'FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');

% 设置打印选项
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'InvertHardcopy', 'off');
%% 图片输出
figW = figureWidth;
figH = figureHeight;
set(figureHandle,'PaperUnits',figureUnits);
set(figureHandle,'PaperPosition',[0 0 figW figH]);
fileout = '耗散项的分布';
print(figureHandle,[fileout,'.png'],'-r600','-dpng');


%% 画漂移项系数之间的联合后验分布（RSBL）
map = colormap(nclCM(232)); % color包里选颜色
map = flipud(map);

A = struct();
A(1).mu2 = Drift2_X(2); 
A(1).sigma2 = sqrt(Drift_var2_X(2));
A(2).mu2 = Drift2_X(4); 
A(2).sigma2 = sqrt(Drift_var2_X(4));
A(3).mu2 = Drift2_Y(2); 
A(3).sigma2 = sqrt(Drift_var2_Y(2));
A(4).mu2 = Drift2_Y(4); 
A(4).sigma2 = sqrt(Drift_var2_Y(4));
A(5).mu2 = Drift2_Y(12); 
A(5).sigma2 = sqrt(Drift_var2_Y(12));
A(6).mu2 = Drift2_Z(6); 
A(6).sigma2 = sqrt(Drift_var2_Z(6));
A(7).mu2 = Drift2_Z(8); 
A(7).sigma2 = sqrt(Drift_var2_Z(8));

formatTickLabels = @(ticks) arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
close all;
figureUnits = 'centimeters';
figureWidth = 28; 
figureHeight = 16;
figureHandle = figure;
set(gcf, 'Units', figureUnits, 'Position', [0 0 figureWidth figureHeight], 'Color', [1 1 1]);

% 定义子图的位置和标签
positions = [
    0.08, 0.62, 0.24, 0.31; % 第一张图
    0.40, 0.62, 0.24, 0.31; % 第二张图
    0.72, 0.62, 0.24, 0.31; % 第三张图
    0.22, 0.14, 0.24, 0.31; % 第四张图
    0.55, 0.14, 0.24, 0.31; % 第五张图
];

xLabels = {'$\theta_1(x)$', '$\theta_4(x)$', '$\theta_4(x)$', '$\theta_5(y)$', '$\theta_8(z)$'};
yLabels = {'$\theta_2(y)$', '$\theta_5(y)$', '$\theta_6(xz)$', '$\theta_6(xz)$', '$\theta_9(xy)$'};

% 绘制子图
for i = 1:5
    ax = axes('Position', positions(i, :)); 
    hold on;
    
    % 生成数据
    if i==1
        x1 = normrnd(A(1).mu2, A(1).sigma2, 2000, 1);
        y1 = normrnd(A(2).mu2, A(2).sigma2, 2000, 1);
    elseif i==2
        x1 = normrnd(A(3).mu2, A(3).sigma2, 2000, 1);
        y1 = normrnd(A(4).mu2, A(4).sigma2, 2000, 1);
    elseif i==3
        x1 = normrnd(A(3).mu2, A(3).sigma2, 2000, 1);
        y1 = normrnd(A(5).mu2, A(5).sigma2, 2000, 1);
    elseif i==4
        x1 = normrnd(A(4).mu2, A(4).sigma2, 2000, 1);
        y1 = normrnd(A(5).mu2, A(5).sigma2, 2000, 1);
    elseif i==5
        x1 = normrnd(A(6).mu2, A(6).sigma2, 2000, 1);
        y1 = normrnd(A(7).mu2, A(7).sigma2, 2000, 1);
    end

    % 绘制散点图
    data1 = [x1, y1];
    density_2D_1 = density2D_KD(data1(:, 1:2), 0.01);
    scatter(data1(:, 1), data1(:, 2), 18, density_2D_1, 'filled');
    
    % 设置坐标轴标签
    xlabel(xLabels{i}, 'FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal', 'Interpreter', 'latex');
    ylabel(yLabels{i}, 'FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal', 'Interpreter', 'latex');
    
    colormap(map);
    colorbar;

    % 设置刻度
    xTicks = linspace(min(data1(:, 1)), max(data1(:, 1)), 3);
    yTicks = linspace(min(data1(:, 2)), max(data1(:, 2)), 3);
    xTickLabels = formatTickLabels(xTicks);
    yTickLabels = formatTickLabels(yTicks);
    
    set(gca, 'Box', 'off', 'LineWidth', 1, ...
             'XGrid', 'off', 'YGrid', 'off', ...
             'TickDir', 'out', 'TickLength', [.005 .005], ...
             'XMinorTick', 'off', 'YMinorTick', 'off', ...
             'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
             'XTick', xTicks, 'XTickLabel', xTickLabels, ...
             'XLim', [min(data1(:, 1))-0.01 max(data1(:, 1))+0.01], ...
             'YTick', yTicks, 'YTickLabel', yTickLabels, ...
             'YLim', [min(data1(:, 2))-0.001 max(data1(:, 2))+0.001], ...
             'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'normal');
  
end
%% 图片输出
figW = figureWidth;
figH = figureHeight;
set(figureHandle,'PaperUnits',figureUnits);
set(figureHandle,'PaperPosition',[0 0 figW figH]);
fileout = '漂移项系数联合后验分布';
print(figureHandle,[fileout,'.png'],'-r900','-dpng');
%% 画空间误差图
map = colormap(nclCM(232)); 
map = flipud(map);
close all;
figureUnits = 'centimeters';
figureWidth = 28; 
figureHeight = 20;
figureHandle = figure;
set(gcf, 'Units', figureUnits, 'Position', [2 5 figureWidth figureHeight], 'Color', [1 1 1]);

% ========== 第一个子图 ==========
subplot(3, 3, 1); 
hold on;
fun1 = @(a,b) -a + b;
fun_1 = @(a,b) Drift2_X(2)*a + Drift2_X(4)*b;
fmesh(fun1, 'LineWidth', 1, 'MeshDensity', 25);
hXLabel1 = xlabel('$x$', 'Interpreter', 'latex');
hYLabel1 = ylabel('$y$', 'Interpreter', 'latex');
hZLabel1 = zlabel('$f_1(x,y)$', 'Interpreter', 'latex');
view(-20, 10);
colormap(map);
colorbar;
set(gca, 'Box', 'off', 'LineWidth', 1, 'GridLineStyle', '-', ...
    'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', ...
    'TickDir', 'out', 'TickLength', [.015 .015], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'ZMinorTick', 'off', ...
    'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], 'ZColor', [.1 .1 .1]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 13);
set([hXLabel1, hYLabel1], 'FontName', 'Times New Roman', 'FontSize', 17);
set(hZLabel1, 'FontName', 'SimSun', 'FontSize', 17);


% ========== 第二个子图 ==========
subplot(3, 3, 2); 
hold on;
fmesh(fun_1, 'LineWidth', 1, 'MeshDensity', 25);
hXLabel2 = xlabel('$x$', 'Interpreter', 'latex');
hYLabel2 = ylabel('$y$', 'Interpreter', 'latex');
hZLabel2= zlabel('$\hat{f}_1(x,y)$', 'Interpreter', 'latex');
view(-20, 10);
colormap(map);
colorbar;
set(gca, 'Box', 'off', 'LineWidth', 1, 'GridLineStyle', '-', ...
    'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', ...
    'TickDir', 'out', 'TickLength', [.015 .015], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'ZMinorTick', 'off', ...
    'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], 'ZColor', [.1 .1 .1]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 13);
set([hXLabel2, hYLabel2], 'FontName', 'Times New Roman', 'FontSize', 17);
set(hZLabel2, 'FontName', 'SimSun', 'FontSize', 17);


% ========== 第三个子图 ==========
subplot(3, 3, 3);
hold on;
fun_error = @(a,b) abs(fun1(a,b) - fun_1(a,b));
fmesh(fun_error, 'LineWidth', 1, 'MeshDensity', 25);
hXLabel3 = xlabel('$x$', 'Interpreter', 'latex');
hYLabel3 = ylabel('$y$', 'Interpreter', 'latex');
hZLabel3= zlabel('$|f_1 -\hat{f}_1|$', 'Interpreter', 'latex');
view(-20, 10);
colormap(map);
colorbar;
set(gca, 'Box', 'off', 'LineWidth', 1, 'GridLineStyle', '-', ...
    'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', ...
    'TickDir', 'out', 'TickLength', [.015 .015], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'ZMinorTick', 'off', ...
    'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], 'ZColor', [.1 .1 .1], ...
    'ZTick', 0:0.02:0.04, 'ZLim', [0 0.04]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 13);
set([hXLabel3, hYLabel3], 'FontName', 'Times New Roman', 'FontSize', 17);
set(hZLabel3, 'FontName', 'SimSun', 'FontSize', 17);

% ========== 第四个子图 ==========
subplot(3, 3, 4); 
hold on;
x_ax = -2.5:0.1:2.5;
[X_pic, Y_pic, Z_pic] = meshgrid(x_ax, x_ax, x_ax);
fun2 =  0.5*X_pic - Y_pic - X_pic.*Z_pic;
fun_2 =  Drift2_Y(2)*X_pic + Drift2_Y(4)*Y_pic + Drift2_Y(12)*X_pic.*Z_pic;
scatter3(X_pic(:), Y_pic(:), Z_pic(:), 1, fun2(:), 'filled')
hXLabel4 = xlabel('$x$', 'Interpreter', 'latex');
hYLabel4 = ylabel('$y$', 'Interpreter', 'latex');
hZLabel4= zlabel('$z$', 'Interpreter', 'latex');
view(-40, 12);
colormap(map);
colorbar;
set(gca, 'Box', 'off', 'LineWidth', 1, 'GridLineStyle', '-', ...
    'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', ...
    'TickDir', 'out', 'TickLength', [.015 .015], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'ZMinorTick', 'off', ...
    'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], 'ZColor', [.1 .1 .1]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 13);
set([hXLabel4, hYLabel4,hZLabel4], 'FontName', 'Times New Roman', 'FontSize', 17);

% ========== 第五个子图 ==========
subplot(3, 3, 5); 
hold on;
scatter3(X_pic(:), Y_pic(:), Z_pic(:), 1, fun_2(:), 'filled')
hXLabel5 = xlabel('$x$', 'Interpreter', 'latex');
hYLabel5 = ylabel('$y$', 'Interpreter', 'latex');
hZLabel5= zlabel('$z$', 'Interpreter', 'latex');
view(-40, 12);
colormap(map);
colorbar;
set(gca, 'Box', 'off', 'LineWidth', 1, 'GridLineStyle', '-', ...
    'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', ...
    'TickDir', 'out', 'TickLength', [.015 .015], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'ZMinorTick', 'off', ...
    'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], 'ZColor', [.1 .1 .1]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 13);
set([hXLabel5, hYLabel5,hZLabel5], 'FontName', 'Times New Roman', 'FontSize', 17);


% ========== 第六个子图 ==========
subplot(3, 3, 6); 
hold on;
fun_error2 = abs(fun2 - fun_2);
scatter3(X_pic(:), Y_pic(:), Z_pic(:), 1, fun_error2(:), 'filled')
hXLabel6 = xlabel('$x$', 'Interpreter', 'latex');
hYLabel6 = ylabel('$y$', 'Interpreter', 'latex');
hZLabel6= zlabel('$z$', 'Interpreter', 'latex');
view(-40, 12);
colormap(map);
colorbar;
set(gca, 'Box', 'off', 'LineWidth', 1, 'GridLineStyle', '-', ...
    'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', ...
    'TickDir', 'out', 'TickLength', [.015 .015], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'ZMinorTick', 'off', ...
    'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], 'ZColor', [.1 .1 .1]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 13);
set([hXLabel6, hYLabel6,hZLabel6], 'FontName', 'Times New Roman', 'FontSize', 17);

% ========== 第七个子图 ==========
subplot(3, 3, 7); 
hold on;
fun3 =  (-8/3)*Z_pic + X_pic.*Y_pic ;
fun_3 =  Drift2_Z(6)*Z_pic + Drift2_Z(8)*X_pic.*Y_pic;
scatter3(X_pic(:), Y_pic(:), Z_pic(:), 1, fun3(:), 'filled')
hXLabel7 = xlabel('$x$', 'Interpreter', 'latex');
hYLabel7 = ylabel('$y$', 'Interpreter', 'latex');
hZLabel7= zlabel('$z$', 'Interpreter', 'latex');
view(-40, 12);
colormap(map);
colorbar;
set(gca, 'Box', 'off', 'LineWidth', 1, 'GridLineStyle', '-', ...
    'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', ...
    'TickDir', 'out', 'TickLength', [.015 .015], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'ZMinorTick', 'off', ...
    'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], 'ZColor', [.1 .1 .1]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 13);
set([hXLabel7, hYLabel7,hZLabel7], 'FontName', 'Times New Roman', 'FontSize', 17);

% ========== 第八个子图 ==========
subplot(3, 3, 8); 
hold on;
scatter3(X_pic(:), Y_pic(:), Z_pic(:), 1, fun_3(:), 'filled')
hXLabel8 = xlabel('$x$', 'Interpreter', 'latex');
hYLabel8 = ylabel('$y$', 'Interpreter', 'latex');
hZLabel8= zlabel('$z$', 'Interpreter', 'latex');
view(-40, 12);
colormap(map);
colorbar;
set(gca, 'Box', 'off', 'LineWidth', 1, 'GridLineStyle', '-', ...
    'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', ...
    'TickDir', 'out', 'TickLength', [.015 .015], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'ZMinorTick', 'off', ...
    'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], 'ZColor', [.1 .1 .1]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 13);
set([hXLabel8, hYLabel8,hZLabel8], 'FontName', 'Times New Roman', 'FontSize', 17);

% ========== 第九个子图 ==========
subplot(3, 3, 9);
hold on;
fun_error3 = abs(fun3 - fun_3);
scatter3(X_pic(:), Y_pic(:), Z_pic(:), 1, fun_error3(:), 'filled')
hXLabel9 = xlabel('$x$', 'Interpreter', 'latex');
hYLabel9 = ylabel('$y$', 'Interpreter', 'latex');
hZLabel9= zlabel('$z$', 'Interpreter', 'latex');
view(-40, 12);
colormap(map);
colorbar;
set(gca, 'Box', 'off', 'LineWidth', 1, 'GridLineStyle', '-', ...
    'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', ...
    'TickDir', 'out', 'TickLength', [.015 .015], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'ZMinorTick', 'off', ...
    'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], 'ZColor', [.1 .1 .1]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 13);
set([hXLabel9, hYLabel9,hZLabel9], 'FontName', 'Times New Roman', 'FontSize', 17);

%% 图片输出
figW = figureWidth;
figH = figureHeight;
set(figureHandle,'PaperUnits',figureUnits);
set(figureHandle,'PaperPosition',[0 0 figW figH]);
fileout = '漂移项三维误差图像';
print(figureHandle,[fileout,'.png'],'-r100','-dpng');

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



