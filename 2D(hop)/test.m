clc
clear
close all
%
rng(100, 'twister');
alpha =0.5;
epsilon = 0.05;
num_N = 10000;
dt = 0.001;
T = dt * num_N;
c_alpha = alpha * gamma((1+alpha)/2) / (2^(1-alpha) * pi^(0.5) * gamma(1-alpha/2));
Lt1 = generate_Lt(alpha,epsilon,num_N,dt,c_alpha);
Lt2 = generate_Lt(alpha,epsilon,num_N,dt,c_alpha);

Xzero = 1;
X_simul = zeros(num_N + 1,1);
Yzero = -1;
Y_simul = zeros(num_N + 1,1);
X_simul(1) = Xzero;
Y_simul(1) = Yzero;
u = 1.5;
w = 0.5;
a = 1;
b = 1;
sigma2 = 4; 
for i = 2:num_N + 1
        X_simul(i) = X_simul(i-1) + ( u*X_simul(i-1) - a*X_simul(i-1).^3 - w*Y_simul(i-1) - b*Y_simul(i-1).^3 -a*X_simul(i-1)* Y_simul(i-1)^2 -b*X_simul(i-1).^2*Y_simul(i-1) ) * dt + sqrt(sigma2) * (Lt1(i) - Lt1(i-1));
        Y_simul(i) = Y_simul(i-1) + (w*X_simul(i-1) + b*X_simul(i-1).^3 + u*Y_simul(i-1) -a*Y_simul(i-1).^3 + b*X_simul(i-1)*Y_simul(i-1).^2 - a*X_simul(i-1).^2*Y_simul(i-1) ) * dt + sqrt(sigma2) * (Lt2(i) - Lt2(i-1));
end

any(isinf(Y_simul))   % 检查是否有Inf
%%
figure;
T = dt * num_N;
plot(0:dt:T, Y_simul, 'DisplayName', 'Y\_simul(t)', 'Color', 'b');
hold on;
plot(0:dt:T, X_simul, 'DisplayName', 'X\_simul(t)', 'Color', 'r');
xlabel('Time (t)');
ylabel('X(t)');
title('随机的');
grid on;
legend('show');

%%
X = X_simul(1:end-1);
Y = Y_simul(1:end-1);
FX_Y =(Y_simul(2:end)-Y_simul(1:end-1))./dt; 
FX_X =(X_simul(2:end)-X_simul(1:end-1))./dt; 

order = 3;  
[deg_funcs, deg_strs] = deg_func_xy(order);
num_funcs = length(deg_funcs);
Phi = zeros(length(X), num_funcs); 
for k = 1:num_funcs
    func = deg_funcs{k}; % 获取第 k 个函数句柄
    Phi(:, k) = func(X , Y); % 计算函数值并存储
end
Phi_old = Phi;

indices_FX_X = find(abs(FX_X) > 100);
indices_FX_Y = find(abs(FX_Y) > 100);
total_indices = unique([indices_FX_X; indices_FX_Y]);
FX_X(total_indices) = [];
FX_Y(total_indices) = [];
X(total_indices) = [];
Y(total_indices) = [];
Phi(total_indices, :) = []; 

[Drift1_X,Drift_var1_X] = MSBL(Phi, FX_X, 1, 1);
[Drift1_Y,Drift_var1_Y] = MSBL(Phi, FX_Y, 1, 1);
p = 0.5;
[~,Drift2_X,Drift_var2_X] = StepwiseRSBL(Phi, FX_X, 1, 1, p);
[~,Drift2_Y,Drift_var2_Y] = StepwiseRSBL(Phi, FX_Y, 1, 1, p);
a1 = (abs(Drift2_X(4))+abs(Drift2_X(9))+abs(Drift2_Y(7))+abs(Drift2_Y(11)))/4;
b1 = (abs(Drift2_X(7))+abs(Drift2_X(11))+abs(Drift2_Y(4))+abs(Drift2_Y(9)))/4;
u1 = (abs(Drift2_X(2))+abs(Drift2_Y(5)))/2;
w1 = (abs(Drift2_X(5))+abs(Drift2_Y(2)))/2;

var_alpha = 2 * dt * c_alpha /(2 - alpha);
XG =  X_simul(1:end-1);
YG = Y_simul(1:end-1);
GX1_X = ((X_simul(2:end,:)- X_simul(1:end-1,:)-Phi_old*Drift1_X*dt)).^2/ var_alpha;
GY1_Y = ((Y_simul(2:end,:) - Y_simul(1:end-1,:) - Phi_old*Drift1_Y*dt)).^2/ var_alpha;

GX2_X = ((X_simul(2:end,:)- X_simul(1:end-1,:)-(u1*X_simul(i-1) - a1*X_simul(i-1).^3 - w1*Y_simul(i-1) - b1*Y_simul(i-1).^3 - a1*X_simul(i-1)* Y_simul(i-1)^2 - b1*X_simul(i-1).^2*Y_simul(i-1) )*dt)).^2/ var_alpha;
GY2_Y = ((Y_simul(2:end,:) - Y_simul(1:end-1,:) - (w1*X_simul(i-1) + b1*X_simul(i-1).^3 + u1*Y_simul(i-1) - a1*Y_simul(i-1).^3 + b1*X_simul(i-1)*Y_simul(i-1).^2 - a1*X_simul(i-1).^2*Y_simul(i-1) )* dt)).^2/ var_alpha;

indices_GX_X2 = find(abs(GX2_X) >100);
indices_GY_Y2 = find(abs(GY2_Y) > 100);
total_indices2 = unique([indices_GX_X2; indices_GY_Y2]);
GX1_X(total_indices2) = [];
GY1_Y(total_indices2) = [];
GX2_X(total_indices2) = [];
GY2_Y(total_indices2) = [];
XG(total_indices2) = [];
YG(total_indices2) = [];
Phi_old(total_indices2, :) = []; 

Phi_new = [ones(length(GX2_X),1)]; 
p2 =0.1;
[~,Diff1_X,Diff_var1_X] = MSBL(Phi_new, GX1_X*10, 1, 1);
[~,Diff1_Y,Diff_var1_Y] = MSBL(Phi_new, GY1_Y*10, 1, 1);
[~,Diff2_X,Diff_var2_X] = StepwiseRSBL(Phi_new, GX2_X*10, 1, 1, p2);
[~,Diff2_Y,Diff_var2_Y] = StepwiseRSBL(Phi_new, GY2_Y*10, 1, 1, p2);
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
         'XTick', 0:3:10.1, ...
         'XLim', [0 10.1], ...
         'YTick', -2:1:2.5, ...
         'YLim', [-2 2.5]);

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
         'XTick', 0:3:10.1, ...
         'XLim', [0 10.1], ...
         'YTick', -2.5:2:2, ...
         'YLim', [-2.5 2]);
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

%% 画漂移项的函数库情况（sbl和rsbl）

labels = deg_strs;
C1 = addcolorplus(252); 
C2 = addcolorplus(119); 
data = Drift1_X;
data1 = Drift2_X;
data2 = Drift1_Y;
data3 = Drift2_Y;

close all
figureUnits = 'centimeters';% 图片尺寸设置（单位：厘米）
figureWidth = 28;
figureHeight = 18;
figureHandle = figure;% 窗口设置
set(gcf, 'Units', figureUnits, 'Position', [0 0 figureWidth figureHeight]);

% 子图布局参数
topMargin = 0.04;       % 上边距
bottomMargin = 0.15;    % 下边距（为标题留空间）
leftMargin = 0.06;
rightMargin = 0.04;
subplotGap = 0.16;      % 子图间距

% 计算子图高度
subplotHeight = (1 - topMargin - bottomMargin - subplotGap)/2;

% ====================== 第一个子图 ======================
subplot('Position', [leftMargin, bottomMargin+subplotHeight+subplotGap, 1-leftMargin-rightMargin, subplotHeight]); 
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
axis([0 16 -2 2]) % 调整坐标轴范围
for i = 1:length(data)
    [maxVal, maxIdx] = max([abs(data(i)), abs(data1(i))]);
    if maxIdx == 1
        labelPos = data(i);
    else
        labelPos = data1(i);
    end
    if maxVal >= 0.01 
        if labelPos > 0
            text(i, labelPos + 0.1, labels{i}, ... % y 方向增加 0.1 的偏移量
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
    'Location', 'northeast', ... % 图例位置在右上角
    'FontSize', 15, ...
    'Box', 'off','FontWeight' , 'normal', 'FontName', 'Times New Roman'); % 去掉图例边框
hXLabel = xlabel('Library functions');
hYLabel = ylabel('Coefficients (1st component)');
set(gca, 'Box', 'off', ... % 边框
    'LineWidth', 1, ... % 线宽
    'XGrid', 'off', 'YGrid', 'off', ... % 网格
    'TickDir', 'out', 'TickLength', [.01 .01], ... % 刻度
    'XMinorTick', 'off', 'YMinorTick', 'off', ... % 小刻度
    'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], 'FontName', 'Times New Roman', 'FontSize', 15,'FontWeight' , 'bold') % 坐标轴颜色
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
set( [hXLabel,hYLabel], 'FontSize', 15, 'FontName', 'Times New Roman','FontWeight' , 'normal')

subplot('Position', [leftMargin, bottomMargin, 1-leftMargin-rightMargin, subplotHeight]); 
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
axis([0 16 -2 2]) % 调整坐标轴范围
for i = 1:length(data2)
    [maxVal, maxIdx] = max([abs(data2(i)), abs(data3(i))]);
    if maxIdx == 1
        labelPos = data2(i);
    else
        labelPos = data3(i);
    end
    if maxVal >= 0.01 
        if labelPos > 0
            text(i, labelPos + 0.1, labels{i}, ... % y 方向增加 0.1 的偏移量
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
    'Location', 'northeast', ... % 图例位置在右上角
    'FontSize', 15, ...
    'Box', 'off','FontWeight' , 'normal', 'FontName', 'Times New Roman'); % 去掉图例边框
hXLabel = xlabel('Library functions');
hYLabel = ylabel('Coefficients (2nd component)');
ax = gca;
ax.XLabel.Position(2) = ax.XLabel.Position(2) - 0.5; % 下移x轴标签

set(gca, 'Box', 'off', ... % 边框
    'LineWidth', 1, ... % 线宽
    'XGrid', 'off', 'YGrid', 'off', ... % 网格
    'TickDir', 'out', 'TickLength', [.01 .01], ... % 刻度
    'XMinorTick', 'off', 'YMinorTick', 'off', ... % 小刻度
    'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], 'FontName', 'Times New Roman', 'FontSize', 15,'FontWeight' , 'bold') % 坐标轴颜色
set(gca, 'Position', [ax.Position(1), ax.Position(2)+0.03, ax.Position(3), ax.Position(4)-0.03]);

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

%% 5个参数的平均估计
C2 = addcolorplus(119); 
close all;
% 计算均值
a1 = (abs(Drift2_X(4))+abs(Drift2_X(9))+abs(Drift2_Y(7))+abs(Drift2_Y(11)))/4;
b1 = (abs(Drift2_X(7))+abs(Drift2_X(11))+abs(Drift2_Y(4))+abs(Drift2_Y(9)))/4;
u1 = (abs(Drift2_X(2))+abs(Drift2_Y(5)))/2;
w1 = (abs(Drift2_X(5))+abs(Drift2_Y(2)))/2;
sigma1 = (abs(Diff2_X(1))+abs(Diff2_Y(1)))/2;

% 计算方差
a1_var = (Drift_var2_X(4)+Drift_var2_X(9)+Drift_var2_Y(7)+Drift_var2_Y(11))/16;
b1_var = (Drift_var2_X(7)+Drift_var2_X(11)+Drift_var2_Y(4)+Drift_var2_Y(9))/16;
u1_var = (Drift_var2_X(2)+Drift_var2_Y(5))/4;
w1_var = (Drift_var2_X(5)+Drift_var2_Y(2))/4;
sigma1_var = (Diff_var2_X(1)+Diff_var2_Y(1))/4;

data = struct();

data(1).mu = a1; 
data(1).sigma = sqrt(a1_var);
data(2).mu = b1; 
data(2).sigma = sqrt(b1_var);
data(3).mu = u1 ; 
data(3).sigma = sqrt(u1_var);
data(4).mu = w1; 
data(4).sigma = sqrt(w1_var);
data(5).mu = sigma1; 
data(5).sigma = sqrt(sigma1_var);

figureUnits = 'centimeters';
figureWidth = 28; % 宽度增加以适应 4 列
figureHeight = 18; % 高度增加以适应 3 行
figureHandle = figure;
set(gcf, 'Units', figureUnits, 'Position', [2 5 figureWidth figureHeight]);
titles = {'a', 'b', '\mu', 'w', '\sigma'};

% 子图布局参数
topMargin = 0.01;       % 上边距
bottomMargin = 0.15;    % 下边距（为标题留空间）
leftMargin = 0.06;
rightMargin = 0.04;
subplotGap = 0.16;      % 子图间距

% 定义子图的位置
positions = [
    0.07, 0.65, 0.27, 0.32; % 第一张图
    0.39, 0.65, 0.27, 0.32; % 第二张图
    0.7, 0.65, 0.27, 0.32; % 第三张图
    0.21, 0.15, 0.27, 0.32; % 第四张图（居中，第一张和第二张之间）
    0.55, 0.15, 0.27, 0.32; % 第五张图（居中，第二张和第三张之间）
];

% 绘制子图
for i = 1:5
    ax = axes('Position', positions(i, :)); % 设置子图位置
    hold on;
    x = linspace(data(i).mu - 4 * data(i).sigma, data(i).mu + 4 * data(i).sigma, 3000); 
    pdf = normpdf(x, data(i).mu, data(i).sigma); 
    area(x, pdf, 'LineWidth', 1.5, 'FaceColor', C2, 'EdgeColor', C2, ...
        'FaceAlpha', .3, 'EdgeAlpha', 1);
    line([data(i).mu data(i).mu], [0 normpdf(data(i).mu, data(i).mu, data(i).sigma)], ...
        'LineStyle', '--', 'Color', C2, 'LineWidth', 1.5);
    text(data(i).mu, normpdf(data(i).mu, data(i).mu, data(i).sigma)+2 , ...
            sprintf('N(%.4f, %.4f)', data(i).mu, data(i).sigma^2), ...
            'HorizontalAlignment', 'center', 'FontSize', 15, 'Color', 'k', 'Interpreter', 'tex', 'FontName', 'Times New Roman', 'FontWeight', 'normal');
    set(gca, 'Box', 'off', 'XGrid', 'off', 'YGrid', 'off', ...
             'TickDir', 'out', 'TickLength', [.01 .01], ...
             'XMinorTick', 'off', 'YMinorTick', 'off', ...
             'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
             'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'normal', ...
             'YTick', 0:round((max(pdf)+2)/3):max(pdf)+5,'YLim', [0 max(pdf)+5]);
      xlabel(sprintf('$\\hat{%s}$', titles{i}), 'FontSize', 17, ...
           'FontName', 'Times New Roman', ...
           'FontWeight', 'normal', ...
           'Interpreter', 'latex');
    
    if i == 1 || i == 4
        ylabel('PDF', 'FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
    end
end

set(gcf, 'Color', [1 1 1]);
%% 图片输出
figW = figureWidth;
figH = figureHeight;
set(figureHandle,'PaperUnits',figureUnits);
set(figureHandle,'PaperPosition',[0 0 figW figH]);
fileout = '五个参数的分布';
print(figureHandle,[fileout,'.png'],'-r600','-dpng');


%% 画漂移项联合后验分布(一维)
map = colormap(nclCM(232)); % color包里选颜色
map = flipud(map);
close all;
A = struct();
index_A =[2,4,5,7,9,11];
for i = 1:length(index_A)
    A(i).mu = abs(Drift2_X(index_A(i))); 
    A(i).sigma = sqrt(Drift_var2_X(index_A(i)));
end
figureUnits = 'centimeters';
figureWidth = 28; 
figureHeight = 22;
figureHandle = figure;
set(gcf, 'Units', figureUnits, 'Position', [0 0 figureWidth figureHeight]);
% 子图布局参数


subplot(3, 3, 1); 
hold on;
x1 = normrnd(A(1).mu, A(1).sigma, 3000, 1);
y1 = normrnd(A(2).mu, A(2).sigma, 3000, 1);
data1 = [x1, y1];
density_2D_1 = density2D_KD(data1(:, 1:2), 0.01);
scatter(data1(:, 1), data1(:, 2), 15, density_2D_1, 'filled');
xlabel('$\mu(x)$', 'Interpreter', 'latex', 'FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
ylabel('$a(x^3)$', 'Interpreter', 'latex', 'FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
colormap(map);
colorbar;   
set(gca, 'Box', 'off','LineWidth', 1,'XGrid', 'off', 'YGrid', 'off','TickDir', 'out', 'TickLength', [.005 .005],'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', min(data1(:, 1))-0.05:0.5:max(data1(:, 1)+0.05),'XLim', [min(data1(:, 1))-0.05 max(data1(:, 1))+0.05], ...
         'YTick', min(data1(:, 2))-0.05:0.2:max(data1(:, 2)+0.05),'YLim', [min(data1(:, 2))-0.05 max(data1(:, 2))+0.05], ...
         'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'normal');
xticks = get(gca, 'XTick');
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
yticks = get(gca, 'YTick');
yticklabels = arrayfun(@(y) sprintf('%.2f', y), yticks, 'UniformOutput', false);
set(gca, 'XTickLabel', xticklabels,'YTickLabel', yticklabels); 

subplot(3, 3, 2); 
hold on;
x1 = normrnd(A(1).mu, A(1).sigma, 3000, 1);
y1 = normrnd(A(3).mu, A(3).sigma, 3000, 1);
data1 = [x1, y1];
density_2D_1 = density2D_KD(data1(:, 1:2), 0.01);
scatter(data1(:, 1), data1(:, 2), 15, density_2D_1, 'filled');
xlabel('$\mu(x)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
ylabel('$\omega(y)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
colormap(map);
colorbar;   
set(gca, 'Box', 'off','LineWidth', 1,'XGrid', 'off', 'YGrid', 'off','TickDir', 'out', 'TickLength', [.005 .005],'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', min(data1(:, 1))-0.05:0.5:max(data1(:, 1)+0.05), ...
         'XLim', [min(data1(:, 1))-0.05 max(data1(:, 1))+0.05], ...
         'YTick', min(data1(:, 2))-0.05:0.3:max(data1(:, 2)+0.05), ...
         'YLim', [min(data1(:, 2))-0.05 max(data1(:, 2))+0.05], ...
         'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'normal');
xticks = get(gca, 'XTick');
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
yticks = get(gca, 'YTick');
yticklabels = arrayfun(@(y) sprintf('%.2f', y), yticks, 'UniformOutput', false);
set(gca, 'XTickLabel', xticklabels,'YTickLabel', yticklabels); 

subplot(3, 3, 3); 
hold on;
x1 = normrnd(A(1).mu, A(1).sigma, 3000, 1);
y1 = normrnd(A(4).mu, A(4).sigma, 3000, 1);
data1 = [x1, y1];
density_2D_1 = density2D_KD(data1(:, 1:2), 0.01);
scatter(data1(:, 1), data1(:, 2), 15, density_2D_1, 'filled');
xlabel('$\mu(x)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
ylabel('$b(y^3)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
colormap(map);
colorbar;   
set(gca, 'Box', 'off','LineWidth', 1,'XGrid', 'off', 'YGrid', 'off','TickDir', 'out', 'TickLength', [.005 .005],'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', min(data1(:, 1))-0.05:0.5:max(data1(:, 1)+0.05), ...
         'XLim', [min(data1(:, 1))-0.05 max(data1(:, 1))+0.05], ...
         'YTick', min(data1(:, 2))-0.05:0.15:max(data1(:, 2)+0.05), ...
         'YLim', [min(data1(:, 2))-0.05 max(data1(:, 2))+0.05], ...
         'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'normal');
xticks = get(gca, 'XTick');
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
yticks = get(gca, 'YTick');
yticklabels = arrayfun(@(y) sprintf('%.2f', y), yticks, 'UniformOutput', false);
set(gca, 'XTickLabel', xticklabels,'YTickLabel', yticklabels); 

subplot(3, 3, 4); 
hold on;
x1 = normrnd(A(2).mu, A(2).sigma, 3000, 1);
y1 = normrnd(A(3).mu, A(3).sigma, 3000, 1);
data1 = [x1, y1];
density_2D_1 = density2D_KD(data1(:, 1:2), 0.01);
scatter(data1(:, 1), data1(:, 2), 15, density_2D_1, 'filled');
xlabel('$a(x^3)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
ylabel('$\omega(y)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
colormap(map);
colorbar;   
set(gca, 'Box', 'off','LineWidth', 1,'XGrid', 'off', 'YGrid', 'off','TickDir', 'out', 'TickLength', [.005 .005],'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', min(data1(:, 1))-0.05:0.25:max(data1(:, 1)+0.05), ...
         'XLim', [min(data1(:, 1))-0.05 max(data1(:, 1))+0.05], ...
         'YTick', min(data1(:, 2))-0.05:0.3:max(data1(:, 2)+0.05), ...
         'YLim', [min(data1(:, 2))-0.05 max(data1(:, 2))+0.05], ...
         'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'normal');

xticks = get(gca, 'XTick');
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
yticks = get(gca, 'YTick');
yticklabels = arrayfun(@(y) sprintf('%.2f', y), yticks, 'UniformOutput', false);
set(gca, 'XTickLabel', xticklabels,'YTickLabel', yticklabels); 

subplot(3, 3, 5); 
hold on;
x1 = normrnd(A(2).mu, A(2).sigma, 3000, 1);
y1 = normrnd(A(4).mu, A(4).sigma, 3000, 1);
data1 = [x1, y1];
density_2D_1 = density2D_KD(data1(:, 1:2), 0.01);
scatter(data1(:, 1), data1(:, 2), 15, density_2D_1, 'filled');
xlabel('$a(x^3)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
ylabel('$b(y^3)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
colormap(map);
colorbar;   
set(gca, 'Box', 'off','LineWidth', 1,'XGrid', 'off', 'YGrid', 'off','TickDir', 'out', 'TickLength', [.005 .005],'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', min(data1(:, 1))-0.05:0.25:max(data1(:, 1)+0.05), ...
         'XLim', [min(data1(:, 1))-0.05 max(data1(:, 1))+0.05], ...
         'YTick', min(data1(:, 2))-0.05:0.15:max(data1(:, 2)+0.05), ...
         'YLim', [min(data1(:, 2))-0.05 max(data1(:, 2))+0.05], ...
         'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'normal');

xticks = get(gca, 'XTick');
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
yticks = get(gca, 'YTick');
yticklabels = arrayfun(@(y) sprintf('%.2f', y), yticks, 'UniformOutput', false);
set(gca, 'XTickLabel', xticklabels,'YTickLabel', yticklabels); 

subplot(3, 3, 6); 
hold on;
x1 = normrnd(A(2).mu, A(2).sigma, 3000, 1);
y1 = normrnd(A(5).mu, A(5).sigma, 3000, 1);
data1 = [x1, y1];
density_2D_1 = density2D_KD(data1(:, 1:2), 0.01);
scatter(data1(:, 1), data1(:, 2), 15, density_2D_1, 'filled');
xlabel('$a(x^3)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
ylabel('$a(xy^2)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
colormap(map);
colorbar;   
set(gca, 'Box', 'off','LineWidth', 1,'XGrid', 'off', 'YGrid', 'off','TickDir', 'out', 'TickLength', [.005 .005],'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', min(data1(:, 1))-0.05:0.3:max(data1(:, 1)+0.05), ...
         'XLim', [min(data1(:, 1))-0.05 max(data1(:, 1))+0.05], ...
         'YTick', min(data1(:, 2))-0.05:0.25:max(data1(:, 2)+0.05), ...
         'YLim', [min(data1(:, 2))-0.05 max(data1(:, 2))+0.05], ...
         'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'normal');

xticks = get(gca, 'XTick');
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
yticks = get(gca, 'YTick');
yticklabels = arrayfun(@(y) sprintf('%.2f', y), yticks, 'UniformOutput', false);
set(gca, 'XTickLabel', xticklabels,'YTickLabel', yticklabels); 

subplot(3, 3, 7); 
hold on;
x1 = normrnd(A(3).mu, A(3).sigma, 3000, 1);
y1 = normrnd(A(4).mu, A(4).sigma, 3000, 1);
data1 = [x1, y1];
density_2D_1 = density2D_KD(data1(:, 1:2), 0.01);
scatter(data1(:, 1), data1(:, 2), 15, density_2D_1, 'filled');
xlabel('$\omega(y)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
ylabel('$b(y^3)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
colormap(map);
colorbar;   
set(gca, 'Box', 'off','LineWidth', 1,'XGrid', 'off', 'YGrid', 'off','TickDir', 'out', 'TickLength', [.005 .005],'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', min(data1(:, 1))-0.05:0.4:max(data1(:, 1)+0.05), ...
         'XLim', [min(data1(:, 1))-0.05 max(data1(:, 1))+0.05], ...
         'YTick', min(data1(:, 2))-0.05:0.15:max(data1(:, 2)+0.05), ...
         'YLim', [min(data1(:, 2))-0.05 max(data1(:, 2))+0.05], ...
         'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'normal');

xticks = get(gca, 'XTick');
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
yticks = get(gca, 'YTick');
yticklabels = arrayfun(@(y) sprintf('%.2f', y), yticks, 'UniformOutput', false);
set(gca, 'XTickLabel', xticklabels,'YTickLabel', yticklabels); 

subplot(3, 3, 8); 
hold on;
x1 = normrnd(A(3).mu, A(3).sigma, 3000, 1);
y1 = normrnd(A(5).mu, A(5).sigma, 3000, 1);
data1 = [x1, y1];
density_2D_1 = density2D_KD(data1(:, 1:2), 0.01);
scatter(data1(:, 1), data1(:, 2), 15, density_2D_1, 'filled');

xlabel('$\omega(y)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
ylabel('$a(xy^2)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
colormap(map);
colorbar;   
set(gca, 'Box', 'off','LineWidth', 1,'XGrid', 'off', 'YGrid', 'off','TickDir', 'out', 'TickLength', [.005 .005],'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', min(data1(:, 1))-0.05:0.5:max(data1(:, 1)+0.05), ...
         'XLim', [min(data1(:, 1))-0.05 max(data1(:, 1))+0.05], ...
         'YTick', min(data1(:, 2))-0.05:0.25:max(data1(:, 2)+0.05), ...
         'YLim', [min(data1(:, 2))-0.05 max(data1(:, 2))+0.05], ...
         'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'normal');
xticks = get(gca, 'XTick');
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
yticks = get(gca, 'YTick');
yticklabels = arrayfun(@(y) sprintf('%.2f', y), yticks, 'UniformOutput', false);
set(gca, 'XTickLabel', xticklabels,'YTickLabel', yticklabels); 

subplot(3, 3, 9); 
hold on;
x1 = normrnd(A(4).mu, A(4).sigma, 3000, 1);
y1 = normrnd(A(6).mu, A(6).sigma, 3000, 1);
data1 = [x1, y1];
density_2D_1 = density2D_KD(data1(:, 1:2), 0.01);
scatter(data1(:, 1), data1(:, 2), 15, density_2D_1, 'filled');

xlabel('$b(y^3)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
ylabel('$b(x^2y)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
colormap(map);
colorbar;   
set(gca, 'Box', 'off','LineWidth', 1,'XGrid', 'off', 'YGrid', 'off','TickDir', 'out', 'TickLength', [.005 .005],'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', min(data1(:, 1))-0.05:0.25:max(data1(:, 1)+0.05), ...
         'XLim', [min(data1(:, 1))-0.05 max(data1(:, 1))+0.05], ...
         'YTick', min(data1(:, 2))-0.05:0.2:max(data1(:, 2)+0.05), ...
         'YLim', [min(data1(:, 2))-0.05 max(data1(:, 2))+0.05], ...
         'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'normal');
xticks = get(gca, 'XTick');
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
yticks = get(gca, 'YTick');
yticklabels = arrayfun(@(y) sprintf('%.2f', y), yticks, 'UniformOutput', false);
set(gca, 'XTickLabel', xticklabels,'YTickLabel', yticklabels); 

set(gcf, 'Color', [1 1 1]);
%% 图片输出
figW = figureWidth;
figH = figureHeight;
set(figureHandle,'PaperUnits',figureUnits);
set(figureHandle,'PaperPosition',[0 0 figW figH]);
fileout = '(一维)漂移项系数联合后验分布';
print(figureHandle,[fileout,'.png'],'-r700','-dpng');

%% 画漂移项联合后验分布(二维)
map = colormap(nclCM(232)); % color包里选颜色
map = flipud(map);
close all;
A = struct();
index_A =[2,4,5,7,9,11];
for i = 1:length(index_A)
    A(i).mu = abs(Drift2_Y(index_A(i))); 
    A(i).sigma = sqrt(Drift_var2_Y(index_A(i)));
end
figureUnits = 'centimeters';
figureWidth = 28; 
figureHeight = 22;
figureHandle = figure;
set(gcf, 'Units', figureUnits, 'Position', [0 0 figureWidth figureHeight]);

figureSize = get(gcf, 'Position');

subplot(3, 3, 1); 
hold on;
x1 = normrnd(A(1).mu, A(1).sigma, 3000, 1);
y1 = normrnd(A(2).mu, A(2).sigma, 3000, 1);
data1 = [x1, y1];
density_2D_1 = density2D_KD(data1(:, 1:2), 0.01);
scatter(data1(:, 1), data1(:, 2), 15, density_2D_1, 'filled');
xlabel('$\mu(x)$', 'Interpreter', 'latex', 'FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
ylabel('$a(x^3)$', 'Interpreter', 'latex', 'FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
colormap(map);
colorbar;   
set(gca, 'Box', 'off','LineWidth', 1,'XGrid', 'off', 'YGrid', 'off','TickDir', 'out', 'TickLength', [.005 .005],'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', min(data1(:, 1))-0.05:0.5:max(data1(:, 1)+0.05),'XLim', [min(data1(:, 1))-0.05 max(data1(:, 1))+0.05], ...
         'YTick', min(data1(:, 2))-0.05:0.2:max(data1(:, 2)+0.05),'YLim', [min(data1(:, 2))-0.05 max(data1(:, 2))+0.05], ...
         'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'normal');
xticks = get(gca, 'XTick');
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
yticks = get(gca, 'YTick');
yticklabels = arrayfun(@(y) sprintf('%.2f', y), yticks, 'UniformOutput', false);
set(gca, 'XTickLabel', xticklabels,'YTickLabel', yticklabels); 

subplot(3, 3, 2); 
hold on;
x1 = normrnd(A(1).mu, A(1).sigma, 3000, 1);
y1 = normrnd(A(3).mu, A(3).sigma, 3000, 1);
data1 = [x1, y1];
density_2D_1 = density2D_KD(data1(:, 1:2), 0.01);
scatter(data1(:, 1), data1(:, 2), 15, density_2D_1, 'filled');
xlabel('$\mu(x)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
ylabel('$\omega(y)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
colormap(map);
colorbar;   
set(gca, 'Box', 'off','LineWidth', 1,'XGrid', 'off', 'YGrid', 'off','TickDir', 'out', 'TickLength', [.005 .005],'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', min(data1(:, 1))-0.05:0.5:max(data1(:, 1)+0.05), ...
         'XLim', [min(data1(:, 1))-0.05 max(data1(:, 1))+0.05], ...
         'YTick', min(data1(:, 2))-0.05:0.3:max(data1(:, 2)+0.05), ...
         'YLim', [min(data1(:, 2))-0.05 max(data1(:, 2))+0.05], ...
         'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'normal');
xticks = get(gca, 'XTick');
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
yticks = get(gca, 'YTick');
yticklabels = arrayfun(@(y) sprintf('%.2f', y), yticks, 'UniformOutput', false);
set(gca, 'XTickLabel', xticklabels,'YTickLabel', yticklabels); 

subplot(3, 3, 3); 
hold on;
x1 = normrnd(A(1).mu, A(1).sigma, 3000, 1);
y1 = normrnd(A(4).mu, A(4).sigma, 3000, 1);
data1 = [x1, y1];
density_2D_1 = density2D_KD(data1(:, 1:2), 0.01);
scatter(data1(:, 1), data1(:, 2), 15, density_2D_1, 'filled');
xlabel('$\mu(x)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
ylabel('$b(y^3)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
colormap(map);
colorbar;   
set(gca, 'Box', 'off','LineWidth', 1,'XGrid', 'off', 'YGrid', 'off','TickDir', 'out', 'TickLength', [.005 .005],'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', min(data1(:, 1))-0.05:0.5:max(data1(:, 1)+0.05), ...
         'XLim', [min(data1(:, 1))-0.05 max(data1(:, 1))+0.05], ...
         'YTick', min(data1(:, 2))-0.05:0.15:max(data1(:, 2)+0.05), ...
         'YLim', [min(data1(:, 2))-0.05 max(data1(:, 2))+0.05], ...
         'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'normal');
xticks = get(gca, 'XTick');
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
yticks = get(gca, 'YTick');
yticklabels = arrayfun(@(y) sprintf('%.2f', y), yticks, 'UniformOutput', false);
set(gca, 'XTickLabel', xticklabels,'YTickLabel', yticklabels); 

subplot(3, 3, 4); 
hold on;
x1 = normrnd(A(2).mu, A(2).sigma, 3000, 1);
y1 = normrnd(A(3).mu, A(3).sigma, 3000, 1);
data1 = [x1, y1];
density_2D_1 = density2D_KD(data1(:, 1:2), 0.01);
scatter(data1(:, 1), data1(:, 2), 15, density_2D_1, 'filled');
xlabel('$a(x^3)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
ylabel('$\omega(y)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
colormap(map);
colorbar;   
set(gca, 'Box', 'off','LineWidth', 1,'XGrid', 'off', 'YGrid', 'off','TickDir', 'out', 'TickLength', [.005 .005],'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', min(data1(:, 1))-0.05:0.25:max(data1(:, 1)+0.05), ...
         'XLim', [min(data1(:, 1))-0.05 max(data1(:, 1))+0.05], ...
         'YTick', min(data1(:, 2))-0.05:0.3:max(data1(:, 2)+0.05), ...
         'YLim', [min(data1(:, 2))-0.05 max(data1(:, 2))+0.05], ...
         'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'normal');

xticks = get(gca, 'XTick');
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
yticks = get(gca, 'YTick');
yticklabels = arrayfun(@(y) sprintf('%.2f', y), yticks, 'UniformOutput', false);
set(gca, 'XTickLabel', xticklabels,'YTickLabel', yticklabels); 

subplot(3, 3, 5); 
hold on;
x1 = normrnd(A(2).mu, A(2).sigma, 3000, 1);
y1 = normrnd(A(4).mu, A(4).sigma, 3000, 1);
data1 = [x1, y1];
density_2D_1 = density2D_KD(data1(:, 1:2), 0.01);
scatter(data1(:, 1), data1(:, 2), 15, density_2D_1, 'filled');
xlabel('$a(x^3)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
ylabel('$b(y^3)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
colormap(map);
colorbar;   
set(gca, 'Box', 'off','LineWidth', 1,'XGrid', 'off', 'YGrid', 'off','TickDir', 'out', 'TickLength', [.005 .005],'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', min(data1(:, 1))-0.05:0.25:max(data1(:, 1)+0.05), ...
         'XLim', [min(data1(:, 1))-0.05 max(data1(:, 1))+0.05], ...
         'YTick', min(data1(:, 2))-0.05:0.15:max(data1(:, 2)+0.05), ...
         'YLim', [min(data1(:, 2))-0.05 max(data1(:, 2))+0.05], ...
         'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'normal');

xticks = get(gca, 'XTick');
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
yticks = get(gca, 'YTick');
yticklabels = arrayfun(@(y) sprintf('%.2f', y), yticks, 'UniformOutput', false);
set(gca, 'XTickLabel', xticklabels,'YTickLabel', yticklabels); 

subplot(3, 3, 6); 
hold on;
x1 = normrnd(A(2).mu, A(2).sigma, 3000, 1);
y1 = normrnd(A(5).mu, A(5).sigma, 3000, 1);
data1 = [x1, y1];
density_2D_1 = density2D_KD(data1(:, 1:2), 0.01);
scatter(data1(:, 1), data1(:, 2), 15, density_2D_1, 'filled');
xlabel('$a(x^3)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
ylabel('$a(xy^2)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
colormap(map);
colorbar;   
set(gca, 'Box', 'off','LineWidth', 1,'XGrid', 'off', 'YGrid', 'off','TickDir', 'out', 'TickLength', [.005 .005],'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', min(data1(:, 1))-0.05:0.3:max(data1(:, 1)+0.05), ...
         'XLim', [min(data1(:, 1))-0.05 max(data1(:, 1))+0.05], ...
         'YTick', min(data1(:, 2))-0.05:0.25:max(data1(:, 2)+0.05), ...
         'YLim', [min(data1(:, 2))-0.05 max(data1(:, 2))+0.05], ...
         'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'normal');

xticks = get(gca, 'XTick');
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
yticks = get(gca, 'YTick');
yticklabels = arrayfun(@(y) sprintf('%.2f', y), yticks, 'UniformOutput', false);
set(gca, 'XTickLabel', xticklabels,'YTickLabel', yticklabels); 

subplot(3, 3, 7); 
hold on;
x1 = normrnd(A(3).mu, A(3).sigma, 3000, 1);
y1 = normrnd(A(4).mu, A(4).sigma, 3000, 1);
data1 = [x1, y1];
density_2D_1 = density2D_KD(data1(:, 1:2), 0.01);
scatter(data1(:, 1), data1(:, 2), 15, density_2D_1, 'filled');
xlabel('$\omega(y)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
ylabel('$b(y^3)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
colormap(map);
colorbar;   
set(gca, 'Box', 'off','LineWidth', 1,'XGrid', 'off', 'YGrid', 'off','TickDir', 'out', 'TickLength', [.005 .005],'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', min(data1(:, 1))-0.05:0.4:max(data1(:, 1)+0.05), ...
         'XLim', [min(data1(:, 1))-0.05 max(data1(:, 1))+0.05], ...
         'YTick', min(data1(:, 2))-0.05:0.15:max(data1(:, 2)+0.05), ...
         'YLim', [min(data1(:, 2))-0.05 max(data1(:, 2))+0.05], ...
         'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'normal');

xticks = get(gca, 'XTick');
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
yticks = get(gca, 'YTick');
yticklabels = arrayfun(@(y) sprintf('%.2f', y), yticks, 'UniformOutput', false);
set(gca, 'XTickLabel', xticklabels,'YTickLabel', yticklabels); 

subplot(3, 3, 8); 
hold on;
x1 = normrnd(A(3).mu, A(3).sigma, 3000, 1);
y1 = normrnd(A(5).mu, A(5).sigma, 3000, 1);
data1 = [x1, y1];
density_2D_1 = density2D_KD(data1(:, 1:2), 0.01);
scatter(data1(:, 1), data1(:, 2), 15, density_2D_1, 'filled');

xlabel('$\omega(y)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
ylabel('$a(xy^2)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
colormap(map);
colorbar;   
set(gca, 'Box', 'off','LineWidth', 1,'XGrid', 'off', 'YGrid', 'off','TickDir', 'out', 'TickLength', [.005 .005],'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', min(data1(:, 1))-0.05:0.5:max(data1(:, 1)+0.05), ...
         'XLim', [min(data1(:, 1))-0.05 max(data1(:, 1))+0.05], ...
         'YTick', min(data1(:, 2))-0.05:0.25:max(data1(:, 2)+0.05), ...
         'YLim', [min(data1(:, 2))-0.05 max(data1(:, 2))+0.05], ...
         'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'normal');

xticks = get(gca, 'XTick');
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
yticks = get(gca, 'YTick');
yticklabels = arrayfun(@(y) sprintf('%.2f', y), yticks, 'UniformOutput', false);
set(gca, 'XTickLabel', xticklabels,'YTickLabel', yticklabels); 

subplot(3, 3, 9); 
hold on;
x1 = normrnd(A(4).mu, A(4).sigma, 3000, 1);
y1 = normrnd(A(6).mu, A(6).sigma, 3000, 1);
data1 = [x1, y1];
density_2D_1 = density2D_KD(data1(:, 1:2), 0.01);
scatter(data1(:, 1), data1(:, 2), 15, density_2D_1, 'filled');

xlabel('$b(y^3)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
ylabel('$b(x^2y)$', 'Interpreter', 'latex','FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
colormap(map);
colorbar;   
set(gca, 'Box', 'off','LineWidth', 1,'XGrid', 'off', 'YGrid', 'off','TickDir', 'out', 'TickLength', [.005 .005],'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', min(data1(:, 1))-0.05:0.25:max(data1(:, 1)+0.05), ...
         'XLim', [min(data1(:, 1))-0.05 max(data1(:, 1))+0.05], ...
         'YTick', min(data1(:, 2))-0.05:0.2:max(data1(:, 2)+0.05), ...
         'YLim', [min(data1(:, 2))-0.05 max(data1(:, 2))+0.05], ...
         'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'normal');

xticks = get(gca, 'XTick');
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
yticks = get(gca, 'YTick');
yticklabels = arrayfun(@(y) sprintf('%.2f', y), yticks, 'UniformOutput', false);
set(gca, 'XTickLabel', xticklabels,'YTickLabel', yticklabels); 

set(gcf, 'Color', [1 1 1]);
%% 图片输出
figW = figureWidth;
figH = figureHeight;
set(figureHandle,'PaperUnits',figureUnits);
set(figureHandle,'PaperPosition',[0 0 figW figH]);
fileout = '(二维)漂移项系数联合后验分布';
print(figureHandle,[fileout,'.png'],'-r600','-dpng');


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

