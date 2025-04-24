clc
clear
close all


rng(100, 'twister');

epsilon = 0.05;
num_N = 25000;


alpha1 = 0.5;
c_alpha1 = alpha1 * gamma((1+alpha1)/2) / (2^(1-alpha1) * pi^(0.5) * gamma(1-alpha1/2));
alpha2 = 1;
c_alpha2 = alpha2 * gamma((1+alpha2)/2) / (2^(1-alpha2) * pi^(0.5) * gamma(1-alpha2/2));
alpha3 = 1.5;
c_alpha3 = alpha3 * gamma((1+alpha3)/2) / (2^(1-alpha3) * pi^(0.5) * gamma(1-alpha3/2));

Lt1 = generate_Lt(0.5,epsilon,num_N,0.1,c_alpha1);
Lt2 = generate_Lt(1.5,epsilon,num_N,0.01,c_alpha3);

L1 =  Lt1(2:end) - Lt1(1:end-1);
L2 =  Lt2(2:end) - Lt2(1:end-1);


Lt3 = generate_Lt(0.5,epsilon,num_N,0.1,c_alpha1);
Lt4 = generate_Lt(0.5,epsilon,num_N,0.01,c_alpha1);
Lt5 = generate_Lt(1,epsilon,num_N,0.1,c_alpha2);
Lt6 = generate_Lt(1,epsilon,num_N,0.01,c_alpha2);
Lt7 = generate_Lt(1.5,epsilon,num_N,0.1,c_alpha3);
Lt8 = generate_Lt(1.5,epsilon,num_N,0.01,c_alpha3);


L3 =  Lt3(2:end) - Lt3(1:end-1);
L4 =  Lt4(2:end) - Lt4(1:end-1);
L5 =  Lt5(2:end) - Lt5(1:end-1);
L6 =  Lt6(2:end) - Lt6(1:end-1);
L7 =  Lt7(2:end) - Lt7(1:end-1);
L8 =  Lt8(2:end) - Lt8(1:end-1);

%%
figureUnits = 'centimeters';
figureWidth = 24; % 宽度增加以适应两个子图
figureHeight = 20;
close all

figureHandle = figure;
set(gcf, 'Units', figureUnits, 'Position', [2 5 figureWidth figureHeight]);

% 子图布局参数
topMargin = 0.01;       % 上边距
bottomMargin = 0.05;    % 下边距（为标题留空间）
leftMargin = 0.09;
rightMargin = 0.02;
subplotGap = 0.1;      % 子图间距

% 定义子图的位置
positions = [
    0.09, 0.65, 0.37, 0.34; % 第一张图
    0.54, 0.65, 0.37, 0.34; % 第三张图
    0.09, 0.18, 0.37, 0.34; % 第四张图（居中，第一张和第二张之间）
    0.54, 0.18, 0.37, 0.34; % 第五张图（居中，第二张和第三张之间）
];

subplot('Position', positions(1, :)); 
hold on;
T = 0.1:0.1:0.1 * num_N;
plot(T, L1);
hXLabel1 = xlabel('T');
hYLabel1 = ylabel('L(t)');
set(gca, 'Box', 'off', ...
         'LineWidth', 1, ...
         'XGrid', 'off', 'YGrid', 'off', ...
         'TickDir', 'out', 'TickLength', [.005 .005], ...
         'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', 0:1000:2500, ...
         'XLim', [0 2500], ...
         'YTick', -2:1:2, ...
         'YLim', [-2 2]);
legend('$\alpha=0.5,\triangle t =0.1$', 'Interpreter', 'latex', ...
    'Location', 'northeast', ... % 图例位置在右上角
    'FontSize', 13, ...
    'Box', 'off','FontWeight' , 'normal','FontName', 'Times New Roman'); 
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15,'FontWeight' , 'normal')
set([hXLabel1, hYLabel1], 'FontSize', 15, 'FontName', 'Times New Roman','FontWeight' , 'normal')
pos = get(gca, 'Position');

subplot('Position', positions(3, :)); 
hold on;
T1 = 0.01:0.01:0.01 * num_N;
plot(T1, L2);
hXLabel1 = xlabel('T');
hYLabel1 = ylabel('L(t)');
set(gca, 'Box', 'off', ...
         'LineWidth', 1, ...
         'XGrid', 'off', 'YGrid', 'off', ...
         'TickDir', 'out', 'TickLength', [.005 .005], ...
         'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'XTick', 0:150:250, ...
         'XLim', [0 250], ...
         'YTick', -1.5:1:1.5, ...
         'YLim', [-1.5 1.5]);
legend('$\alpha=1.5,\triangle t =0.01$', 'Interpreter', 'latex', ...
    'Location', 'northeast', ... % 图例位置在右上角
    'FontSize', 13, ...
    'Box', 'off','FontWeight' , 'normal','FontName', 'Times New Roman'); 
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15,'FontWeight' , 'normal')
set([hXLabel1, hYLabel1], 'FontSize', 15, 'FontName', 'Times New Roman','FontWeight' , 'normal')


C1 = addcolorplus(119); 
C2 = addcolorplus(252); 
C3 = addcolorplus(238); 


[pdf_L3, x_L3] = ksdensity(L3, 'Bandwidth', 0.1, 'NumPoints', 2000);
[pdf_L5, x_L5] = ksdensity(L5, 'Bandwidth', 0.1, 'NumPoints', 2000);
[pdf_L7, x_L7] = ksdensity(L7, 'Bandwidth', 0.1, 'NumPoints', 2000);

subplot('Position', positions(2, :)); 
hold on;
plot(x_L3, pdf_L3, 'Color', C1, 'LineWidth', 1.5);
plot(x_L5, pdf_L5, 'Color', C2, 'LineWidth', 1.5);
plot(x_L7, pdf_L7, 'Color', C3, 'LineWidth', 1.5);
set(gca, 'Box', 'off', 'XGrid', 'off', 'YGrid', 'off', ...
         'TickDir', 'out', 'TickLength', [.01 .01], ...
         'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'normal', ...
         'XTick', -1:0.5:1, ...
         'XLim', [-1 1], ...
         'YTick', 0:2.5:5, ...
         'YLim', [0 5]);
legend({'$\alpha=0.5,\triangle t =0.1$', '$\alpha=1,\triangle t =0.1$','$\alpha=1.5,\triangle t =0.1$'},  'Interpreter', 'latex',...
    'Location', 'northeast', ... % 图例位置在右上角
    'FontSize', 13, ...
    'Box', 'off','FontWeight' , 'normal','FontName', 'Times New Roman'); % 去掉图例边框
hXLabel1 = xlabel('L(t)');
hYLabel1 = ylabel('PDF');
set([hXLabel1, hYLabel1], 'FontSize', 15, 'FontName', 'Times New Roman','FontWeight' , 'normal')


[pdf_L4, x_L4] = ksdensity(L4, 'Bandwidth', 0.035, 'NumPoints', 2000);
[pdf_L6, x_L6] = ksdensity(L6, 'Bandwidth', 0.035, 'NumPoints', 2000);
[pdf_L8, x_L8] = ksdensity(L8, 'Bandwidth', 0.035, 'NumPoints', 2000);
subplot('Position', positions(4, :)); 
hold on;
plot(x_L4, pdf_L4, 'Color', C1, 'LineWidth', 1.5);
plot(x_L6, pdf_L6, 'Color', C2, 'LineWidth', 1.5);
plot(x_L8, pdf_L8, 'Color', C3, 'LineWidth', 1.5);

set(gca, 'Box', 'off', 'XGrid', 'off', 'YGrid', 'off', ...
         'TickDir', 'out', 'TickLength', [.01 .01], ...
         'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
         'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'normal', ...
         'XTick', -0.5:0.25:0.5, ...
         'XLim', [-0.5 0.5], ...
         'YTick', 0:7:16, ...
         'YLim', [0 16]);
legend({'$\alpha=0.5,\triangle t =0.01$', '$\alpha=1,\triangle t =0.01$','$\alpha=1.5,\triangle t =0.01$'},  'Interpreter', 'latex',...
    'Location', 'northeast', ... % 图例位置在右上角
    'FontSize', 13, ...
    'Box', 'off','FontWeight' , 'normal','FontName', 'Times New Roman'); % 去掉图例边框
hXLabel1 = xlabel('L(t)');
hYLabel1 = ylabel('PDF');
set([hXLabel1, hYLabel1], 'FontSize', 15, 'FontName', 'Times New Roman','FontWeight' , 'normal')

set(gcf, 'Color', [1 1 1]); % 设置背景为白色


%% 图片输出
figW = figureWidth;
figH = figureHeight;
set(figureHandle,'PaperUnits',figureUnits);
set(figureHandle,'PaperPosition',[0 0 figW figH]);
fileout = 'levy随机数生成';
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


