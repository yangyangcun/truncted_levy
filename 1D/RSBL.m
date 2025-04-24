function [X, X_VAR,gamma_ind, gamma_est, count, gamm] = RSBL(Phi, Y, lambda, Learn_Lambda, p )
% 鲁棒改进版MSBL：通过动态样本加权抑制异常值影响
% 主要改进：在每次迭代中根据残差计算样本权重，加权Phi和Y
%
% 输入参数：
%   Phi: N x M 字典矩阵
%   Y: N x L 观测矩阵
%   lambda: 正则化参数
%   Learn_Lambda: 是否学习lambda（1=是，0=否）
%   varargin: 可选参数（如剪枝阈值、最大迭代次数等）
%
% 输出参数：
%   X: M x L 估计的解矩阵
%   gamma_ind: 非零gamma的索引
%   gamma_est: 最终的gamma值
%   count: 迭代次数
%   gamm: 第一次迭代的gamma值（调试用）

% ==== 参数解析 ====
[N, M] = size(Phi);
[N, L] = size(Y);

% 默认参数
PRUNE_GAMMA = 1e-4;  % 剪枝阈值
EPSILON = 1e-8;      % 停止条件
MAX_ITERS = 700;     % 最大迭代次数
PRINT = 0;           % 是否显示进度
HUBER_DELTA_FACTOR = 1.345;  % Huber权重阈值系数

% ==== 初始化 ====
gamma = 100 * rand(M, 1);  % 初始化gamma
keep_list = (1:M)';        % 初始保留列表
m = length(keep_list);
mu = zeros(M, L);          % 初始化解矩阵
count = 0;                 % 迭代计数器

weights = ones(N, 1);      % 初始化样本权重

% ==== 主循环 ====
while (1)
    % ==== 剪枝：移除gamma过小的项 ====
    if (min(gamma) < PRUNE_GAMMA)
        index = find(gamma > PRUNE_GAMMA);
        gamma = gamma(index);
        Phi = Phi(:, index);
        keep_list = keep_list(index);
        m = length(gamma);
    end
   
    % 保存第一次迭代的gamma值（调试用）
    if count == 1
        gamm = gamma;
    end
    

    % ==== 动态样本加权 ====
    if count > 1  % 首次迭代后开始计算残差
        residual = Y_weighted - Phi_weighted * mu;  % 残差矩阵（N x L）
        residual_norm = sqrt(sum(residual.^2, 2));  % 每个样本的残差范数（N x 1）
        
        % 计算自适应阈值
        mad = median(abs(residual_norm - median(residual_norm))) / 0.6745;
        delta = HUBER_DELTA_FACTOR * mad;
        %p = 1;%F
        %p = 0.01;%G
        weights = ones(N, 1);
        idx = residual_norm > delta;
        weights(idx) = delta ./ (residual_norm(idx)).^p;
   
    end

    % ==== 构建加权后的Phi和Y ====
    W_sqrt = sqrt(weights);  % 权重对角矩阵

    Phi_weighted = W_sqrt .* Phi;   % 加权Phi
    Y_weighted = W_sqrt .* Y;       % 加权Y

    % ==== 更新解矩阵mu ====
    mu_old = mu;
    Gamma = diag(gamma);
    G = diag(sqrt(gamma));

    [U, S, V] = svd(Phi_weighted * G, 'econ');  % 对加权Phi进行SVD

    % 计算Xi和mu
    if size(S, 1) > 1
        diag_S = diag(S);
    else
        diag_S = S(1);
    end
    Xi = G * V * diag((diag_S ./ (diag_S.^2 + lambda + 1e-16))) * U';
    mu = Xi * Y_weighted;  % 使用加权观测值更新mu

    Xi_var = diag((eye(m) - Xi * Phi_weighted) * Gamma);

    % ==== 更新gamma ====
    gamma_old = gamma;
    mu2_bar = sum(abs(mu).^2, 2) / L;
    
    % 修正点：Sigma_w_diag 的计算

    Sigma_w_diag = real(gamma - (sum(Xi' .* (Phi_weighted * Gamma)))');
    gamma = mu2_bar + Sigma_w_diag;
    

    
    % ==== 更新lambda（如果启用） ====
    if Learn_Lambda == 1
        lambda = (norm(Y_weighted - Phi_weighted * mu, 'fro')^2 / L) / (N - m + sum(Sigma_w_diag ./ gamma_old));
    end

    % ==== 检查停止条件 ====
    count = count + 1;
    if (PRINT)
        fprintf('iters: %d, num coeffs: %d, gamma change: %.4f\n', count, m, max(abs(gamma - gamma_old)));
    end
    if (count >= MAX_ITERS)
        break;
    end
    if (size(mu) == size(mu_old))
        dmu = max(max(abs(mu_old - mu)));
        if (dmu < EPSILON)  break;  end;
    end;
end

% ==== 整理输出结果 ====
gamma_ind = sort(keep_list);  % 非零gamma的索引
gamma_est = zeros(M, 1);
gamma_est(keep_list, 1) = gamma;  % 最终的gamma值
X = zeros(M, L);
X(keep_list, :) = mu;  % 最终的解矩阵

X_VAR = zeros(M,L);
X_VAR(keep_list,:) = Xi_var; 
if (PRINT)
    fprintf('\nFinish running Robust MSBL ...\n');
end
return;