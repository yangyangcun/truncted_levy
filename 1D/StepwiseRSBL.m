function [selected_terms, final_model,final_model_var] = StepwiseRSBL(Phi, Y, lambda, Learn_Lambda, p)
    % 初始调用RSBL算法
    [X, ~, gamma_ind, gamma_est, ~, ~] = RSBL(Phi, Y, lambda, Learn_Lambda, p);
    
    % 获取非零gamma的索引和对应的gamma值
    non_zero_indices = find(abs(gamma_est) > 0);
    non_zero_gammas = gamma_est(non_zero_indices);
    
    % 按gamma值排序，选择贡献最大的函数项
    [~, sorted_indices] = sort(non_zero_gammas, 'descend');
    sorted_indices = non_zero_indices(sorted_indices);
    
    % 初始化选择的函数项和模型
    selected_terms = [];
    final_model = zeros(size(X));
    final_model_var = zeros(size(X));
    best_bic = inf;
    
    % 逐步添加函数项
    for i = 1:length(sorted_indices)
        % 添加当前函数项
        current_term = sorted_indices(i);
        selected_terms = [selected_terms; current_term];
        
        % 使用选择的函数项重新拟合模型
        Phi_selected = Phi(:, selected_terms);
        [X_new, X_var, ~, ~, ~, ~] = RSBL(Phi_selected, Y, lambda, Learn_Lambda, p);
        
        % 计算BIC
        residual = Y - Phi_selected * X_new;
        n = size(Y, 1);
        k = length(selected_terms);
        rss = sum(residual(:).^2); % 残差平方和
        bic = n * log(rss / n) + k * log(n); % BIC公式
        
        % 如果BIC没有显著改善，则停止添加
        if bic >= best_bic
            selected_terms = selected_terms(1:end-1); % 移除最后添加的项
            break;
        else
            best_bic = bic;
            final_model(selected_terms, :) = X_new;
            final_model_var(selected_terms, :) = X_var;
        end
    end
end