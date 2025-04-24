function [func_list, str_list] = deg_func_xy(order)
    % 生成关于x和y的函数库
    % 输入:
    %   order: 阶数
    % 输出:
    %   func_list: 函数句柄列表
    %   str_list: 对应的函数表达式字符串列表

    % 初始化
    func_list = {};
    str_list = {};

    % 生成x的函数：x^0 到 x^order
    for i = 0:order
        func_list{end+1} = @(x, y) x.^i;
        str_list{end+1} = sprintf('x.^%d', i);
    end

    % 生成y的函数：y^0 到 y^order
    for j = 1:order
        func_list{end+1} = @(x, y) y.^j;
        str_list{end+1} = sprintf('y.^%d', j);
    end

    % 生成exp函数：exp(i*x) 和 exp(j*y)
    %for i = 1:order
        %func_list{end+1} = @(x, y) exp(i * x);
        %str_list{end+1} = sprintf('exp(%d * x)', i);
    %end
    %for j = 1:order
        %func_list{end+1} = @(x, y) exp(j * y);
        %str_list{end+1} = sprintf('exp(%d * y)', j);
    %end

    % 生成分数函数：1/(1 + x^i) 和 1/(1 + y^j)
    %for i = 1:order
        %func_list{end+1} = @(x, y) 1 ./ (1 + x.^i);
        %str_list{end+1} = sprintf('1./(1 + x.^%d)', i);
    %end
    %for j = 1:order
        %func_list{end+1} = @(x, y) 1 ./ (1 + y.^j);
        %str_list{end+1} = sprintf('1./(1 + y.^%d)', j);
    %end

    % 生成x和y的乘积函数：x^i * y^j
    for i = 1:order
        for j = 1:order
            func_list{end+1} = @(x, y) x.^i .* y.^j;
            str_list{end+1} = sprintf('x.^%d .* y.^%d', i, j);
        end
    end
end