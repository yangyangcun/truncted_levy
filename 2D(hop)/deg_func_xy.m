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
        str_list{end+1} = sprintf('x^%d', i);
    end

    % 生成y的函数：y^0 到 y^order
    for j = 1:order
        func_list{end+1} = @(x, y) y.^j;
        str_list{end+1} = sprintf('y^%d', j);
    end

    % 生成x和y的乘积函数：x^i * y^j
    for i = 1:order
        for j = 1:order
            func_list{end+1} = @(x, y) x.^i .* y.^j;
            str_list{end+1} = sprintf('x^%dy^%d', i, j);
        end
    end
end