function [func_list, str_list] = deg_func_xy(order)

    func_list = {};
    str_list = {};

    % 生成x的函数：x^0 到 x^order
    for i = 0:order
        func_list{end+1} = @(x, y, z) x.^i;
        str_list{end+1} = sprintf('x^%d', i);
    end

    % 生成y的函数：y^0 到 y^order
    for j = 1:order
        func_list{end+1} = @(x, y, z) y.^j;
        str_list{end+1} = sprintf('y^%d', j);
    end
    
    for j = 1:order
        func_list{end+1} = @(x, y, z) z.^j;
        str_list{end+1} = sprintf('z^%d', j);
    end


    for i = 1:order
        for j = 1:order
            func_list{end+1} = @(x, y, z) x.^i .* y.^j;
            str_list{end+1} = sprintf('x^%dy^%d', i, j);
        end
    end

    % 生成 x 和 z 的乘积函数：x^i * z^k
    for i = 1:order
        for k = 1:order
            func_list{end+1} = @(x, y, z) x.^i .* z.^k;
            str_list{end+1} = sprintf('x^%dz^%d', i, k);
        end
    end

    % 生成 y 和 z 的乘积函数：y^j * z^k
    for j = 1:order
        for k = 1:order
            func_list{end+1} = @(x, y, z) y.^j .* z.^k;
            str_list{end+1} = sprintf('y^%dz^%d', j, k);
        end
    end
end