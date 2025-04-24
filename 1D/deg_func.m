function [func_list, str_list] = deg_func(order)
    
    % 生成x^0 到 x^order
    x_funcs = {};
    x_strs = {};
    for i = 0:order
        x_funcs{end+1} = @(x) x.^i;
        x_strs{end+1} = sprintf('x^%d', i);
    end
    %生成exp
    for i = 1:order
        x_funcs{end+1} = @(x) exp(i * x);
        x_strs{end+1} = sprintf('exp(%dx)',i);
    end
    
    % 生成分数 1./(1 + 基函数)
    frac_funcs = {};
    frac_strs = {};
    for i = 2:length(x_funcs)
        frac_funcs{end+1} = @(x) 1 ./ ( 1 + x_funcs{i}(x));
        frac_strs{end+1} = sprintf('1/(1+%s)',x_strs{i});    
    end

    % 合并only_x和mul_x
    func_list = [x_funcs,frac_funcs];
    str_list = [ x_strs, frac_strs];
    
     %生成三角函数
     for i = 1:order
        func_list{end+1} = @(x) sin(i * x);
        str_list{end+1} = sprintf('sin(%dx)',i);
     end

    for i = 1:order
        func_list{end+1} = @(x) cos(i * x);
        str_list{end+1} = sprintf('cos(%dx)',i);
    end

end