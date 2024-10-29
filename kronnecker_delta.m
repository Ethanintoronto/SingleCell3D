function newExp = kronnecker_delta(exp, new)
old = {'i', 'j'};
    for v = symvar(exp)
        new_var = v; 
        for i = 1:numel(old)
            new_var = regexprep(char(new_var), old(i), new(i)); 
        end
        exp = subs(exp, v, sym(new_var));     % Substitute old with new
    end
    newExp = exp;
end 