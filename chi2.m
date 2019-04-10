function [tbl, chi2stat, pval, x, y] = chi2(n, N)
    
    y = [];
    x = [];
    for i = 1:length(N)
        y = [y repmat(i, 1, N(i))];
        x = [x repmat(1, 1, n(i))];
        x = [x repmat(2, 1, N(i) - n(i))];
    end

    [tbl, chi2stat, pval] = crosstab(x, y);
