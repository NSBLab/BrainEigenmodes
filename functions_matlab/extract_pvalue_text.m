function pvalue_text = extract_pvalue_text(pval, is_corrected, correction_subscript)

if nargin<3
    correction_subscript = 'FWER';
end
if nargin<2
    is_corrected = 0;
end
    
if ~is_corrected
    if pval>=0.01
        pvalue_text = sprintf('p = %.3f', pval);
    elseif pval<0.01 && pval>=0.001
        pvalue_text = 'p < 0.01';
    elseif pval<0.001
        pvalue_text = 'p < 0.001';
    end
else
    if pval>=0.01
        pvalue_text = sprintf('p_{%s} = %.3f', correction_subscript, pval);
    elseif pval<0.01 && pval>=0.001
        pvalue_text = sprintf('p_{%s} < 0.01', correction_subscript);
    elseif pval<0.001
        pvalue_text = sprintf('p_{%s} < 0.001', correction_subscript);
    end
end