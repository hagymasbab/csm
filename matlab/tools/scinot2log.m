function l = scinot2log(coeff,expo)
    s = sign(coeff);
    if s == -1
        error('logarithm of negative number');
    end
    act_log10 = expo + log10(abs(coeff));
    l = s * act_log10 / log10(exp(1));
end
