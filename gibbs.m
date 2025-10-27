function [g,h,cp,s] = gibbs(inputs,T,coeffs,P,sigma_k,sigma_total)
    R = inputs.R;
    P_ref = inputs.P_ref;
    if iscell(coeffs)
        coeffs = coeffs{1};
    end
    if sigma_k == 0
        g = 0;
        h = 0;
        cp = 0;
        s = 0;
    else
        h = enthalpy(inputs,T,coeffs);
        s = entropy(inputs,T,coeffs);
        cp = specificHeat(inputs,T,coeffs);
        g = h - T*s + R*T*log(sigma_k/sigma_total) + R*T*log(P/P_ref);
    end
end