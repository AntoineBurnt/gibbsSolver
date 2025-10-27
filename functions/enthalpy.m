function h = enthalpy(inputs,T,coeffs)
    R = inputs.R;

    if iscell(coeffs)
        coeffs = coeffs{1};
    end
    n_intervals = length(coeffs(:,1));

    for i = 1:n_intervals
        interval = i;
        if (coeffs(interval,2) > T)
            break

        end
    end

    a = coeffs(interval,3:10);

    h = (R*T*(-a(1)*T^-2 + a(2)*log(T)/T + a(3) + a(4)*T/2 + a(5)*T^2/3 + ...
        a(6)*T^3/4 + a(7)*T^4/5 + a(8)/T));

end