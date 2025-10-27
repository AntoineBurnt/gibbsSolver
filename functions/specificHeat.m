function cp = specificHeat(inputs,T,coeffs)
    R = inputs.R;
    if iscell(coeffs)
        coeffs = coeffs{1};
    end
    n_intervals = length(coeffs(:,1));

for i = 1:n_intervals
        interval = i;
        if (coeffs(i,2) > T)
            break

        end
    end

    a = coeffs(interval,3:9);

    cp = R*(a(1)*T^-2 + a(2)*T^-1 + a(3) + a(4)*T + a(5)*T^2 + a(6)*T^3 + a(7)*T^4);
end