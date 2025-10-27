  function fo_st = STbalance(etta,active_a_dict,active_s_dict,fuel,sigma_f,oxid,sigma_o)
%% Stoichiometric F/O ratio
%This section determines the stoichiometric fuel to oxidizer ratio of the
%fuel and oxidizer mixtures defined in the Input section. 
%a(Fuel) + b(Oxidizer) = cH2O + dN2 + eCO2

%Number of atoms per "unit" of fuel or oxidizer mixture
combustion_atoms = ["H" "O" "N" "C"];
oxidizer_atoms = zeros(1,4);
fuel_atoms = zeros(1,4);

for i = 1:4
    if isKey(active_a_dict,combustion_atoms(i))
        for ii = 1:length(fuel)
            fuel_atoms(i) = sum(etta(active_a_dict(combustion_atoms(i)),active_s_dict(fuel)).*sigma_f);
            oxidizer_atoms(i) = sum(etta(active_a_dict(combustion_atoms(i)),active_s_dict(oxid)).*sigma_o,2);
        end
    else
        fuel_atoms(i) = 0;
        oxidizer_atoms(i) = 0;
    end
end

%Combustion reaction balancing
%aFuel + bOxidizer = cH2O + dN2 + eCO2
balance = [fuel_atoms(1)  oxidizer_atoms(1) -2 0 0;
           fuel_atoms(2) oxidizer_atoms(2) -1 0 -2;
           fuel_atoms(3) oxidizer_atoms(3) 0 -2 0;
           fuel_atoms(4) oxidizer_atoms(4) 0 0 -1];

[balance,pivots] = rref(balance);
col = setdiff(1:5,pivots);

%Stoichiometric F/O ratio
fo_st = balance(1,col)/balance(2,col);
end