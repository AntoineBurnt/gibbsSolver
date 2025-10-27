clear variables
close all
clc
%% Inputs

%Fuel Mixture
fuel = ["H2"];
sigma_f = [1]; %[kmol]

%Oxidizer Mixture
oxid = ["O2"];
sigma_o = [1]; %[kmol]


%Fuel/Oxidizer mixture temperature
T_initial = 298.15; %[K]

%Combustion Pressure
P = 101325*10; %[Pa]

%Target Temperature
T_target = 3383; %[K]

%Gibbs minimization tolerance
tol = 1*10^-7;

%Species data source (0 is textbook species, 1 is from NASA database)
source = 1;

%Fuel/Oxidizer mixtures check
if length(fuel) ~= length(sigma_f) || length(oxid) ~= length(sigma_o)
    error("Incorrect fuel or oxidizer mixture!")
end

%Initialize Constants
T_ref = 298.15; %[K]
R = 8314.46261815324; %[J/(kmol*K)]
P_ref = 100000; %[Pa]

%Identify reactants
reactants = unique(horzcat(fuel,oxid),'stable');
%% Species Data Initialization
[species,atoms] = speciesLoader('data/species.dat',reactants,source);

n_a = length(atoms);
atoms_dict = dictionary(atoms,1:n_a);

n_s = length(species);
species_dict = dictionary([species.name],1:n_s);

%Construct etta matrix
etta = zeros(n_a,n_s);

for i = 1:n_s
    spec = species(i);

    for ii = 1:length(spec.atoms)
        etta(atoms_dict(spec.atoms(ii)),i) = spec.etta(ii);
    end
end

%Molecular Weight Vector
M = [species.weight]; %[kg/kmol]


%Construct coefficients matrix
coeffs = cell(1,n_s);

for i = 1:n_s
    coeffs(i) = {horzcat(species(i).intervals,species(i).coefficients)};
end

% Each row of the matrix represents a temperature interval. The first two
% columns represent the temperature interval over which the coefficients
% are valid. The remaining columns represents coefficients used for
% determining the thermodynamic properties. Each row follow this format:
% [lower higher a_1 a_2 a_3 a_4 a_5 a_6 a_7 b_1 b_2] 

%Initialize constants struct
constants = struct();
constants.T_ref = T_ref;
constants.P_ref = P_ref;
constants.R = R;
constants.P = P;
constants.n_a = n_a;
constants.n_s = n_s;
constants.tol = tol;

%Calculate stoiciometric F/O ratio
fo_st = STbalance(etta,atoms_dict,species_dict,fuel,sigma_f,oxid,sigma_o);
fo = fo_st/10;

fo_prev = 0;
T_eq_prev = T_initial;


%% Log setup information
disp('* Atoms present *'  )
disp(join(atoms,', '))

disp('* Species under consideration *')
disp(join([species.name],', '))

disp(['* Pre-combustion mixture temperature: ',num2str(T_initial) ,'K *'])
disp(['* Combustion Pressure: ', num2str(P/P_ref),'atm *'])
disp(['* Target post-combustion temperature: ',num2str(T_target) ,'K *'])
disp('*******************************************************************')

T_eqs = [];
equivalences = [];

phi_s = [0 0];
sigma_s = zeros(2,n_s);

equivalence = 0;
control = 0;

while control < 3
    %Mixture mass
    m = 0;

    for i = 1:length(fuel)
        m = m + sigma_f(i)*fo*M(species_dict(fuel(i)))
    end

    for i = 1:length(oxid)
        m = m + sigma_o(i)*M(species_dict(oxid(i)))
    end

 sigma_pre = horzcat(fo*sigma_f,sigma_o)/m;

%Pre-combustion conditions
H = 0; %[J]
B = zeros(1,n_a); %[kmol]

reactant_ids = species_dict(reactants);

for x = 1:length(reactant_ids)
    %Calculate B vector
    for i = 1:n_a
        B(i) = B(i) + etta(i,reactant_ids(x))*sigma_pre(x);
    end

    %Calculate enthalpy
    H = H + enthalpy(constants,T_initial,coeffs{reactant_ids(x)})*sigma_pre(x);
end

%Call solver
[sigma,sigma_total,T_eq,iteration] = solve(H,B,etta,coeffs,constants);

equivalence = fo/fo_st;

%Flow control logic
if control == 0
    if equivalence >= 3.1
        control = 1;
        fo = fo_st;
        if fo == fo_st
            if max(T_eqs) < T_target
                error("T_eq at fo_st was lower than T_target, lower T_target and try again")
            end
        end
    else
        T_eqs = horzcat(T_eqs,T_eq);
        equivalences = horzcat(equivalences,equivalence);
        fo = fo + 0.1*fo_st;
    end
else
    if (abs(T_target - T_eq) < tol)
        
        phi_s(control) = equivalence;
        sigma_s(control,:) = sigma;

        control = control + 1;
        fo = fo*5;
    else
        dfodT = (fo_prev-fo)/(T_eq_prev-T_eq);

        fo_prev = fo;
        T_eq_prev = T_eq;

        fo = max(fo + clip(dfodT*(T_target - T_eq),-fo*0.25,fo*0.25),0);

    end
end

end
%% Post-process

%Molar mass of products
M_prod = sum(sigma.*M)/sum(sigma);

%Density of products
rho_prod = P/(sigma_total*R*T_eq);

%% Results
disp('*******************************************************************')
disp('Equivalence ratios')
disp(phi_s)

for x = 1:n_s
    active_species(x) = species(x).name;
end

sigma_fraction = sigma/sigma_total;
mass_fraction = zeros(2,n_s);

mass_fraction(1,:) = (sigma_s(1,:).*M)/(sum(sigma_s(1,:).*M));
mass_fraction(2,:) = (sigma_s(2,:).*M)/(sum(sigma_s(2,:).*M));

mass_fraction = mass_fraction*100;

sig_indexs = [];

for x = 1:n_s
    if mass_fraction(1,x) > 0.01 ||  mass_fraction(2,x) > 0.01
        sig_indexs = horzcat(sig_indexs,x);
    end
end

sig_mass_fraction = mass_fraction(:,sig_indexs);
sig_species = active_species(sig_indexs);


disp(sig_species)
disp(sig_mass_fraction)

%Species Mass Fraction
figure('Name','Species Mass Fraction');
bar(sig_species,sig_mass_fraction);
grid on
ylim([0,100]);
legend(['\phi_1 = ', num2str(round(phi_s(1),2))], ...
    ['\phi_2 = ', num2str(round(phi_s(2),2))],'Location','northwest');

ylabel("Mass Fraction [%]")
xlabel("Species")
%title('Species Mass Fraction');
set(gca, 'fontsize', 14)
set(gcf, 'Position',  [0, 100, 700, 600])

%t_EQ VS PHI
figure('Name','Equilibrium Temperature');
plot(equivalences,T_eqs,"LineWidth",2)
hold on
plot(phi_s,T_target,"Marker","x","MarkerSize",10);
grid on
ylim([0,4000]);
xlim([0 max(equivalences)+0.1])
ylabel("Equilibrium Temperature, T_{eq} [K]")
xlabel("Equivalence Ratio, \phi")
%title('CHange this');
set(gca, 'fontsize', 14)
set(gcf, 'Position',  [700, 100, 700, 600])