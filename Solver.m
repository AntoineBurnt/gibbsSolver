clear variables
close all
clc
addpath('functions')
%% Inputs

%Fuel Mixture
fuel = ["C12H24"];
sigma_f = [1]; %[kmol]

%Oxidizer Mixture
oxid = ["N2O"];
sigma_o = [1]; %[kmol]

%Fuel/Oxidizer mixture temperature
T_initial = 298.15; %[K]

%F/O Ratio
%fo = 1/18;
equivalence = 1;

%Combustion Pressure
P = 101325*10; %[Pa]

%Target Temperature
T_target = 3383; %[K]

%Gibbs minimization tolerance
tol = 1*10^-6;

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

%Initial pre-combustion mixture vector
sigma_pre = horzcat(equivalence*fo_st*sigma_f,sigma_o);

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

%% Log setup information
disp('* Atoms present *'  )
disp(join(atoms,', '))

disp('* Species under consideration *')
disp(join([species.name],', '))

disp(['* Pre-combustion mixture temperature: ',num2str(T_initial) ,'K *'])
disp(['* Combustion Pressure: ', num2str(P/P_ref),'atm *'])
disp(['* Total enthaly: ', num2str(H),'[J] *'])
disp('*******************************************************************')

%Call solver
[sigma,sigma_total,T_eq,iteration] = solve(H,B,etta,coeffs,constants);


%% Results
disp('*******************************************************************')

for x = 1:n_s
    active_species(x) = [species(x).name];
    M_active(x) = M(x);
end

sigma_fraction = sigma/sigma_total;


mass_fraction= (sigma.*M_active)/(sum(sigma.*M_active));
mass_fraction = mass_fraction*100;

sig_indexs = [];

for x = 1:n_s
    if mass_fraction(x) > 0.01
        sig_indexs = horzcat(sig_indexs,x);
    end
end

sig_mass_fraction = mass_fraction(sig_indexs);
sig_species = active_species(sig_indexs);


disp(sig_species)
disp(sig_mass_fraction)

%Species Mass Fraction
figure('Name','Species Mass Fraction');
bar(sig_species,sig_mass_fraction);
grid on
ylim([0,100]);
%legend(['\phi = ', num2str(round(equivalence,2))],'Location','northwest');

ylabel("Mass Fraction [%]")
xlabel("Species")
title('Species Mass Fraction');
set(gca, 'fontsize', 14)
set(gcf, 'Position',  [100, 100, 750, 600])
