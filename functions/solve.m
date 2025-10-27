function [sigma,sigma_total,T_eq,iteration] = solve(H,B,etta,coeffs,constants)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

tol = constants.tol;
n_s = constants.n_s;
n_a = constants.n_a;
R = constants.R;

constants_vect = repmat(constants,1,n_s);


%Iteration Initial Conditions
sigma = (1/n_s)*ones(1,n_s);
sigma_total = sum(sigma);
T_eq = 4000;

converged = 0;
iteration = 0;

%% Inner Loop Start
while ~converged
    %Determine h, c_p, g of each species
    T_eq_vect = T_eq*ones(1,n_s);
    P_vect = constants.P*ones(1,n_s);
    sigma_total_vect = sigma_total*ones(1,n_s);


    [g_vect,h_vect,cp_vect,~] = arrayfun(@gibbs,constants_vect,T_eq_vect,coeffs,P_vect,sigma,sigma_total_vect);

    %% A Matrix
    A = zeros(n_a+2,n_a+2);

    for i = 1:n_a
        for ii = 1:n_a
            A(i,ii) = sum(etta(i,:).*sigma.*etta(ii,:));
        end
        A(i,n_a+1) = sum(etta(i,:).*sigma);
        A(i,n_a+2) = sum(etta(i,:).*sigma.*(h_vect./(R*T_eq)));
    end

    for ii = 1:n_a
        A(n_a + 1,ii) = sum(etta(ii,:).*sigma.*(h_vect./(R*T_eq)));
        A(n_a + 2,ii) = sum(etta(ii,:).*sigma);
    end

    A(n_a+1,n_a+1) = sum(sigma.*h_vect./(R*T_eq));
    A(n_a+1,n_a+2) = sum(sigma.*((cp_vect./R)+(h_vect./(R*T_eq)).^2));
    A(n_a+2,n_a+1) = sum(sigma) - (sigma_total);
    A(n_a+2,n_a+2) = A(n_a+1,n_a+1);

    %% C vector
    C = zeros(n_a+2,1);

    for i = 1:n_a
        C(i) = B(i) + sum(etta(i,:).*sigma.*((g_vect./(R*T_eq))-1));
    end

    C(n_a+1) = (H/(R*T_eq)) + sum(sigma.*(h_vect/(R*T_eq)).*((g_vect./(R*T_eq))-1));
    C(n_a+2) = sigma_total + sum(sigma.*((g_vect./(R*T_eq))-1));

    %Compute solution vector
    X = A^-1*C;

    %Calculate delta_sigma vector

    delta_ln_sigma = zeros(1,n_s);
    for i = 1:n_s
        delta_ln_sigma(i) = X(n_a+1) + ((h_vect(i)/(R*T_eq))*X(n_a+2)) - (g_vect(i)/(R*T_eq));
        for ii = 1:n_a
            delta_ln_sigma(i) = delta_ln_sigma(i) + X(ii)*etta(ii,i);
        end
    end

    %Under relaxation parameters

    beta_1 = ones(1,n_s);
    beta_2 = ones(1,n_s);

    for x = 1:n_s
        if (sigma(x) > 0.00000001)
            beta_1(x) = 2/(max([abs(X(n_a+2)) abs(X(n_a+1)) abs(delta_ln_sigma(x))]));
        elseif (delta_ln_sigma(x) > 0)
            beta_2(x) = abs((log(0.0001)-log(sigma(x))+log(sigma_total))./(delta_ln_sigma(x)-X(n_a+1)));
        end
    end
    beta = min([1 beta_1 beta_2]);


    %Compute values for next iteration
    T_new =  exp(log(T_eq) + beta*X(n_a+2));
    sigma_total_new = exp(log(sigma_total) + beta*X(n_a+1));
    sigma_new = exp(log(sigma) + beta*delta_ln_sigma);

    T_eq = T_new;
    sigma_total = sigma_total_new;
    sigma = sigma_new;

    %Convergence Tests
    sigma_test = (sum(sigma.*abs(delta_ln_sigma))/sum(sigma) <= 0.5*10^-5);
    sigma_total_test = (sigma_total*X(n_a+1)/sum(sigma) <= 0.5*10^-5);
    T_test = X(n_a+2) <= 1*10^-4;
    Element_test = all((abs(transpose(B) - sum(etta.*sigma,2)) <= max(B) * 1*10^-6));
    
    if (sigma_test && sigma_total_test && T_test && Element_test)
        converged = 1;
    end

    iteration = iteration + 1;
    if iteration > 100
        error("Failed to converge");
    end
end
%% Thermodynamic Derivatives

% Logarithm of temperature with constant Pressure
A_t = zeros(n_a + 1,n_a + 1);

for i = 1:n_a
    for ii = 1:n_a
        A_t(i,ii) = sum(etta(i,:).*sigma.*etta(ii,:));
    end
end

for i = 1:n_a
    A_t(n_a+1,i) = sum(etta(i,:).*sigma);
end

for i = 1:n_a
    A_t(i,n_a+1) = sum(etta(i,:).*sigma);
end

C_t = zeros(n_a + 1,1);

for i = 1:n_a
    C_t(i) = -sum(etta(i,:).*sigma.*h_vect/(R*T_eq));
end 

C_t(n_a+1) = -sum(sigma.*h_vect/(R*T_eq));


X_t = A_t^(-1)*C_t;

dalphasdlnT = X_t(1:n_a);

dlnsigma_totaldlnT = X_t(n_a+1);

for i = 1:n_s
    dlnsigmasdlnT(i) = h_vect(i)/(R*T_eq) + dlnsigma_totaldlnT + sum(etta(:,i).*dalphasdlnT);
end

dlnVdlnT = 1 + dlnsigma_totaldlnT;


%Logarithm of pressure with constant temperature
A_p = zeros(n_a + 1,n_a + 1);

for i = 1:n_a
    for ii = 1:n_a
        A_p(i,ii) = sum(etta(i,:).*sigma.*etta(ii,:));
    end
end

for i = 1:n_a
    A_p(n_a+1,i) = sum(etta(i,:).*sigma);
end

for i = 1:n_a
    A_p(i,n_a+1) = sum(etta(i,:).*sigma);
end

C_p = zeros(n_a + 1,1);

for i = 1:n_a
    C_p(i) = sum(etta(i,:).*sigma);
end 

C_p(n_a+1) = sum(sigma);

X_p = A_p^(-1)*C_p;

dalphasdlnP = X_p(1:n_a);

dlnsigma_totaldlnP = X_p(n_a+1);

for i = 1:n_s
    dlnsigmasdlnP(i) = -1 + dlnsigma_totaldlnP + sum(etta(:,i).*dalphasdlnP);
end

dlnVdlnP = -1 + dlnsigma_totaldlnP;


%Specific heat @ constant Pressure
c_pf = sum(cp_vect.*sigma);

c_pr = sum(sigma.*h_vect.*dlnsigmasdlnT/T_eq);
%Specific heat @ constant Pressure
c_pe = c_pf + c_pr;

%Specific heat @ constant volume
c_v = c_pe + (sigma_total*R*dlnVdlnT^2)/(dlnVdlnP);

gamma = c_pe/c_v;
gamma_s = -gamma/dlnVdlnP;

%Speed of sound


a = sqrt(sigma_total*R*T_eq*gamma_s);

rho = constants.P/(sigma_total*R*T_eq);

U = H - sigma_total*R*T_eq;

G = sum(g_vect.*sigma);
S  = (H - G)/T_eq;
end

