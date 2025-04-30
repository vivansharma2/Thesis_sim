% This function computes the utility for the case without DTC depending on the vector x of input taxes, the initial values for the productivity of the clean and dirty sector and the inital environmental quality.
function U = mysimopttaxnew2noDTC(x, Ac0, Ad0, Aa0, S0)
global rho sigma psi phi alpha gamma eta_d eta_c eta_a qsi epsilon delta numsim psi1 psi2 
%%% Setting vectors' sizes
A_c = zeros(numsim,1);
A_d = zeros(numsim,1);
A_a = zeros(numsim,1);
icl = zeros(numsim,1); % indicator for clean being more productive
idr = zeros(numsim,1); % indicator for dirty being more productive

Yc = zeros(numsim,1);
Yd = zeros(numsim,1);
Ya = zeros(numsim,1);
numerc = zeros(numsim,1); % numerator, Yc
numerd = zeros(numsim,1); % numerator, Yd 
numera = zeros(numsim,1); % numerator, Ya 
denomc = zeros(numsim,1); % denominator, Yc
denomd = zeros(numsim,1); % denominator, Yd 
denoma = zeros(numsim,1); % denominator, Ya

xcit = zeros(numsim,1); % usage of machines (for calculating C)
xdit = zeros(numsim,1);
xait = zeros(numsim,1);

Y = zeros(numsim,1);
C = zeros(numsim,1);
S = zeros(numsim,1); 
s_c = eta_c/(eta_c+eta_d+eta_a)*ones(numsim,1); % scientist allocations
s_d = eta_d/(eta_c+eta_d+eta_a)*ones(numsim,1);
s_a = (1 - s_c - s_d).*ones(numsim,1);
tau_d = x(1:numsim);                    % first half = dirty input tax
tau_a = x(numsim + 1 : 2 * numsim);     % second half = adaptive input tax 
tempinc = zeros(numsim,1);
Omega = zeros(numsim,1);

% Constraints
Yamax = zeros(numsim,1);
Xamax = zeros(numsim,1);

%%% Initial values
% Environmental variables
S(1) = S0;   
tempinc(1) = max(0, 3 * log2(max((1120 - S(1)), 1e-10) / 280));
Omega(1) = 1- (psi1*tempinc(1) + psi2*(tempinc(1)^2));

% Productivities
A_c(1)=(1+gamma*eta_c* s_c(1)) * Ac0;
A_d(1)=(1+gamma*eta_d* s_d(1)) * Ad0;

icl(1) = A_c(1) >= A_d(1);
idr(1) = A_d(1) > A_c(1);

A_a(1) = (1 + gamma * eta_a * s_a(1)) * Aa0;

% Intermediate good production evolution
numerc(1) = (A_c(1)^(-phi) + (1 + tau_d(1))^(1 - epsilon) * A_d(1)^(-phi) + ...
         (1 + tau_a(1))^(1 - epsilon) * Omega(1)^(1 - epsilon) * A_a(1)^(-phi));

denomc(1) = A_c(1)^(-phi) + (1 + tau_d(1))^(-epsilon) * A_d(1)^(-phi) + ...
        (1 + tau_a(1))^(-epsilon) * Omega(1)^(1 - epsilon) * A_a(1)^(-phi);

Yc(1) = (alpha / psi)^(alpha / (1 - alpha)) * Omega(1)^(1 / (1 - alpha)) * ...
     ((numerc(1))^(-alpha / phi) / denomc(1)) * (A_c(1))^(1 - alpha - phi);

numerd(1) = (A_c(1)^(-phi) + (1 + tau_d(1))^(1 - epsilon) * A_d(1)^(-phi) + ...
         (1 + tau_a(1))^(1 - epsilon) * Omega(1)^(1 - epsilon) * A_a(1)^(-phi));

denomd(1) = (1 + tau_d(1))^(epsilon) * A_c(1)^(-phi) + A_d(1)^(-phi) + ...
        ((1 + tau_d(1)) / (1 + tau_a(1)))^(epsilon) * Omega(1)^(1 - epsilon) * A_a(1)^(-phi);

Yd(1) = (alpha / psi)^(alpha / (1 - alpha)) * Omega(1)^(1 / (1 - alpha)) * ...
     ((numerd(1))^(-alpha / phi) / denomd(1)) * (A_d(1))^(1 - alpha - phi);

numera(1) = (Omega(1)^(-(1 - epsilon)) * A_c(1)^(-phi) + ...
         (1 + tau_d(1))^(1 - epsilon) * Omega(1)^(-(1 - epsilon)) * A_d(1)^(-phi) + ...
         (1 + tau_a(1))^(1 - epsilon) * (A_a(1))^(-phi));

denoma(1) = (1 + tau_a(1))^(epsilon) * Omega(1)^(-(1 - epsilon)) * A_c(1)^(-phi) + ...
        ((1 + tau_d(1)) / (1 + tau_a(1)))^(epsilon) * Omega(1)^(-(1 - epsilon)) * A_d(1)^(-phi) + ...
        A_a(1)^(-phi);

Yamax(1) = ((1 - Omega(1)) / Omega(1)) * (Yc(1) + Yd(1));

Ya(1) = (alpha / psi)^(alpha / (1 - alpha)) * ...
     ((numera(1))^(-alpha / phi) / denoma(1)) * A_a(1)^(1 - alpha - phi);

if Ya(1) > Yamax(1)
    Ya(1) = Yamax(1);
end

% Final good production evolution
Y(1) = (Yc(1)^((epsilon - 1) / epsilon) + Yd(1)^((epsilon - 1) / epsilon) + Ya(1)^((epsilon - 1) / epsilon))^(epsilon / (epsilon - 1));

% Machine production evolution
xcit(1) = (alpha / psi)^(alpha / (1 - alpha)) * Omega(1)^(1 / (1 - alpha)) * ...
     ((numerc(1))^(-1/ phi) / denomc(1)) * A_c(1)^(- phi);
xdit(1) = (alpha / psi)^(alpha / (1 - alpha)) * Omega(1)^(1 / (1 - alpha)) * ...
     ((numerd(1))^(-1/ phi) / denomd(1)) * A_d(1)^(- phi);
xait(1) = (alpha / psi)^(alpha / (1 - alpha)) * Omega(1)^(1 / (1 - alpha)) * ...
     ((numera(1))^(-1 / phi) / denoma(1)) * A_a(1)^(- phi);
Xamax(1) = ((1 - Omega(1)) / Omega(1)) * (xcit(1) + xdit(1));

if xait(1) > Xamax(1)
    xait(1) = Xamax(1);
end

% Consumption evolution
C(1) = Y(1) - psi*(xcit(1) + xdit(1) + xait(1));

%%%%%%%%%%%%%%%%
%%% Simulations
%%%%%%%%%%%%%%%%

for n = 2:numsim

% Update productivity
A_c(n) = (1 + gamma * eta_c * s_c(n)) * A_c(n-1);
A_d(n) = (1 + gamma * eta_d * s_d(n)) * A_d(n-1);

% Create indicators
icl(n) = A_c(n) >= A_d(n);
idr(n) = A_d(n) > A_c(n);

% Adaptive productivity 
A_a(n) = (1 + gamma * eta_a * s_a(n)) * A_a(n-1);

% Environment evolution
S(n) = min(max(0.00000000000000001,-qsi * (Yd(n-1) + idr(n) * Ya(n-1)) + (1 + delta) * S(n-1)), 1120); %notice that S should not be very close to zero

% Temperature and damage
tempinc(n) = 3 * log2(max(1e-6, (1120 - S(n))) / 280);
Omega(n) = max(1e-6, 1 - (psi1 * tempinc(n) + psi2 * tempinc(n)^2));
Omega(n) = max(1e-6, Omega(n));

% Intermediate outputs
numerc(n) = (A_c(n)^(-phi) + (1 + tau_d(n))^(1 - epsilon) * A_d(n)^(-phi) + ...
                 (1 + tau_a(n))^(1 - epsilon) * Omega(n)^(1 - epsilon) * A_a(n)^(-phi));
denomc(n) = A_c(n)^(-phi) + (1 + tau_d(n))^(-epsilon) * A_d(n)^(-phi) + ...
                (1 + tau_a(n))^(-epsilon) * Omega(n)^(1 - epsilon) * A_a(n)^(-phi);
Yc(n) = (alpha / psi)^(alpha / (1 - alpha)) * Omega(n)^(1 / (1 - alpha)) * ...
            ((numerc(n))^(-alpha / phi) / denomc(n)) * A_c(n)^(1 - alpha - phi);

numerd(n) = (A_c(n)^(-phi) + (1 + tau_d(n))^(1 - epsilon) * A_d(n)^(-phi) + ...
                 (1 + tau_a(n))^(1 - epsilon) * Omega(n)^(1 - epsilon) * A_a(n)^(-phi));
denomd(n) = (1 + tau_d(n))^(epsilon) * A_c(n)^(-phi) + A_d(n)^(-phi) + ...
                ((1 + tau_d(n)) / (1 + tau_a(n)))^(epsilon) * Omega(n)^(1 - epsilon) * A_a(n)^(-phi);
Yd(n) = (alpha / psi)^(alpha / (1 - alpha)) * Omega(n)^(1 / (1 - alpha)) * ...
            ((numerd(n))^(-alpha / phi) / denomd(n)) * A_d(n)^(1 - alpha - phi);

Yamax(n) = ((1 - Omega(n)) / Omega(n)) * (Yc(n) + Yd(n));

numera(n) = (Omega(n)^(-(1 - epsilon)) * A_c(n)^(-phi) + ...
                 (1 + tau_d(n))^(1 - epsilon) * Omega(n)^(-(1 - epsilon)) * A_d(n)^(-phi) + ...
                 (1 + tau_a(n))^(1 - epsilon) * A_a(n)^(-phi));
denoma(n) = (1 + tau_a(n))^(epsilon) * Omega(n)^(-(1 - epsilon)) * A_c(n)^(-phi) + ...
                ((1 + tau_d(n)) / (1 + tau_a(n)))^(epsilon) * Omega(n)^(-(1 - epsilon)) * A_d(n)^(-phi) + ...
                A_a(n)^(-phi);
Ya(n) = (alpha / psi)^(alpha / (1 - alpha)) * ...
            ((numera(n))^(-alpha / phi) / denoma(n)) * A_a(n)^(1 - alpha - phi);

if Ya(n) > Yamax(n)
    Ya(n) = Yamax(n);
end

% Final output
Y(n) = (Yc(n)^((epsilon - 1) / epsilon) + ...
    Yd(n)^((epsilon - 1) / epsilon) + ...
    Ya(n)^((epsilon - 1) / epsilon))^(epsilon / (epsilon - 1));

% Machine inputs
xcit(n) = (alpha / psi)^(alpha / (1 - alpha)) * Omega(n)^(1 / (1 - alpha)) * ...
            ((numerc(n))^(-1/ phi) / denomc(n)) * A_c(n)^(- phi);
xdit(n) = (alpha / psi)^(alpha / (1 - alpha)) * Omega(n)^(1 / (1 - alpha)) * ...
            ((numerd(n))^(-1 / phi) / denomd(n)) * A_d(n)^(- phi);
xait(n) = (alpha / psi)^(alpha / (1 - alpha)) * ...
            ((numera(n))^(-1 / phi) / denoma(n)) * A_a(n)^(- phi);
Xamax(n) = ((1 - Omega(n)) / Omega(n)) * (xcit(n) + xdit(n));

if xait(n) > Xamax(n)
    xait(n) = Xamax(n);
end

% Consumption
C(n) = Y(n) - psi * (xcit(n) + xdit(n) + xait(n));
end

Teste1 = (1+repmat(rho,numsim,1));
Teste2 = Teste1.^((0:numsim-1)');
Teste3 = 1./Teste2;

Util = zeros(numsim,1);

for j = 1:numsim
    try
        % Environmental utility modifier
        env_impact = 2 * sqrt(5 * min(S(j), 5)) - min(S(j), 5);
        utility_input = env_impact * C(j);

        % Robust check
        if utility_input > 1e-10 && isreal(utility_input)
            Util(j) = -(1 / (1 - sigma)) * Teste3(j) * utility_input^(1 - sigma) ...
                      - tax_penalty;
        else
            Util(j) = -1e10;
        end
    end
end
U=sum(Util(1:end));




