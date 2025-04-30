%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR: Vivan Sharma
% TASK: Thesis model simulation
% DATE MODIFIIED: 30 Apr 2025
%%%%%%%%%%%%%%%%%%%%%%%%%

%********************************************************************
% 1. INITIAL PARAMETERS
%********************************************************************
global dt rho sigma psi alpha gamma eta_d eta_c eta_a qsi epsilon delta numsim phi S_bar lambda t_disaster max_t psi1 psi2 

%==== Time and discounting ====%
dt          = 5; % number of years in a period
rho         = 0.001 * dt; % discount rate
epsilon     = 10; % elasticity of substitution (>1)
sigma       = 2; % relative risk aversion 

%==== Production and innovation parameters ====%
alpha = 1/3;% 1/3; % share of machines in production
psi= alpha^2; % cost of machines
gamma = 1; % size of innovation
eta_d = 0.02*dt; % probability of success in the dirty sector
eta_c = 0.02*dt; % probability of success in the the clean sector
eta_a = 0.02*dt; % probability of success in the the adaptive sector

%==== Simulation horizon ====%
numsim=80; % number of periods

%==== Production (used for calculating initial productivities) ====%
Yc0=307.77; % production of non fossil fuel energy in the world primary supply of energy from 2002 to 2006 source in Quadrillion of Btu 
Yd0=1893.25; % production of fossil fuel energy in the world primary supply of energy from 2002 to 2006 source in Quadrillion of Btu

%==== Emissions and damage ====%
emission0=17.48251782; % world emissions in ppm from 2002 to 2006 (converted from original data in GtCO2)
t_disaster=6;% disaster temperature
S_bar = 280*(2^(t_disaster/3)-1); % corresponding value for S_bar
S0 = S_bar - 99; % 99pm increase in CO2 concentration since preindustrial times
max_t=3; % temperature up to which the utility damage function is matched with Nordhaus's one
lambda = fmincon(@(l)damage_calibS(l),0.35,[],[],[],[],0.00001,0.999999,[],optimset('Tolfun',1e-11));
delta =.5*emission0/S0; % half of CO2 emissions are absorbed at current atmospheric levels

%==== Damage function - change as you wish ====%
psi1 = 0.05;
psi2 = 0.04;
tempinc0 =  3*log2((1120 - S0)/ 280); % increase in temperature since preindustrial times 
loss0 = psi1*tempinc0 + psi2*(tempinc0^2); % 0.58% of output is lost this period
Omega0 = 1-loss0;

%==== Derived parameters and productivities ====%
phi=(1-alpha)*(1-epsilon);
qsi = emission0/Yd0;
Ac0=(alpha/psi)^(-alpha/(1-alpha))*Yc0*(1+(Yc0/Yd0)^(1/epsilon-1))^(alpha/phi+1); % initial value unchanged for simplicity
Ad0 =(alpha/psi)^(-alpha/(1-alpha))*Yd0*(1+(Yd0/Yc0)^(1/epsilon-1))^(alpha/phi+1); % and for Ad0
Aa0 = (Ac0+Ad0)/4; % initial assumption, change as you wish

%********************************************************************
% 2. OPTIONS
%********************************************************************
mode = 1; % 0 for no IDTC, 1 for IDTC
delay = 30 * mode;  % delay must be 0 if mode == 0

%==== Initializing delay stock vectors ====%
StockAc = [];
StockAd = [];
StockAa = [];
StockC  = [];
StockQ  = [];
StockS  = [];
StockOmega = [];
StockYc = [];
StockYd = []; 
StockYa = [];
Stocktau_d_all = [];  
Stocktau_a_all = [];
StockS_c = [];
StockS_d = [];
StockS_a = [];
utcomp=0;

numsimleft=numsim-delay;
%********************************************************************
% NO INTERVENTION NEEDED FURTHER
%********************************************************************

%********************************************************************
% 3. OPTIMIZATION
%********************************************************************

%==== Settings ====%
display_iter=0;    %=1 if you want to see iterations, ~1 if otherwise
display_diag=1;    %=1 if you want to see diagnostics, ~1 if otherwise

if display_iter==1
    if display_diag==1
    options = optimset('Display','iter', 'Diagnostics','on','FunValCheck','on','TolFun',1e-11,'TolX',1e-9);
    else
    options = optimset('Display','iter', 'Diagnostics','off', 'FunValCheck','on','TolFun',1e-11,'TolX',1e-9);
    end
else
    if display_diag==1
    options = optimset('Display','off', 'Diagnostics','on','FunValCheck','on','TolFun',1e-11,'TolX',1e-9);
    else 
    options = optimset('Display','off', 'Diagnostics','off', 'FunValCheck','on','TolFun',1e-11,'TolX',1e-9);
    end
end   

%==== Optimization ====%
init = [
    0.24;
    0.17;
    0.1;
    0.04;
    zeros(76, 1)
]; % Copying original

if mode == 0 
    td0=init; % initial guess for the input tax on the dirty input
    ta0=init; % initial guess for the input tax on the adaptive input
    x0 = [td0; ta0]; % stacking the taxes
    lb = zeros(2 * numsim, 1);                  % lower bounds 
    ub = 100000 * ones(2 * numsim, 1);          % upper bounds 
else 
    sc0=ones(numsim,1); % initial guess for the share of scientists in clean research
    sa0=zeros(numsim,1); % initial guess for the share of scientists in adaptive research
    sd0=zeros(numsim,1); % initial guess for the share of scientists in dirty research
    x0 = [sc0; sd0; sa0];
    lb = zeros(3 * numsim, 1);                  % lower bounds 
    ub = 100000 * ones(3 * numsim, 1);          % upper bounds
end 

%==== Scientist constraint - must sum to 1 ====%
n_vars = 3 * numsim;

Aeq = zeros(numsim, n_vars);
for t = 1:numsim
    Aeq(t, t) = 1;
    Aeq(t, t + numsim) = 1;
    Aeq(t, t + 2*numsim) = 1;
end
beq = ones(numsim, 1);

%==== With delay ====%
if delay>0
    numsim = delay;
    comp = [ ...
        eta_c/(eta_c+eta_d+eta_a)*ones(delay,1);   % clean scientist share 
        eta_d/(eta_c+eta_d+eta_a)*ones(delay,1);   % dirty scientist share 
        eta_a/(eta_c+eta_d+eta_a)*ones(delay,1);   % adaptive scientist share 
        zeros(delay,1);    
        zeros(delay,1); % stacked: [tau_d; tau_a]
    ]; 
    % comp stacks the allocation of scientists to the clean sector and the initial carbon tax
    RespC=noidtc(comp, Ac0, Ad0, Aa0, S0); % compute the relevant parameters for the economy
    StockAc=RespC.Ac;
    StockAd=RespC.Ad;
    StockAa=RespC.Aa;
    StockC=RespC.C;

    StockYc = RespC.Yc;
    StockYd = RespC.Yd;
    StockYa = RespC.Ya;

    StockS_c = (eta_c/(eta_c+eta_d+eta_a)) * ones(delay,1);
    StockS_d = (eta_d/(eta_c+eta_d+eta_a)) * ones(delay,1);
    StockS_a = (eta_a/(eta_c+eta_d+eta_a)) * ones(delay,1);
    Stocktau_d_all = zeros(delay,1);
    Stocktau_a_all=zeros(delay,1);

    Stockshare_Yc = RespC.share_Yc;
    Stockshare_Yd = RespC.share_Yd;
    Stockshare_Ya = RespC.share_Ya;
    Stocknumerc = RespC.numerc;
    Stocknumerd = RespC.numerd;
    Stocknumera = RespC.numera;
    Stockdenomc = RespC.denomc;
    Stockdenomd = RespC.denomd;
    Stockdenoma = RespC.denoma;

    StockOmega = RespC.Omega;
    StockQ=zeros(delay,1);
    StockS=RespC.S;
    utcomp=-sum(RespC.Util(1:end));
    S0 = RespC.S(delay);
    Ac0=StockAc(delay);
    Ad0=StockAd(delay);
    Aa0=StockAa(delay);
end

%==== Update ====%
numsim=numsimleft;


%==== Run ====%
if mode == 0
    [x,fval,exitflag] = fmincon(@(x)optnoidtc(x, Ac0, Ad0, Aa0, S0),x0,[],[],[],[],lb,ub,[],options); % no IDTC
    x = [ ...
        eta_c/(eta_c+eta_d+eta_a)*ones(numsim,1);   % clean scientist share 
        eta_d/(eta_c+eta_d+eta_a)*ones(numsim,1);   % dirty scientist share 
        eta_a/(eta_c+eta_d+eta_a)*ones(numsim,1);   % adaptive scientist share 
        0.1*ones(numsim,1);    
        0.1*ones(numsim,1); % stacked: [tau_d; tau_a]
    ]; % x now combines the (constant) allocation of scientists with the tax rate
else 
    [x,fval,exitflag] = fmincon(@(x)optidtc(x, Ac0, Ad0, Aa0, S0), ...
                            x0, [], [], Aeq, beq, lb, ub, [], options);
end 

%********************************************************************
% 3. RESULTS
%********************************************************************

%==== Stack results ====%
Util=utcomp-1/(1+rho)^delay*fval;
if mode == 0
    Resp = noidtc(x, Ac0, Ad0, Aa0, S0);
    
    % Combine past and simulated values for plotting or saving
    S_c_all = [StockS_c; Resp.S_c];
    S_d_all = [StockS_d; Resp.S_d];
    S_a_all = [StockS_a; Resp.S_a];   % adaptive (no stock assumed)
    
    Acc = [StockAc; Resp.Ac];         % clean productivity
    Add = [StockAd; Resp.Ad];         % dirty productivity
    Aaa = Resp.Aa;                    % adaptive productivity (no stock assumed)
    
    Yc_all = Resp.Yc;
    Yd_all = Resp.Yd;
    Ya_all = Resp.Ya;
    
    Omega = Resp.Omega;
    
    tau_d_all = [Stocktau_d_all; Resp.tau_d]; % dirty tax
    tau_a_all = [Stocktau_a_all; Resp.tau_a];  % adaptive tax
    
    Cc = [StockC; Resp.C];            % consumption
    Ss = [StockS; Resp.S];            % environmental quality

else 
    Resp = idtc(x, Ac0, Ad0, Aa0, S0);
    
    % Combine past and simulated values for plotting or saving
    S_c_all = [StockS_c; Resp.S_c];    % share of scientists in clean
    S_d_all = [StockS_d; Resp.S_d];    % dirty (no stock assumed)
    S_a_all = [StockS_a; Resp.S_a];    % adaptive (no stock assumed)
    
    Acc = [StockAc; Resp.Ac];         % clean productivity
    Add = [StockAd; Resp.Ad];         % dirty productivity
    Aaa = [StockAa; Resp.Aa];         % adaptive productivity (no stock assumed)
    
    Omega = [StockOmega; Resp.Omega];
    share_Yc = [Stockshare_Yc; Resp.share_Yc];
    share_Yd = [Stockshare_Yd; Resp.share_Yd];
    share_Ya = [Stockshare_Ya; Resp.share_Ya];
    tau_d_all = [Stocktau_d_all; Resp.tau_d]; % dirty tax
    tau_a_all = [Stocktau_a_all; Resp.tau_d]; % adaptive tax

    Numerc = [Stocknumerc; Resp.numerc];
    Numerd = [Stocknumerd; Resp.numerd];
    Numera = [Stocknumera; Resp.numera]; 

    Denomc = [Stockdenomc; Resp.denomc];
    Denomd = [Stockdenomd; Resp.denomd];
    Denoma = [Stockdenoma; Resp.denoma];

    Qq = [StockQ; Resp.Q];
    Cc = [StockC; Resp.C];            % consumption
    Ss = [StockS; Resp.S];            % environmental quality

end 

%********************************************************************
% 5. PLOTTING
%********************************************************************
numsimc = numsim + delay;
if mode == 0
figure;

% 1. Productivity and consumption
subplot(2,2,1)
plot(1:numsimc, Acc(1:numsimc), 'LineWidth', 1.5); hold on;
plot(1:numsimc, Add(1:numsimc), 'LineWidth', 1.5);
plot(1:numsimc, Aaa(1:numsimc), 'LineWidth', 1.5);
xlabel('Period')
ylabel('Economy Metrics')
legend('A_c', 'A_d', 'A_a')
title('Productivity')
grid on

% 2. Temperature path
subplot(2,2,2)
plot(1:numsimc, Omega(1:numsimc), 'r', 'LineWidth', 1.5);
xlabel('Period')
ylabel('Proportion of output saved')
legend('\Delta')
title('Remaining output')
grid on

% 3. Production
subplot(2,2,[3 4]) % use the bottom full row
plot(1:numsimc, Yc_all(1:numsimc), 'b', 'LineWidth', 1.5); hold on;
plot(1:numsimc, Yd_all(1:numsimc), 'r', 'LineWidth', 1.5);
plot(1:numsimc, Ya_all(1:numsimc), 'g', 'LineWidth', 1.5);
plot(1:numsimc, Cc(1:numsimc), 'k--', 'LineWidth', 1.5);
xlabel('Period')
ylabel('Intermediate Outputs')
legend('Y_c', 'Y_d', 'Y_a', 'C')
title('Sector Output Over Time')
grid on

else
figure;

% 1. Productivity
subplot(3,2,1)
plot(1:numsimc, Acc(1:numsimc), 'LineWidth', 1.5); hold on;
plot(1:numsimc, Add(1:numsimc), 'LineWidth', 1.5);
plot(1:numsimc, Aaa(1:numsimc), 'LineWidth', 1.5);
xlabel('Period')
ylabel('Productivity')
legend('A_c', 'A_d', 'A_a')
title('Productivity')
grid on

% 2. Temperature/Omega
subplot(3,2,2)
plot(1:numsimc, Omega(1:numsimc), 'r', 'LineWidth', 1.5);
xlabel('Period')
ylabel('Remaining Output')
legend('\Omega')
title('Damage Function (1 - Loss)')
grid on

% 3. Production
subplot(3,2,3)
plot(1:numsimc, share_Yc(1:numsimc), 'b', 'LineWidth', 1.5); hold on;
plot(1:numsimc, share_Yd(1:numsimc), 'r', 'LineWidth', 1.5);
plot(1:numsimc, share_Ya(1:numsimc), 'g', 'LineWidth', 1.5);
xlabel('Period')
ylabel('Output')
legend('Y_c', 'Y_d', 'Y_a')
title('Share of Sector Output')
grid on

% 4. Consumption
subplot(3,2,4)
plot(1:numsimc, Cc(1:numsimc), 'k--', 'LineWidth', 1.5);
xlabel('Period')
ylabel('Consumption')
title('Final Consumption Over Time')
grid on

% 5. Scientist shares
subplot(3,2,5)
plot(1:numsimc, S_c_all(1:numsimc), 'b', 'LineWidth', 1.5); hold on;
plot(1:numsimc, S_d_all(1:numsimc), 'r', 'LineWidth', 1.5);
plot(1:numsimc, S_a_all(1:numsimc), 'g', 'LineWidth', 1.5);
xlabel('Period')
ylabel('Scientist Share')
legend('S_c', 'S_d', 'S_a')
title('Scientist Allocation')
grid on

% 6. Subsidy Q
subplot(3,2,6)
plot(1:numsimc, Qq(1:numsimc), 'm', 'LineWidth', 1.5);
xlabel('Period')
ylabel('Subsidy to Clean')
title('Subsidy Path')
legend('Q')
grid on
end 