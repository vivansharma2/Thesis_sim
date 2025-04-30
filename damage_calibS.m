% this function computes the difference according to the norm L2 between Nordhaus's function of environmental quality and ours for different values of lambda, on the interval of temperature increase from 0 to t_max.
function d = damage_calibS(lam)
global t_disaster max_t S_bar
S_min = 280 * 2^(t_disaster/3) - 280*2^(max_t/3);
dd=0;
pas=(S_bar-S_min)/1000;
for S = S_bar - pas/2:-pas:S_min+pas/2
    t = 3 * log((280 * 2^(t_disaster/3) - S)/280)/log(2);
    dd=dd+pas*(1/(1+0.00284*t^2)-(((t_disaster-t)^lam - lam * (t_disaster ^(lam -1)) * (t_disaster -t))/((1-lam)*(t_disaster^lam))))^2; 
end
d=dd;