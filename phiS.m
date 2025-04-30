% ******This is the function phi(S) *****
function pphiS = phiS(S)
global lambda t_disaster
t = 3 * log((280 * 2^(t_disaster/3) - S)/280)/log(2);
pphiS = ((t_disaster - t)^lambda - lambda * t_disaster ^(lambda -1) * (t_disaster -t))/((1-lambda)*t_disaster ^lambda);