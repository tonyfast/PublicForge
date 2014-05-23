function AK=AIC(Res,Coeff)
% Inputs are residuals and coefficients from a regressions session. Calculates Akaike Information Criterion.
n=size(Res,1);
RSS=sum(Res.^2);
K=length(Coeff);
AK=n*log(RSS/n)+2*K;
