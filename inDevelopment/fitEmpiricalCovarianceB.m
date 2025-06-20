function [EMP_COV,psi,sigma2,l]=fitEmpiricalCovarianceB(psi,EMP_COV)
%% Set coefficients
sigma2=1;
l=0.001;
for k=1:250
% Compute the function
C2=sigma2*exp(-(psi.^2)/(2*l^2));    
C2_dsigma2=exp(-(psi.^2)/(2*l^2));
C2_dl=(psi.^2)./(l^3).*sigma2.*exp(-(psi.^2)/(2*l^2));
% Compute best fitted coefficients
if k==1
sigma2=C2_dsigma2\EMP_COV;
else
Sol=[C2_dsigma2(1:5),C2_dl(1:5)]\(EMP_COV(1:5)-C2(1:5));
sigma2=abs(sigma2+0.1*Sol(1));
l=l+0.2*Sol(2);
end
figure(16)
clf
hold on
plot(psi,EMP_COV,'*')
plot(psi,C2)
title([num2str(sigma2),' ',num2str(l),' ',num2str(std(((EMP_COV-C2)))),' ',num2str(k)])
end
end