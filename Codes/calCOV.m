function K=calCOV(X1,X2,T,L,Kernel)
%This function finds the covariance matrix using the given coordinates and
%parameters.
%
%INPUTS
% X1: test or incuding point coordinates (nm)
% X2: test or inducing point coordinates (nm)
% T: correlation function parameter
% L: correlation function parameter
% Kernel: string specifying kernel type
%
%OUTPUTS
% K: correlation matrix
%
%Author:
%   Mohamadreza Fazel, Presse lab, 2024
%

if nargin < 5
    Kernel = 'Exponential'; 
end
Dist = pdist2(X1,X2);
if strcmp(Kernel,'Exponential')
    K =(T^2).*exp(-Dist.^2/L^2/2);
elseif strcmp(Kernel,'Quadratic')
    K =(T^2).*(1+Dist.^2/L^2/2);
end

end
