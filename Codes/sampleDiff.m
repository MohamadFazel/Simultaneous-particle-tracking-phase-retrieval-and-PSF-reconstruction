function D = sampleDiff(Chain,DelX,Dt)
%This function samples diffusion coefficient using direct sampling
%
%INPUTS
% Chain: Chain of samples (for decription see "runPhaseRetrieval" help)
% DelX: axial location difference between the planes (first plane is
%       reference) (nm)
% Dt: frame exposure time (ns)
%
%OUTPUTS
% D: sampled diffusion coefficient
%
%Author:
%   Mohamadreza Fazel, Presse lab, 2024
%

AlphaPrior = 0.2;
BetaPrior = 20;
X = Chain.X;
Y = Chain.Y;
Z = Chain.Z;
R = size(DelX,1);

Alpha = 3*R*(length(X)-1)+AlphaPrior;
Beta = BetaPrior + R*sum((X(2:end)-X(1:end-1)).^2 ...
        +(Y(2:end)-Y(1:end-1)).^2 + (Z(2:end)-Z(1:end-1)).^2)/(2*Dt); 
D = 1./gamrnd(Alpha,1/Beta);

end