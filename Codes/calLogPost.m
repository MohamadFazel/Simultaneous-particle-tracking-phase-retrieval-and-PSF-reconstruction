function [LogPost,LogLike] = calLogPost(Data,Chain,PSFstack,Struct,Chol_A,Chol_Phi,BNP)
%This function calculates log-posterior
%
%INPUTS
% Data: the input data is a 4D array with the 3rd abd 4th dimension being
%       the frames and the planes.
% Chain: Chain of samples (for decription see "runPhaseRetrieval" help)
% PSFstack: the current model 4D array with a size similar to the data
% Struct: Structure containing parameters of the experiment (see PSFstruct
%         in "runPhaseRetrieval" help)
% Chol_A: Choleski decomposition of magnitude correlation kernel
% Chol_Phi: Choleski decomposition of phase correlation kernel
% BNP: structure containing parameters used in the algorithm (for decription
%      see "runPhaseRetrieval" help)
% 
%OUTPUTS
% LogPost: log-posterior
% LogLike: log-likelihood
%
%Author:
%   Mohamadreza Fazel, Presse lab, 2024
%


Mag = Chain.Mag;
Phase = Chain.Phase;
X = Chain.X;
Y = Chain.Y;
Z = Chain.Z;
SigX = Struct.NPix*Struct.PixelSize/2;
D = Chain.D; %BNP.D;
Dt = BNP.Dt;
Sig_D = sqrt(2*D*Dt);

%Poisson likelihood with Stirling approximation
LogLike = sum(Data(:).*log(PSFstack(:))-(Data(:)+0.5).*log(Data(:))+Data(:)-PSFstack(:)-0.9189);

LogPriorMag = -0.5*sum(((log(Mag(:))-0)'/Chol_A).^2) ...
        - sum(log(diag(Chol_A))) - length(Mag(:))*log(2*pi)/2;
LogPriorPhase = -0.5*sum(((Phase(:)-0)'/Chol_Phi).^2) ...
        - sum(log(diag(Chol_Phi))) - length(Phase(:))*log(2*pi)/2;

LogPrior1X = sum(log(normpdf(X(1),0,SigX)));
LogPrior1Y = sum(log(normpdf(Y(1),0,SigX)));
LogPrior1Z = sum(log(normpdf(Z(1),0,SigX)));

LogPriorX = sum(log(normpdf(X(2:end),X(1:end-1),Sig_D)));
LogPriorY = sum(log(normpdf(Y(2:end),Y(1:end-1),Sig_D)));
LogPriorZ = sum(log(normpdf(Z(2:end),Z(1:end-1),Sig_D)));

LogPost = LogLike + LogPriorMag + LogPriorPhase + LogPrior1X + LogPrior1Y + ...
    LogPrior1Z + LogPriorX + LogPriorY + LogPriorZ;

end