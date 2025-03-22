function [Bg,PSFstack,AcceptBg] = sampleBg(Data,Chain,PSFstack,Norm1_PSFstack,DelX,AcceptBg,CMOS_Noise)
%This function samples background (number of photons)
%
%INPUTS
% Data: the input data is a 4D array with the 3rd abd 4th dimension being
%       the frames and the planes.
% Chain: Chain of samples (for decription see "runPhaseRetrieval" help)
% PSFstack: the current model 4D array with a size similar to the data
% Norm1_PSFstack: PSF model without background and particle intensity
% DelX: axial location difference between the planes (first plane is
%       reference) (nm)
% AcceptBg: number of accepted proposals for background
% CMOS_Noise: pixel-map of CMOS camera noise (zero)
% 
%OUTPUTS
% Bg: updated background
% PSFstack: updated model using new background
% AcceptBg: updated number of accepted proposed Bg
%
%Author:
%   Mohamadreza Fazel, Presse lab, 2024
%

Bg = Chain.Bg;
I = Chain.I;
Alpha_Prop = 25000;
Alpha = 1;
Beta = 120;

tBg = gamrnd(Alpha_Prop,Bg/Alpha_Prop);
tPSF = [];
for ii = 1:size(DelX,1)
    BG = CMOS_Noise(:,:,ii)+tBg(ii);
    tPSF = cat(4,tPSF,I(ii)*Norm1_PSFstack(:,:,:,ii)+BG);
end
DLogL = sum(Data(:).*(log(tPSF(:))-log(PSFstack(:)))-(tPSF(:)-PSFstack(:)));
DLogL = DLogL/10;

DLogPrior = sum(log(gampdf(tBg,Alpha,Beta)) - log(gampdf(Bg,Alpha,Beta)));

DLogProp = sum(log(gampdf(Bg,Alpha_Prop,tBg/Alpha_Prop)) ...
    - log(gampdf(tBg,Alpha_Prop,Bg/Alpha_Prop)));

DLogPost = DLogL + DLogPrior + DLogProp;

if DLogPost > log(rand())
    Bg = tBg;
    PSFstack = tPSF;
    AcceptBg = AcceptBg + 1;
end

end