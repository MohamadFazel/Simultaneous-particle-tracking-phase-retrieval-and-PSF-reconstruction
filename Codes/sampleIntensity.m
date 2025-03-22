function [I,PSFstack,AcceptI]=sampleIntensity(Data,Chain,PSFstack,Norm1_PSFstack,DelX,AcceptI,CMOS_Noise)
%This function samples particle intensity (number of photons)
%
%INPUTS
% Data: the input data is a 4D array with the 3rd abd 4th dimension being
%       the frames and the planes.
% Chain: Chain of samples (for decription see "runPhaseRetrieval" help)
% PSFstack: the current model 4D array with a size similar to the data
% Norm1_PSFstack: PSF model without background and particle intensity
% DelX: axial location difference between the planes (first plane is
%       reference) (nm)
% AcceptI: number of accepted proposals for particle intensity
% CMOS_Noise: pixel-map of CMOS camera noise (zero)
% 
%OUTPUTS
% I: updated particle intensity
% PSFstack: updated model using new particle intnsity
% AcceptI: updated number of accepted proposed I
%
%Author:
%   Mohamadreza Fazel, Presse lab, 2024
%

Bg = Chain.Bg;
I = Chain.I;
Alpha_Prop = 20000;
Alpha = 5;
Beta = 2000;

Ind = randi(size(I,1));
tI = I;
tI(Ind) = gamrnd(Alpha_Prop,I(Ind)/Alpha_Prop);
tPSF = [];
for ii = 1:size(DelX,1)
    BG = CMOS_Noise(:,:,ii) + Bg(ii);
    tPSF = cat(4,tPSF,tI(ii)*Norm1_PSFstack(:,:,:,ii)+BG);
end
DLogL = sum(Data(:).*(log(tPSF(:))-log(PSFstack(:)))-(tPSF(:)-PSFstack(:)));
DLogL = DLogL/10;

DLogPrior = sum(log(gampdf(tI,Alpha,Beta)) - log(gampdf(I,Alpha,Beta)));

DLogProp = sum(log(gampdf(I,Alpha_Prop,tI/Alpha_Prop)) ...
    - log(gampdf(tI,Alpha_Prop,I/Alpha_Prop)));

DLogPost = DLogL + DLogPrior + DLogProp;

if DLogPost > log(rand())
    I = tI;
    PSFstack = tPSF;
    AcceptI = AcceptI + 1;
end

end