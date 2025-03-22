function [SigConv,PSFstack,AcceptSig] = sampleSigConv(Data,Chain,PSFstack,...
    DefocusK,Mask,DelX,AcceptSig,XOffsetPhase,YOffsetPhase,tmpPhase,...
    SubPixelZeros,StartInd,EndInd,SubPixel,Tform,SigConv,CMOS_Noise)
%This function samples the sigma of the Gaussian filter used to smoothen
%the PSF (see supplementary Note 1).
%
%INPUTS
% Data: the input data is a 4D array with the 3rd abd 4th dimension being
%       the frames and the planes.
% Chain: Chain of samples (for decription see "runPhaseRetrieval" help)
% PSFstack: the current model 4D array with a size similar to the data
% DefoucusK: phase due to one nanometer offset with respect to focus
% Mask: one over the frequency range passed by the objective and zero otherwise
% DelX: axial location difference between the planes (first plane is
%       reference) (nm)
% AcceptSig: number of accepted proposals for this parameter
% XOffsetPhase: phase due to one nanometer movement along the x-axis
% YOffsetPhase: phase due to one nanometer movement along the y-axis
% tmpPhase: pupil phase after subracting the location contribution
% SubPixelZeros: a zero frame with the subpixel size used in zero-padding
% StartInd: starting index used in zero-padding
% EndInd: end index used in zero padding
% SubPixel: number of model subpixels within a data pixel
% Tform: a 3x3 matrix of affine transform used in plane registeration
% SigConv: Sigma of the Gaussian to smotthen the model
% CMOS_Noise: pixel-map of CMOS camera noise (zero)
% 
%OUTPUT
% SigConv: updated value for SigConv
% PSFstack: updated model using the new value for SigConv
% AcceptSig: number of accepted proposals
%
%Author:
%   Mohamadreza Fazel, Presse lab, 2024
%

Mag = Chain.Mag;
Bg = Chain.Bg;
I = Chain.I;
X = Chain.X;
Y = Chain.Y;
Z = Chain.Z;

% proposing new value
Alpha = 10000;
tmpSigConv = gamrnd(Alpha,SigConv/Alpha);

%generating model based on the proposed value above
tPSF = [];
for ii = 1:size(DelX,1)
    tPSF = cat(4,tPSF,findPSF(Mag,tmpPhase,Bg(ii),I(ii),DefocusK,Z+DelX(ii,3),Mask,...
        SubPixelZeros,StartInd,EndInd,SubPixel,X+DelX(ii,1),Y+DelX(ii,2),...
        XOffsetPhase,YOffsetPhase,Tform(ii),tmpSigConv,CMOS_Noise(:,:,ii)));
end
%caclulating log-likelihood ratio
DLogL = sum(Data(:).*(log(tPSF(:))-log(PSFstack(:)))-(tPSF(:)-PSFstack(:)));
%calculating log-prior ratio
PropLogL = log(gampdf(SigConv,Alpha,tmpSigConv/Alpha)) - log(gampdf(tmpSigConv,Alpha,SigConv/Alpha));
%calculating log proposal dist. ratio
PLogL = log(gampdf(tmpSigConv,1.5,1)) - log(gampdf(SigConv,1.5,1));

if DLogL + PLogL + PropLogL > log(rand())
    SigConv = tmpSigConv;
    PSFstack = tPSF;
    AcceptSig = AcceptSig + 1;
end

end