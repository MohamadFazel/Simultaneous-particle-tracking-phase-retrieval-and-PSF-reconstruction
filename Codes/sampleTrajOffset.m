function [X,Y,Z,PSFstack,Norm1_PSFstack,AcceptX] = sampleTrajOffset(Data,Chain,PSFstack,BNP,...
              Struct,DefocusK,Mask,DelX,AcceptX,XOffsetPhase,YOffsetPhase,...
              tmpPhase,SubPixelZeros,StartInd,EndInd,SubPixel,Tform,SigConv,Temp,CMOS_Noise)
%This function samples the x-trajectory, y-trajectory and z-trajectory
%
%INPUTS
% Data: the input data is a 4D array with the 3rd abd 4th dimension being
%       the frames and the planes.
% Chain: Chain of samples (for decription see "runPhaseRetrieval" help)
% PSFstack: the current model 4D array with a size similar to the data
% BNP: structure containing parameters used in the algorithm (for decription
%      see "runPhaseRetrieval" help)
% Struct: Structure containing parameters of the experiment (see PSFstruct
%         in "runPhaseRetrieval" help)
% DefoucusK: phase due to one nanometer offset with respect to focus
% Mask: one over the frequency range passed by the objective and zero otherwise
% DelX: axial location difference between the planes (first plane is
%       reference) (nm)
% AcceptX: number of accepted proposals for trajectories
% XOffsetPhase: phase due to one nanometer movement along the x-axis
% YOffsetPhase: phase due to one nanometer movement along the y-axis
% tmpPhase: pupil phase after subracting the location contribution
% SubPixelZeros: a zero frame with the subpixel size used in zero-padding
% StartInd: starting index used in zero-padding
% EndInd: end index used in zero padding
% SubPixel: number of model subpixels within a data pixel
% Tform: a 3x3 matrix of affine transform used in plane registeration
% SigConv: Sigma of the Gaussian to smotthen the model
% Temp: temperature for temporal sampling
% CMOS_Noise: pixel-map of CMOS camera noise (zero)
%
%OUTPUS:
% X: updated x-trajectory
% Y: updated y-trajectory
% Z: updated z-trajectory
% PSFstack: updated model using new phase and magnitude
% Norm1_PSFstack: Model without particle intensity and background
% AcceptX: number of accepted proposals
%
%Author:
%   Mohamadreza Fazel, Presse lab, 2024
%

Mag = Chain.Mag;
Bg = Chain.Bg;
I = Chain.I;
Sig = 1;
X = Chain.X;
Y = Chain.Y;
Z = Chain.Z;
D = Chain.D;
Dt = BNP.Dt;
Sig_D = sqrt(2*D*Dt);

tmp = rand();
if tmp < 1/3
    tX = X;
    tY = Y;
    tZ = Z + Sig*randn();
else %tmp> 1/4 && tmp < 1/2
    Ind = ceil(size(Data,3)*rand());
    tX = X; tX(Ind) = X(Ind) + 10*Sig*randn();
    tY = Y; tY(Ind) = Y(Ind) + 10*Sig*randn();
    tZ = Z; tZ(Ind) = Z(Ind) + 5*Sig*randn();
% elseif tmp > 3/4
%     Ind = ceil(size(Data,3)*rand());
%     tX = X; tX(Ind) = X(Ind) + 40*Sig*randn();
%     tY = Y; 
%     tZ = Z; 
% else
%     Ind = ceil(size(Data,3)*rand());
%     tZ = Z; tZ(Ind) = Z(Ind) + 20*Sig*randn();
%     tY = Y;
%     tX = X;
end

tPSF = [];
Norm1_PSFstack = [];
for ii = 1:size(DelX,1)
    [BgPSF,tNorm1_PSF]=findPSF(Mag,tmpPhase,Bg(ii),I(ii),DefocusK,tZ+DelX(ii,3),Mask,...
        SubPixelZeros,StartInd,EndInd,SubPixel,tX+DelX(ii,1),tY+DelX(ii,2),...
        XOffsetPhase,YOffsetPhase,Tform(ii),SigConv,CMOS_Noise(:,:,ii));
    tPSF = cat(4,tPSF,BgPSF);
    Norm1_PSFstack = cat(4,Norm1_PSFstack,tNorm1_PSF);
end
DLogL = sum(Data(:).*(log(tPSF(:))-log(PSFstack(:)))-(tPSF(:)-PSFstack(:)));
DLogL = DLogL/Temp;

SigX = Struct.NPix*Struct.PixelSize/2;
DLPrior1X = sum(log(normpdf(tX(1),0,SigX)) ...
    -log(normpdf(X(1),0,SigX)));
DLPrior1Y = sum(log(normpdf(tY(1),0,SigX)) ...
    -log(normpdf(Y(1),0,SigX)));
DLPrior1Z = sum(log(normpdf(tZ(1),0,SigX)) ...
    -log(normpdf(Z(1),0,SigX)));

DLPriorX = sum(log(normpdf(tX(2:end),tX(1:end-1),Sig_D)) ...
    -log(normpdf(X(2:end),X(1:end-1),Sig_D)));
DLPriorY = sum(log(normpdf(tY(2:end),tY(1:end-1),Sig_D)) ...
    -log(normpdf(Y(2:end),Y(1:end-1),Sig_D)));
DLPriorZ = sum(log(normpdf(tZ(2:end),tZ(1:end-1),Sig_D)) ...
    -log(normpdf(Z(2:end),Z(1:end-1),Sig_D)));

DLogPost = DLogL + DLPrior1X + DLPrior1Y + DLPrior1Z + DLPriorX + DLPriorY + DLPriorZ;

if DLogPost > log(rand())
    X = tX;
    Y = tY;
    Z = tZ;
    PSFstack = tPSF;
    AcceptX = AcceptX + 1;
end

DEBUG = 0;
if DEBUG
   figure(210);hold
   plot(Struct.X,'o-')
   plot(X,'o-');ylabel('X')
   figure(211);hold
   plot(Struct.Y,'o-')
   plot(Y,'o-');ylabel('Y')
   figure(212);hold
   plot(Struct.Z,'o-')
   plot(Z,'o-');ylabel('Z')
   pause(0.2)
end          
          
end 