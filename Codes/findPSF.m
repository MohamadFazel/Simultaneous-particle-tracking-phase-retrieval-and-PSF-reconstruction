function [BgPSF,PSF] = findPSF(Pupil_Mag,Pupil_Phase,Bg,I,DefocusK,Z,Mask,...
    SubPixelZeros,StartInd,EndInd,SubPixel,X,Y,XOffsetPhase,YOffsetPhase,Tform,SigConv,CMOS_Noise)
%findPSF generates PSF using the input pupil phase and magnitude, trajectory 
% and other paraneters.
%
%INPUTS:
% Pupil_Mag: pupil magnitude
% Pupil_Phase: pupil phase
% Bg: uniform background photon count for each plane
% I: number of photons per particle in each plane
% DefoucusK: phase due to one nanometer offset with respect to focus
% Z: z-trajectory (nm)
% Mask: one over the frequency range passed by the objective and zero otherwise
% SubPixelZeros: a zero frame with the subpixel size used in zero-padding
% StartInd: starting index used in zero-padding
% EndInd: end index used in zero padding
% SubPixel: number of model subpixels within a data pixel
% X: x-trajectory (nm)
% Y: y-trajectory (nm)
% XOffsetPhase: phase due to one nanometer movement along the x-axis
% YOffsetPhase: phase due to one nanometer movement along the y-axis
% Tform: a 3x3 matrix of affine transform used in plane registeration
% SigConv: Sigma of the Gaussian to smotthen the model
% CMOS_Noise: pixel-map of CMOS camera noise (zero)
%
%OUTPUTS:
% BgPSF: PSF model at different frames
% PSF: PSF model normalized to one at each frame
%
%Author:
%   Mohamadreza Fazel, Presse lab, 2024
%
NPix = size(Mask,1);
PSF = zeros(NPix,NPix,length(Z),'single');
BgPSF = zeros(NPix,NPix,length(Z),'single');

if nargin < 12
    for zz = 1:length(Z)
        DefocusPhase = DefocusK*Z(zz);
        Phase = Mask.*(DefocusPhase+Pupil_Phase);
        %Optical transfer function (propagator)
        OTF = Mask.*Pupil_Mag.*exp(1i*Phase);
        
        %Parseval Normalizaton
        Norm = sqrt(sum(sum(abs(OTF))))*NPix;
        PSF(:,:,zz) = single((abs(fftshift(fft2(OTF/Norm)))).^2);
        PSF(:,:,zz) = PSF(:,:,zz)/sum(sum(PSF(:,:,zz)));
        BgPSF(:,:,zz) = I(zz)*PSF(:,:,zz)+Bg;
    end
else
    for zz = 1:length(Z)
        DefocusPhase = DefocusK*Z(zz);
        Phase = Mask.*(DefocusPhase+Pupil_Phase ...
            +XOffsetPhase*X(zz)+YOffsetPhase*Y(zz));
        %Optical transfer function (propagator)
        OTF = Mask.*Pupil_Mag.*exp(1i*Phase);
        %Parseval Normalizaton
        Norm = sqrt(sum(sum(abs(OTF))))*NPix;
        %Padding zeros for subpixel resolution
        SubPixelZeros(StartInd:EndInd,StartInd:EndInd) = OTF;
        tmpPSF = single((abs(fftshift(fft2(SubPixelZeros/Norm)))).^2);
        tPSF = blockresample(tmpPSF,1/SubPixel)/SubPixel;
        PSF(:,:,zz) = tPSF;
        PSF(:,:,zz) = PSF(:,:,zz)/sum(sum(PSF(:,:,zz)));
    end
   
    PSF = imwarp(PSF,Tform,'OutputView',imref2d(size(tPSF)));
    PSF = imgaussfilt(PSF,SigConv);
    BG = CMOS_Noise + Bg;
    BgPSF = I*PSF+BG;
end

end