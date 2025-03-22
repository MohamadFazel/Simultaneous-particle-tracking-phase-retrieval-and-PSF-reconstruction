function [Data,PSFstruct,Data1]=genDynamicParticle(PSFstruct,PSFtype,D,Dt,NFrame,I,Bg,Scale)

if nargin < 8
    Scale = 1;
elseif nargin < 7
    error('There must be at least 7 inputs'); 
end
if ~isfield(PSFstruct,'SubPixel')
    PSFstruct.SubPixel = 4;
end

ImSZ = PSFstruct.NPix * PSFstruct.PixelSize;
X0 = ImSZ*rand()/3 - ImSZ/6;
Y0 = ImSZ*rand()/3 - ImSZ/6;
Z0 = -300 + 600*rand();
if ~isfield(PSFstruct,'Z') || isempty(PSFstruct.Z)
    PSFstruct.Z = Z0+cumsum(sqrt(2*D*Dt)*randn(NFrame,1)); %Z-samples (nm)
end
if ~isfield(PSFstruct,'X') || isempty(PSFstruct.X)
    PSFstruct.X = X0+cumsum(sqrt(2*D*Dt)*randn(NFrame,1)); %X-samples (nm)
end
if ~isfield(PSFstruct,'Y') || isempty(PSFstruct.Y)
    PSFstruct.Y = Y0+cumsum(sqrt(2*D*Dt)*randn(NFrame,1)); %Y-samples (nm)
end

[PSF,PSFstruct] = genPSF(PSFstruct,PSFtype,Scale);
Data1 = poissrnd(I*PSF+Bg(1));

DelX = PSFstruct.DelX;
Data = zeros([size(Data1),size(DelX,1)+1]);
Data(:,:,:,1) = Data1;
PSFtype = 'Given';
tPSFstruct = PSFstruct;

% optimizer = registration.optimizer.RegularStepGradientDescent;
% metric = registration.metric.MeanSquares;
for ii = 1:size(DelX,1)
    tPSFstruct.Z = PSFstruct.Z+DelX(ii,3); %Z-samples channel2 (nm)
    tPSFstruct.X = PSFstruct.X+DelX(ii,1); %Z-samples channel2 (nm)
    tPSFstruct.Y = PSFstruct.Y+DelX(ii,2); %Z-samples channel2 (nm)

    % PSFtype.Mask = PSFstruct.Pupil_Mag;
    % PSFtype.Phase = PSFstruct.Pupil_Phase;
    [tPSF,~] = genPSF(tPSFstruct,PSFtype,Scale);
    if isfield(PSFstruct,'Tform')
        tPSF = imwarp(tPSF,PSFstruct.Tform(ii+1),'OutputView',imref2d(size(tPSF)));
    end
    
    Data2 = poissrnd(I*tPSF+Bg(ii+1));
    Data(:,:,:,ii+1) = Data2;
end

end