function Out=blockresample(In,Zm)
%blockresample() changes the pixel size of the input dataset.
%   blockresample(), take an image and based on the value of the input Zm,
%   changes the image. For the values larger than one, it changes a pixel
%   into a block of identical pixels. For values smaller than one, it take  
%   a block of pixels and makes them a single pixel, where the value of 
%   that pixel is the sum of all the pixel values inside that block.
%
% INPUT:
%   In:   The input image. (Default = No default)
%   Zm:   The zoom-factor, it can be any positive value. (Pixels)
%
% OUTPUT:
%   Out:  The output image.
%
% REQUIRES: 
%   MATLAB 2014a or later versions.
%
% CITATION:
%   Mohamadreza Fazel, Michael J. Wester, Hanieh Mazloom-Farsibaf,
%   Marjolein M.B.M. Meddens and Keith A. Lidek, "Bayesian Multiple Emitter
%   Fitting using Reversible Jump Markov Chain Monte Carlo".
%
% Created by:
%   Mohamadreza Fazel (Presse Lab 2022)
%

SZIn = size(In);
if numel(SZIn) < 3
    SZIn(3) = 1;
end
Out = zeros(SZIn(1)*Zm,SZIn(2)*Zm,SZIn(3));
if Zm > 1
    pixelsX=SZIn(1);
    pixelsY=SZIn(2);
    In=squeeze(In);
    for pp = 1:SZIn(3)
        out=repmat(reshape(In(:,:,pp),[1 pixelsX*pixelsY]),[Zm 1]);
        out=reshape(out,[pixelsX*Zm pixelsY]);
        out=(shiftdim(out,1));
        out=repmat(reshape(out,[1 pixelsX*pixelsY*Zm]),[Zm 1]);
        out=reshape(out,[pixelsY*Zm pixelsX*Zm]);
        Out(:,:,pp)=(shiftdim(out,1));
    end
elseif Zm < 1
    if floor(1/Zm) < Zm
       error('Sorry! 1/ZoomFactor must be an integer number.'); 
    end
    if SZIn(1)~=SZIn(2)
        error('Sorry! The input image must be a squre image.');
    end
    if floor(SZIn(1)*Zm)~=SZIn(1)*Zm
        error('Sorry! The size of the image must be proportional to the zoom-factor');
    end
  %for pp = 1:SZIn(3)
    %Image = In(:,:,pp);
    Image = In;
    A = reshape(Image,[1/Zm SZIn(1)*Zm 1/Zm SZIn(2)*Zm]);
    B = permute(A,[2 4 1 3]);
    C = reshape(B,[SZIn(1)*Zm SZIn(2)*Zm (1/Zm)^2]);
    %Out(:,:,pp)=sum(C,3);
    Out = sum(C,3);
  %end  
else 
    Out = In;
end
end
