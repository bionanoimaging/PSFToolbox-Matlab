% weight=IterateCoefficients(subpix,CutoffK,TotalSize,RollOff,Iterations) : Calculates a matrix of coefficients optimized for interpolation in Fourier space. To this aim an IFTA algorithm is used.
% subpix : array of sizes for number of subdivisions between pixels
% CutoffK : number of pixels to use (in both directions) as a kernel
% TotalSize : Size of the multidimensional FFT array
% RollOff : Sizes to use as border pixels
% Iterations : Number of iterations to use during optimization
%
% Example:
% map=IterateCoefficients(10,10,256,20,1000);  % Generates an interpolation map in 0.1 pixel distances for one fft direction
% map=IterateCoefficients([10 10],[10 10],[256 256],[10 10],100);  % 2D interpolation kernels for a 2D FFT with only 10 pixels rolloff everywhere
%
% Copyright (C) 2020, Rainer Heintzmann,  All rights reserved.
%_________________________________________________________________________
%
function weight=IterateCoefficients(subpix,CutoffK,TotalSize,RollOff,Iterations)
if nargin < 5
    Iterations=1000;
end
if TotalSize<=1
    fprintf('Warning z-size is one. Replacing iteration matrix with 1\n');
    weight=1;
    return;
end
MySize=[TotalSize subpix];
if length(TotalSize)==1
    allramps=xx(MySize,'freq') .* ramp(MySize,2,'corner')/(subpix(1)-1);
    toprocess=[1 0];
    fmask=~ (abs(xx(MySize))<=CutoffK(1));
    rmask=(abs(xx(MySize))<(floor(TotalSize(1)/2)-RollOff));
elseif length(TotalSize)==2
    allramps=xx(MySize,'freq') .* ramp(MySize,3,'corner')/(subpix(1)-1)+yy(MySize,'freq') .* ramp(MySize,4,'corner')/(subpix(2)-1);
    toprocess=[1 1 0 0];
    fmask=~((abs(xx(MySize))<=CutoffK(1)) & (abs(yy(MySize))<=CutoffK(2)));
    rmask=(abs(xx(MySize))<(floor(TotalSize(1)/2)-RollOff(1))) & (abs(yy(MySize))<(floor(TotalSize(2)/2)-RollOff(2)));
elseif length(TotalSize)==3
    allramps=xx(MySize,'freq') .* ramp(MySize,5,'corner')/(subpix(1)-1)+yy(MySize,'freq') .* ramp(MySize,5,'corner')/(subpix(2)-1)+zz(TotalSize,'freq') .* ramp(MySize,6,'corner')/(subpix(3)-1);
    toprocess=[1 1 1 0 0 0];
    fmask=~((abs(xx(MySize))<=CutoffK(1)) & (abs(yy(MySize))<=CutoffK(2)) & (abs(zz(MySize))<=CutoffK(3)));
    rmask=(abs(xx(MySize))<(floor(TotalSize(1)/2)-RollOff(1))) & (abs(yy(MySize))<(floor(TotalSize(2)/2)-RollOff(2))) & (abs(zz(MySize))<(floor(TotalSize(3)/2)-RollOff(3)));
else
    error('IterateCoefficients only for a maximum of 3 dimensions');
end
perfsignal=exp(2*pi*i*allramps)/sqrt(size(allramps,1));  % Dim it down to account for Fourier space 

rspace=perfsignal;
meanPrevErr=0;
for n=1:Iterations
    fspace=dip_fouriertransform(rspace,'forward', toprocess);
    fspace(fmask)=0;
    fspace=real(fspace);  % also force these coefficients to be real valued
%    dipshow(1,fspace);

    rspace=dip_fouriertransform(fspace,'inverse', toprocess);
    maxErr=max(abs(rspace(rmask)-perfsignal(rmask)));
    meanErr=mean(abs(rspace(rmask)-perfsignal(rmask)));
    fprintf('Error iteration %d/%d is  : max %g mean %g\n',n,Iterations,maxErr,meanErr);
    if abs((meanErr-meanPrevErr)/meanErr) < 1e-6
        fprintf('no more update\n');
        break;
    end
    meanPrevErr=meanErr;
%    dipshow(2,real(rspace));
    rspace(rmask)=perfsignal(rmask);
%    drawnow;
end
rspace=dip_fouriertransform(fspace,'inverse', toprocess);
fprintf('Error %g is  : max %g mean %g\n',subpix,max(abs(rspace(rmask)-perfsignal(rmask))),mean(abs(rspace(rmask)-perfsignal(rmask))));

% if 1
%     pos = 0.25
%     myposy = floor(size(rspace,2)*pos);
%     toPlot = imag(rspace(:,myposy));
%     toPlot = cat(1,toPlot,toPlot);
%     plot(toPlot);
%     hold on
%     toPlot = imag(perfsignal(:,myposy));
%     toPlot = cat(1,toPlot,toPlot);
%     plot(toPlot,'r');
% end

if length(TotalSize)==1
    weight=fspace(floor(TotalSize(1)/2) - CutoffK(1):floor(TotalSize(1)/2) + CutoffK(1),:);
elseif length(TotalSize)==2
    weight=fspace(floor(TotalSize(1)/2) - CutoffK(1):floor(TotalSize(1)/2) + CutoffK(1),floor(TotalSize(2)/2) - CutoffK(2):floor(TotalSize(2)/2) + CutoffK(2),:,:);
elseif length(TotalSize)==3
    weight=fspace(floor(TotalSize(1)/2) - CutoffK(1):floor(TotalSize(1)/2) + CutoffK(1),floor(TotalSize(2)/2) - CutoffK(2):floor(TotalSize(2)/2) + CutoffK(2),floor(TotalSize(3)/2) - CutoffK(3):floor(TotalSize(3)/2) + CutoffK(3),:,:,:);
end
normfac=repmat(sqrt(sum(abssqr(weight),[],1)),[size(weight,1) 1]);
weight=weight ./ normfac;  % this normalization is needed for energy conservation in projections
weight=real(weight);  % again make sure that the datatype is real.