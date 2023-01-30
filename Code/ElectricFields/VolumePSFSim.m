% [h,amp3d]=VolumePSFSim(ImageParam,PSFParam,InputAmp,doExtend)  : simulates a PSF using the SimLens function and the slice propagation method
% ImageParam is a structure containing the following minimum entries:
% ImageParam=struct('Sampling',[40 40 100];'Size',[256 256 16]);
% PSFParam is a structure containing the following minimum entries:
% PSFParam=struct('NA',1.2,'n',1.33,'lambdaEm',520);
% InputAmp : option argument with the input polarisation as amplitude for Ex and Ey. If omitted: circular polarisation is assumed
%            The InputAmp is a 2D complex field with two components along Z for X and Y input polarisation respectively
% doExtend: otpional flag, to indicate if the calculation is over a padded region that will be cut away (corresponds to value 1) or not padded region (default=0)
% AddParams: this includes the ri and thicknesses of samples, coverslip and immersion medium. See help GenericPSFSim.m to learn more
%         * if AddParams == [] or not added to the function, the code will run as there is no interfaces and no aberrations due to that, n = PSFParam.n will be used troughout the calculation in this case.
% AddPhase : additional phase that can be added into the function // in 2D
% 
% Example1: non-aberrated
% ImageParam=struct('Sampling',[40 40 80],'Size',[128 128 32]);
% PSFParam=struct('NA',1.4,'n',1.518,'MinOtf',1.2e-3,'lambdaEm',520,'Aplanatic',1,'polarization','linearX');
% InputAmp=Make2DPolPhase(ImageParam,PSFParam);
% [h,amp3d]=VolumePSFSim(ImageParam,PSFParam,InputAmp,1);
% 
% Example2: spherical aberration due to ri mismatch and Fresnel coefficients addded 
% ImageParam=struct('Sampling',[40 40 80],'Size',[128 128 32]);
% PSFParam=struct('NA',1.4,'n',1.518,'MinOtf',1.2e-3,'lambdaEm',520,'Aplanatic',1,'polarization','linearX');
% InputAmp=Make2DPolPhase(ImageParam,PSFParam);
% AddParams=struct('ns',1.518,'ng',1.518,'ni',1.518,'ng0',1.518,'ni0',1.518,'ts',0,'tg',1.7e5,'tg0',1.7e5,'wd',1.5e5);
% orderlist={[3,1];[1 -1 3*pi 60]}; % coma aberration, tilt by 60 degrees about the ?-axis: check ZernikePoly.m for more info 
% Zernikephse = ZernikePoly(orderlist,ImageParam,PSFParam);
% [h,amp3d]=VolumePSFSim(ImageParam,PSFParam,InputAmp,1,AddParams,Zernikephse);
% 
% Copyright (C) 2021, Rainer Heintzmann,  All rights reserved.  
% Last modified 23.08.21 by R. Dina 
%_________________________________________________________________________
function [h,amp3d]=VolumePSFSim(ImageParam,PSFParam,InputAmp,doExtend,AddParams,AddPhase)

if nargin<6 || isempty(AddPhase)
    AddPhase=0;
end
if iscell(AddPhase)
    AddPhase=ZernikePoly(AddPhase,ImageParam,PSFParam);
end

% if nargin<5 || isempty(AddParams)
%     AddParams=[];
%     ndes=PSFParam.n;
if nargin>=5 && isstruct(AddParams)
        PSFParam.n=AddParams.ni; % the field is propagating in the immersion medium whose refractive index is given in real condition
        ndes=AddParams.ni0;  % ri in the design condition
else
     AddParams=[];
    ndes=PSFParam.n;
%         error('Please provide the ri and thicknesses parameters in a struct format as the 5th argument of this function or set it empty ([]).\nType help VolumePSFSim on command window the co for more info.')
end

if ~isfield(PSFParam,'Aplanatic')
    PSFParam.Aplanatic = 1;  % 1 means emission the intensity is scaled by cos(theta) and the electric field is scaled by sqrt(cos(theta)), 0: no aplanatic factor, updated today 29.06.21
end

if nargin < 4
    doExtend=0;
end

BorderFrac=0.2;  % 10% of the z-extend is reserved for boarder.
kernelSize=8;  % +/- 8 pixels in Fourier space for the interpolation kernel

global PupilInterpolators;  % contains the interpolation coefficients. If this does not exist, it was generated in Pupil3DPrepare
PupilInterpolators.Newlambda=PSFParam.lambdaEm;
PupilInterpolators.Newpixelsize=ImageParam.Sampling;
PupilInterpolators.NewNA=PSFParam.NA;  % Put a somewhat larger NA, since the jinc below will fix this.
PupilInterpolators.NewRI=PSFParam.n;
PupilInterpolators.ndes=ndes;

mysize=ImageParam.Size;
if (doExtend)
    mysize(3)=ceil(mysize(3)*(1.0+BorderFrac));
end
Pupil3DPrepare(mysize,BorderFrac,kernelSize);  % Fills the global PupilInterpolators structure


if nargin < 3 || isempty(InputAmp)
    InputAmp=Make2DPolPhase(ImageParam,PSFParam);  % generates the polarization dependent on the ImageParam and PSFParam. Specifically using PSFParam.polarization
end

FPlane=mSimLens(InputAmp,PSFParam,ImageParam,PSFParam.Aplanatic,1,AddParams,AddPhase);  % Circular input polarisation, Emission aplanatic factor, jinc-based-Aperture
atf0=newim(mysize,'scomplex');
atfX=FillProjSphere(atf0,squeeze(FPlane(:,:,:,0)),PupilInterpolators.indexList2D,PupilInterpolators.fullIndex3D,PupilInterpolators.factorList);
atfY=FillProjSphere(atf0,squeeze(FPlane(:,:,:,1)),PupilInterpolators.indexList2D,PupilInterpolators.fullIndex3D,PupilInterpolators.factorList);
atfZ=FillProjSphere(atf0,squeeze(FPlane(:,:,:,2)),PupilInterpolators.indexList2D,PupilInterpolators.fullIndex3D,PupilInterpolators.factorList);
amp3d=cat(4,ift(atfX),ift(atfY),ift(atfZ));  % changed for compatibility to dimension 4. RH 26.11.2018
if (doExtend)
    amp3d=extract(amp3d,ImageParam.Size);  % remove extended region
end
h=squeeze(sum(abssqr(amp3d),[],4));

h=h/sum(h); 