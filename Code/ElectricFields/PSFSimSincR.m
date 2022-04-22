% [h,amp3d]=PSFSimSincR(ImageParam,PSFParam,doExtend,AddParams,AddPhase)  : simulates a PSF using the SimLens function and the Sinc(abs(R)) method
% ImageParam is a structure of image parameters (see GenericPSFSim)
% PSFParam: structure of PSF parameters (see GenericPSFSim)
% doExtend: optional flag, to indicate if the calculation is over a padded region that is then cut away. (default=0)
% AddParams: this includes the ri and thicknesses of samples, coverslip and immersion medium. See help GenericPSFSim.m to learn more
%         * if AddParams == [] or not added to the function, the code will run as there is no interfaces and no aberrations due to that, n = PSFParam.n will be used troughout the calculation in this case.
% AddPhase : additional phase that can be added into the function // in 2D
% 
% Example2: spherical aberration due to ri mismatch and Fresnel coefficients addded 
% ImageParam=struct('Sampling',[40 40 80],'Size',[128 128 32]);
% PSFParam=struct('NA',1.4,'n',1.518,'MinOtf',1.2e-3,'lambdaEm',520,'Aplanatic',1,'polarization','linearX');
% InputAmp=Make2DPolPhase(ImageParam,PSFParam);
% AddParams=struct('ns',1.518,'ng',1.518,'ni',1.518,'ng0',1.518,'ni0',1.518,'ts',0,'tg',1.7e5,'tg0',1.7e5,'wd',1.5e5);
% orderlist={[3,1];[1 -1 3*pi 60]}; % coma aberration, tilt by 60 degrees about the ?-axis: check ZernikePoly.m for more info 
% Zernikephse = ZernikePoly(orderlist,ImageParam,PSFParam);
% [h,amp3d]=PSFSimSincR(ImageParam,PSFParam,1,AddParams,Zernikephse);
% 
% Copyright (C) 2021, Rainer Heintzmann,  All rights reserved.  
% Last modified 23.08.21 by R. Dina 
%_________________________________________________________________________
function [h,amp3d]=PSFSimSincR(ImageParam,PSFParam,doExtend,AddParams,AddPhase)

if nargin<5 || isempty(AddPhase)
    AddPhase=0;
end

if nargin<4 || isempty(AddParams)
    AddParams=[];
    ndes=PSFParam.n;
elseif nargin>=4
    if isstruct(AddParams)
        PSFParam.n=AddParams.ni; % the field is propagating in the immersion medium whose refractive index is given in real condition
        ndes=AddParams.ni0;  % ri in the design condition
    else
        error('Please provide the ri and thicknesses parameters in a struct format as the 4th argument of this function or set it empty ([]).\nType help PSFSimSincR or help SimLens on command window the co for more info.')
    end
end

if ~isfield(PSFParam,'Aplanatic')
    PSFParam.Aplanatic = 1;  % 1 means emission the intensity is scaled by cos(theta) and the electric field is scaled by sqrt(cos(theta)), 0: no aplanatic factor, updated today 29.06.21
end

if nargin < 3
    doExtend=1;
end

BorderFraction=0.1;  % not needed any longer
if doExtend
   OrigSize=ImageParam.Size;
   ImageParam.Size=ceil(ImageParam.Size .* (1+BorderFraction));
end

InputPol=Make2DPolPhase(ImageParam,PSFParam);  % generates the polarization dependent on the ImageParam and PSFParam. Specifically using PSFParam.polarization

% Aplanatic for this function is PSFParam.Aplanatic+1 to compensate for the way that the circular symmetry works in the sincR function
FPlane=mSimLens(InputPol,PSFParam,ImageParam,PSFParam.Aplanatic+1,2,AddParams,AddPhase);  % Circular input polarisation, Emission aplanatic factor, jinc-based-Aperture

% %FPlane=SimLens(InputPol,PSFParam,ImageParam,PSFParam.Aplanatic,2);  % this can be changed, if the Fourier-space of sincR is FORCED to be uniform (see Pro3DSincR.m) in the projection
[amp3d,h]=Proj3DSincR(ImageParam,PSFParam,FPlane,BorderFraction);

if doExtend
   ImageParam.Size=OrigSize;
   h=extract(h,ImageParam.Size);
   amp3d=extract(amp3d,ImageParam.Size);
end

h=h/sum(h);

