% [h,amp3d]=mPSFSimSincR(ImageParam,PSFParam,doExtend,AddParams,AddPhase)  : simulates a PSF using the SimLens function and the Sinc(abs(R)) method
% A modified version of PSFSimSincR to avoid wrap around
% ImageParam is a structure of image parameters (see GenericPSFSim)
% PSFParam: structure of PSF parameters (see GenericPSFSim)
% doExtend: optional flag, to indicate if the calculation is over a padded region that is then cut away. (default=0)
% AddParams: additional parameters in a structure format having as fields:
%                 'RIReal' (ri of [sample, coverslip, immersion medium] in this order in real condition); 
%                 'RIDesign' (ri of [coverslip, immersion medium] in this order in design condition); 
%                 'ThicknessReal' (thicknesses in real condition with the same order as 'RIReal';
%                 'ThicknessDesign' (thicknesses in design condition with the same order as 'RIDesign';
%         * RIReal will be assumed the same as RIDesign and ThicknessReal same as ThicknessDesign if missing in the structure 
%         * if AddParams == [] or not added to the function, the code will run as there is no interfaces and no aberrations due to that, n = PSFParam.n will be used troughout the calculation in this case.
% AddPhase : additional phase that can be added into the function // in 2D
% 
% Example2: spherical aberration due to ri mismatch and Fresnel coefficients addded 
% ImageParam=struct('Sampling',[40 40 80],'Size',[128 128 32]);
% PSFParam=struct('NA',1.4,'n',1.518,'MinOtf',1.2e-3,'lambdaEm',520,'Aplanatic',1,'polarization','linearX');
% InputAmp=Make2DPolPhase(ImageParam,PSFParam);
% AddParams=struct('RIReal',[1.33 1.518 1.516],'RIDesign',[1.518 PSFParam.n],'ThicknessReal',[2e3 1.7e5 1.5e5],'ThicknessDesign',[1.7e5 1.5e5]);
% orderlist={[3,1];[1 -1 3*pi 60]}; % coma aberration, tilt by 60 degrees about the ?-axis: check ZernikePoly.m for more info 
% Zernikephse = ZernikePoly(orderlist,ImageParam,PSFParam);
% [h,amp3d]=PSFSimSincR(ImageParam,PSFParam,1,AddParams,Zernikephse);
% 
% Copyright (C) 2020, Rainer Heintzmann,  All rights reserved.  
% Edit 19.08.21 -- R. Dina 
%_________________________________________________________________________
function [h,amp3d]=mPSFSimSincR(ImageParam,PSFParam,doExtend,AddParams,AddPhase)

if nargin<5 || isempty(AddPhase)
    AddPhase=0;
end

% if nargin<4 || isempty(AddParams)
%     AddParams=[];
%     ndes=PSFParam.n;
if nargin>=4 && isstruct(AddParams)
        PSFParam.n=AddParams.ni; % the field is propagating in the immersion medium whose refractive index is given in real condition
%         ndes=AddParams.ni0;  % ri in the design condition
else
    AddParams=[];
%         error('Please provide the ri and thicknesses parameters in a struct format as the 4th argument of this function or set it empty ([]).\nType help PSFSimSincR or help SimLens on command window the co for more info.')

end

if ~isfield(PSFParam,'Aplanatic')
    PSFParam.Aplanatic = -1;  % 1 means emission the intensity is scaled by cos(theta) and the electric field is scaled by sqrt(cos(theta)), 0: no aplanatic factor, updated today 29.06.21
end

if nargin < 3
    doExtend=1;
end

BorderFraction=0.25;  
if doExtend
   OrigSize=ImageParam.Size;
   ImageParam.Size=ceil(ImageParam.Size .* (1+BorderFraction));
end

% Expand the window to avoid FT wrap-around // Edit 03.08.21 
expw=0;
if expw==1 % expand window
    BigImageParam=ImageParam;
    bigsz=2.*ImageParam.Size(1:2);
    BigImageParam.Size=[bigsz ImageParam.Size(3)+4];
    InputPol=Make2DPolPhase(BigImageParam,PSFParam);  % generates the polarization dependent on the ImageParam and PSFParam. Specifically using PSFParam.polarization
else 
    BigImageParam=ImageParam;
    bigsz=ImageParam.Size(1:2);
    BigImageParam.Size=[bigsz ImageParam.Size(3)];
    InputPol=Make2DPolPhase(BigImageParam,PSFParam);  
end
% -- change start -- % Question: why the aplanatic for this function is PSFParam.Aplanatic+1 though?

if iscell(AddPhase)
    orderlist=AddPhase; clear AddPhase
    AddPhase=ZernikePoly(orderlist,BigImageParam,PSFParam,0);
else
    if AddPhase~=0 
        if all(size2d(AddPhase)~=bigsz)
            AddPhase=real(ift(extract(ft(AddPhase),bigsz)));
            AddPhase=AddPhase./max(AddPhase);
        end
    end
end
% % FPlane=SimLens(BFPAperture,PSFParam,ImageParam,aplanar,smoothAperture)  
% FPlane=mSimLens(InputPol,PSFParam,BigImageParam,PSFParam.Aplanatic+1,2,AddParams,AddPhase);  % Circular input polarisation, Emission aplanatic factor, jinc-based-Aperture
FPlane=mSimLens(InputPol,PSFParam,BigImageParam,PSFParam.Aplanatic,2,AddParams,AddPhase);  % Circular input polarisation, Emission aplanatic factor, jinc-based-Aperture
%
% % FPlane=SimLens(InputPol,PSFParam,ImageParam,PSFParam.Aplanatic+2,2);  % compensate for the way that the circular symmetry works in the sincR function
% -- end change --

% %FPlane=SimLens(InputPol,PSFParam,ImageParam,PSFParam.Aplanatic,2);  % this can be changed, if the Fourier-space of sincR is FORCED to be uniform (see Pro3DSincR.m) in the projection
[amp3d,h]=Proj3DSincR(BigImageParam,PSFParam,FPlane,BorderFraction);

if doExtend
   ImageParam.Size=OrigSize;
   h=extract(h,ImageParam.Size);
   amp3d=extract(amp3d,ImageParam.Size);
end

% extract the correct size // Edit 03.08.21
if expw==1
    szamp=size(amp3d);
    szamp(1:3)=ImageParam.Size(1:3);
    amp3d=extract(amp3d,szamp);

    szh=size(h);
    szh(1:3)=ImageParam.Size(1:3);
    h=extract(h,szh);
end
% end edit ---
h=h/sum(h);

