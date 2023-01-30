% [h,amp3d]=SlicePSFSim(ImageParam,PSFParam,AddParams,AddPhase) : simulates a PSF using the SimLens function and the slice propagation method
% ImageParam : structure of image parameters (see GenericPSFSim)
% PSFParam: structure of PSF parameters (see GenericPSFSim)
% AddParams : contains additional parameters about refractive indexes (ri) and thicknesses of the sample, coverlsip and immersion medium in real and design conditions. 
%                       The emission wavelength is function of the medium in which the light is propagating (i.e. ri of immersion medium in real condition) whereas the 
%                       opening aperture is limited by how the optical system was designed (i.e. function of the ri of immersion medium in design condition).
%                       e.g. AddParams=struct('ns',1.33,'ng',1.518,'ni',1.516,'ts',2e3,'tg',1.7e5,'tg0',1.7e5,'wd',1.5e5); 
%                         'ns': ri of sample, can be in a vector format if the sample is in a stratified medium.  First element of ns (if a vector) is the ri of where the emitter is at, it is the ri of the medium which the furthest from the coverslip. The last element of the vector is the ri of the medium closest to the coverslip
%                         'ng': ri of coverslip in real condition
%                         'ni': ri of immersion medium in real condition
%                         'ni0': ri of immersion medium in design condition
%                         'ts': thickness of the sample or position where the emitter is at, this should be with the same length as the vector of ns and within the same order as ns.
%                         'tg': thickness of coverslip in real condition
%                         'tg0': thickness of coverslip in design condition
%                         'wd': working distance which is supposed to be the thickness of the immersion medium in design condition
% AddPhase : additional phase that can be added into the function // in 2D
% 
% Example:
% ImageParam=struct('Sampling',[100 100 100],'Size',[256 256 16]);
% PSFParam=struct('NA',0.6,'n',1.33,'lambdaEm',500);
% AddParams=struct('ns',1.35,'ng',1.518,'ni',1.33,'ng0',1.518,'ni0',1.33,'ts',1e3,'tg',1.7e5,'tg0',1.7e5,'wd',1.5e5);
% orderlist={[3,1];[1 -1 1 60]}; % order of the Zernike phase: coma and tilt aberration
% [TiltComaPhase,rho,phi] = ZernikePoly(orderlist,ImageParam,PSFParam);
% [h,amp3d]=SlicePSFSim(ImageParam,PSFParam,AddParams,TiltComaPhase);
% 
% Copyright (C) 2021, Rainer Heintzmann,  All rights reserved.
% edit 29.06.21: add of ri and thicknesses parameters 
%_________________________________________________________________________
function [h,amp3d]=SlicePSFSim(ImageParam,PSFParam,AddParams,AddPhase)
if nargin<4 || isempty(AddPhase)
    AddPhase=0;
end
if iscell(AddPhase)
    AddPhase=ZernikePoly(AddPhase,ImageParam,PSFParam);
end
if nargin<3 || isempty(AddParams)
    AddParams=[];
end
if ~isfield(PSFParam,'Aplanatic')
    PSFParam.Aplanatic = 1;  % 1 means emission the intensity is scaled by cos(theta) and the electric field is scaled by sqrt(cos(theta)), 0: no aplanatic factor, updated today 29.06.21
end
if ~isfield(PSFParam, 'PosMidZ')
    PSFParam.PosMidZ = 0.0;
end
% mSimLens(BFPAperture,lambda,scales,NA,aplanar,smoothAperture,AddParams)

InputPol=Make2DPolPhase(ImageParam,PSFParam);  % generates the polarization dependent on the ImageParam and PSFParam. Specifically using PSFParam.polarization
FPlane=mSimLens(InputPol,PSFParam,ImageParam,PSFParam.Aplanatic,1,AddParams,AddPhase);  % Circular input polarisation, Emission aplanatic factor, jinc-based-Aperture
[amp3d,h]=mPropagateFieldBlock([],PSFParam.PosMidZ,FPlane,PSFParam,ImageParam,0,AddParams); % constructs the vectorial PSF and its intensity distribution. Does NOT apply the aperture a second time
h=h/sum(h); 