% [Propagated,I,ZoomedPupil]=RWFastProp(OrigPupil,PSFParam,ImageParam,ZDist) : propagtes a (4D) pupil to a volume using 1D FFTs along Z
% OrigPupil : original Pupil (not resampled)
% PSFParam : List of PSF parameters
% ImageParam: List of image parameters
% ZDist: Z-Position from the focus to propagate to. If ImageParam.Size(3) is larger than one, a stack is generated with ZDist being the position of the middle (floor(size/2)) slice.
%
% Example:
% 
% ImageParam=struct('Sampling',[40 40 100],'Size',[256 256 16]);
% PSFParam=struct('NA',1.2,'n',1.33,'lambdaEm',520,'polarization','circular');
% FPlane=exp(1i*phiphi(ImageParam.Size(1:2)));  % simulates with smooth aperture
% [Propagated,ZoomedPupil]=RWFastProp(PSFParam,ImageParam,FPlane)
% 
% Copyright (C) 2020, Rainer Heintzmann,  All rights reserved.
%_________________________________________________________________________
function [I,Propagated]=RWFastProp(PSFParam,ImageParam,PupilPlane)
global use_xyz_cuda; 
global use_ramp_cuda; 
saved_use_xyz_cuda=use_xyz_cuda;
saveduse_ramp_cuda=use_ramp_cuda;

if ~isfield(PSFParam,'Aplanatic')
    PSFParam.Aplanatic = -1;  % -1 means emission, 0: no aplanatic factor, 1: excitation
end

lambda=PSFParam.lambdaEm ./ PSFParam.n;

sz = ImageParam.Size(1:2);
ktotal=2*pi/lambda;
if ~isfield(PSFParam,'useCuda') || PSFParam.useCuda == 0
    mykr = ktotal*rrscale(sz,ImageParam.Sampling);
else
    use_xyz_cuda=1.0;use_ramp_cuda=1.0;
    mykr = ktotal*rrscale(sz,ImageParam.Sampling);
end

% sz=sz*2;
% mykr = ktotal*rrscale(sz,ImageParam.Sampling./2); % oversampling twice

Phi = atan2(ImageParam.Sampling(2)*yy(sz),ImageParam.Sampling(1)*xx(sz));
Phi2 = Phi*2.0;
cos2Phi = cos(Phi2);
sin2Phi = sin(Phi2);
cosPhi = cos(Phi);
sinPhi = sin(Phi);

cosAlpha = sqrt(1.0-abssqr(PSFParam.NA/PSFParam.n));
midZ = floor(ImageParam.Size(3)/2);

Xi = 1.0-(zz([1 1 ImageParam.Size(3)],'freq')+0.5)*(1-cosAlpha); % is cos(theta) 
Xi = Xi - (1.0-Xi(2))/4.5;  % seems to be a good compromize to reduce the error

% Xi=permute(expanddim(dip_image(cos(linspace(0,asin(PSFParam.NA/PSFParam.n), 3.*ImageParam.Size(3)))),3),[3 2 1]);
% Xi=Xi(:,:,end:-1:0);
% % Xi=(Xi-min(Xi)).*((1-cosAlpha)./(max(Xi)-min(Xi)))+cosAlpha;
% XiStep=mean(abs(Xi(:,:,1:end)-Xi(:,:,0:end-1)));
% Xi=linspace(cosAlpha,1,ImageParam.Size(3));
% Xi=permute(expanddim(dip_image(Xi),3),[3 2 1]);

sinTheta = sqrt(1-abssqr(Xi));
sqrtXi = sqrt(Xi);

BesselArg = sinTheta*mykr;   % eqn. 3.2, Richards & Wolf II
if nargin > 2
    error('Pupil-Plane modifications. Not implemented yet.')
    sqrtXi = sqrtXi * ft(PupilPlane)/besselj(0,BesselArg);
end
FT_I0 = sqrtXi * (1+Xi) * besselj(0,BesselArg);
FT_I1 = sqrtXi * sinTheta * besselj(1,BesselArg);
FT_I2 = sqrtXi * (1-Xi) * besselj(2,BesselArg);

% scale = 3*lambda / ImageParam.Sampling(3)/ abs(double(squeeze(Xi(end)-Xi(0))));%(1-cosAlpha);%-XiStep);
% scale=lambda / ImageParam.Sampling(3);
scale = lambda / ImageParam.Sampling(3) / (1-cosAlpha);
I0 = czt_z(FT_I0,scale); I1 = czt_z(FT_I1,scale); I2 = czt_z(FT_I2,scale);   % 1D FFTs along Z.

if ~isfield(PSFParam,'polarization') 
    PSFParam.polarization = 'circular';
end
switch PSFParam.polarization
    case 'linearX'
        Propagated = cat(4,I0+I2*cos2Phi, I2*sin2Phi,-2.0*1i*I1*cosPhi);  % eq. 2.30 Richards & Wolf II
    case 'linearY'
        Propagated = cat(4,I2*sin2Phi,I0+I2*cos2Phi,-2.0*1i*I1*sinPhi);  % eq. 2.30 Richards & Wolf II
    case 'circular'
        Propagated = cat(4,I0+I2*(1i*sin2Phi+cos2Phi), 1i*I0+I2*(1i*cos2Phi +sin2Phi),-I1*2.0*(1i*cosPhi-sinPhi));  % eq. 2.30 Richards & Wolf II
    otherwise
        error('Polarization type not implemented for RWFastProp');
end

if 0
    refAmp = SimLens([],PSFParam,ImageParam,0,1);
    myAmp=ft2d(Propagated(:,:,midZ,:));
    Eps =1e-6;
    corFac = refAmp *conj(myAmp) /(abssqr(myAmp) + Eps);
    Propagated = ift2d(ft2d(Propagated) .* corFac);
end

I=squeeze(sum(abssqr(Propagated),[],4));
myNorm = sum(I(:,:,midZ));
I = I/myNorm;
Propagated = Propagated /sqrt(myNorm);

% single slice cannot be normalized as it would skew the 3D PSF
use_xyz_cuda=saved_use_xyz_cuda;
use_ramp_cuda=saveduse_ramp_cuda;
