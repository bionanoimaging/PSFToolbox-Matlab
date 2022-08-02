% [RPlane,I,FPlane]=mPropagateFieldBlock(APlane,PosMidZ,FPlane,PSFParam,ImageParam,LimitByAperture,AddParams) : Propagates an electical field in XY to a 3D data block
% INPUT:
%   APlane : A 2D cell array (or image_arr) with vectorial components of the field. Can also be a single amplitude.
%   PosMidZ : The Z distance of the middle plane to propagate to (in units of the pixel units)
%   ImageParam, PSFParam : Parameters for image and PSF. See "help GenericPSFSim.m" for documentation
%   LimitByAperture:  If true, an Aperture as given by the PSFParam and ImageParam will be applied. If your data already went through an aperture, please set to zero. default=true
%   FPlane : an optional argument. If given this is used instead of the real-space field argument 1 (i.e. APlane). If empty, the FPlane is computed and returned in the first call.
%                 The field has to be provided as an image array with 3 components
% AddParams : contains additional parameters about refractive indexes (ri) and thicknesses of the sample, coverlsip and immersion medium in real and design conditions. 
%                       The emission wavelength is function of the medium in which the light is propagating (i.e. ri of immersion medium in real condition) whereas the 
%                       opening aperture is limited by how the optical system was designed (i.e. function of the ri of immersion medium in design condition).
% 
% OUTPUT: 
%   RPlane : electric field (in 4D) in the image plane
%   I : Intensity in the image plane
%   FPlane : electric field in the back focal plane
%
% EXAMPLE: Simulate a vectorial PSF, by first creating a proper back focal plane using SimLens
%   ImageParam=struct('Sampling',[40 40 100],'Size',[256 256 16]);
%   PSFParam=struct('NA',1.2,'n',1.33,'lambdaEm',520);
%   FPlane=SimLens([],PSFParam,ImageParam,0,1);
%   [RPlane,I,FPlane]=PropagateFieldBlock([],0.0,FPlane,PSFParam,ImageParam,0)  % calculates a focal stack 
% 
% Copyright (C) 2021, Rainer Heintzmann,  All rights reserved.
% This is a modified version of PropagateFieldBlock [edit: 28.06.21 - R. Dina]
% Aberration due to refractive index mismatch is included
% ________________________________________________________________________
function [RPlane,I,Filtered,Propagator]=mPropagateFieldBlock(APlane,PosMidZ,FPlane,PSFParam,ImageParam,LimitByAperture,AddParams)

if nargin < 6
    LimitByAperture=1;
end

if nargin<7 || isempty(AddParams)  % i.e. if the structure is not given
    AddParams=[]; % for our next reference
    ndes=PSFParam.n; % ri in the immersion medium in a design or ideal condition
    nreal=ndes; % ri in the immersion medium in real condition 
    nlist=[ndes, nreal];
else
    ndes=AddParams.ni0;
    nreal=AddParams.ni;
    PSFParam.n=nreal;  % the field is propagating in the immersion medium whose refractive index is given in real condition
    nlist=[AddParams.ni, AddParams.ni0, AddParams.ng, AddParams.ng0, AddParams.ns];
end
% -- end --

NA=min([PSFParam.NA, nlist]); % the aperture is defined by the minimum aperture
PSFParam.NA=NA;
% NA=min(PSFParam.NA/ndes, PSFParam.NA/nreal); % the aperture is defined by the minimum aperture
lambda=PSFParam.lambdaEm / PSFParam.n;


scales=ImageParam.Sampling;

if numel(ImageParam.Sampling) < 3
    if numel(ImageParam.Size) < 3
        ImageParam.Sampling(3)=0;
    else
        error('ImageParam is missing the z-sampling.')
    end
end

if numel(ImageParam.Size) < 3
    ImageParam.Size(3)=1;
end

if ~isempty(APlane)
    if iscell(APlane) || isa(APlane,'dip_image_array')
        Field=cat(4,APlane{:});
    else
        Field=APlane;
    end
end
if nargin < 3 || isempty(FPlane)
    % FPlane=dip_fouriertransform(Field,'forward',[1 1 1 0]);
    Field=expanddim(Field,4);
    FPlane=ft3d(Field);
end
clear Field
if iscell(FPlane) || isa(APlane,'dip_image_array')
        FPlane=cat(4,FPlane{:});
end
if ndims(FPlane) < 3
    FPlane=expanddim(FPlane,4);
end

sz=size(FPlane);
%kxysqr=(xx(sz,'freq')/scales(1)*xx(sz,'freq')/scales(1)+yy(sz,'freq')/scales(2)*yy(sz,'freq')/scales(2));
kxysqr=(abssqr(xx(sz,'freq')/scales(1))+abssqr(yy(sz,'freq')/scales(2)));
kz=2*pi*sqrt(1/(lambda*lambda)-kxysqr);
kz(1/(lambda*lambda)-kxysqr < 0)=0;

% PSFParamAp=PSFParam; PSFParamAp.NA=NA; % define properly the amount of light going through the system
jpsf=jincPSF(ImageParam,PSFParam,ndes);
% jpsf=jincPSF(ImageParam,PSFParam,ndes);
% jpsf=DampEdge(jpsf,0.05,2,1,3);
if (LimitByAperture)
    Aperture= ft(jpsf);  % use a slightly smooth version of the aperture based on its real space representation
    % Aperture= kxysqr*lambda*lambda < NA;
    Filtered=FPlane .* Aperture;
    clear Aperture;  % not needed any longer
else
    Filtered=FPlane;
end
%%
Propagator=exp(-1i*kz*(ImageParam.Sampling(3)*zz([1 1 ImageParam.Size(3)])+PosMidZ)); % constructs the propagator to one particular plane
% the negative sign on the propagator means that the ASF propagates upwards towards the lense. The k-shere points upwards, when doing fft(asf)
% PlaneDamp=DampEdge(Filtered.*Propagator,0.05,3,1,2);
RPlane = ift2d(Filtered.*Propagator); 
I=abssqr(RPlane(:,:,:,0));
for fc=1:size(RPlane,4)-1
    I=I+abssqr(RPlane(:,:,:,fc)); % other vector field components (Y,Z) if existing
end
