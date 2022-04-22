% [IncoherentRadius,PupilRadius,OTF_KZExtent,EwaltRadii]=AbbeMaxRadiusFromPSF(PSFparam,ImageParam,UseExcitation): Calculates the maximum relative radius (0 to 0.5) of the Abbe limit
% returns maximum aperture radius and maximum incoherence aperture radius
% to get relative radius, devide by floor(sizeX);
% Note that the function returns two radii, one for X and one for Y
% 
% Example:
% ImageParam=struct('Sampling',[20 20 50],'Size',[128 128 32]);
% PSFParam=struct('NA',1.2,'n',1.33,'lambdaEm',520,'Aplanatic',-1);  % detection aplanatic factor
% [IncoherentRadius,PupilRadius,EwaltRadii]=AbbeMaxRadiusFromPSF(PSFParam,ImageParam)
% 
% Copyright (C) 2020, Rainer Heintzmann,  All rights reserved.
% edit 30.06.21 R. Dina
%_________________________________________________________________________
%
function [IncoherentRadius,PupilRadius,OTF_KZExtent,EwaltRadii]=AbbeMaxRadiusFromPSF(PSFParam,ImageParam,UseExcitation,AddParams)
if nargin < 3
    UseExcitation=0;
end

if ~UseExcitation
    lambda = PSFParam.lambdaEm;
else
    lambda = PSFParam.lambdaEx;
end
% change --------
if nargin<4
    ni0=PSFParam.n; % refractive index of the immersion medium in which the objective lens is designed for. The maximum opening aperture which defined the amount of light that can be collected by the lens is defined of this 
    ni=ni0;
else
    if isstruct(AddParams)
        ni0=AddParams.ni0; 
        ni=AddParams.ni;
    else % this is to review // we need to make it uniform
        ni0=AddParams;
        ni=PSFParam.n;
    end
end

lambda=lambda./ni; % the wavelength is defined in the medium where the light is propagating
NA=PSFParam.NA./max([ni0 ni]); % the maximum opening aperture is defined by how it was designed by the manufacturer

PupilRadius = ImageParam.Size(1:2) .* (2*NA/lambda).*ImageParam.Sampling(1:2)/2.0; %aperture radius (along X)
IncoherentRadius = 2*PupilRadius;
 
if nargout > 2
    alpha=asin(NA);  % half opening angle
    EwaltRadii=ImageParam.Size .* (2/lambda).*ImageParam.Sampling/2.0; %aperture radius (along X)
    OTF_KZExtent = EwaltRadii(3) .* (1-cos(alpha));
%    OTF_KZExtent = ImageParam.Size(3) .* (2*PSFparam.NA/PSFparam.lambdaEm).*ImageParam.Sampling(3)/2.0 .* (1-cos(alpha))/sina;
end

% end ------
