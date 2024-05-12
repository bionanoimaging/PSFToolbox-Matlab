% jincPSF(ImageParam,PSFParam)  : calculates a 2D jinc scalar PSF for low NA
% 
% Example:
% ImageParam=struct('Sampling',[20 20 50],'Size',[128 128 32]);
% PSFParam=struct('NA',1.2,'n',1.33,'lambdaEm',520,'Aplanatic',-1);  
% jincPSF(ImageParam,PSFParam)
%
% Copyright (C) 2020, Rainer Heintzmann,  All rights reserved.
% edit 29.06.21 R. Dina
%_________________________________________________________________________

function res=jincPSF(ImageParam,PSFParam,ndes)
% lambda=PSFParam.lambdaEm./ImageParam.Sampling(1:2);
% na=PSFParam.NA;
if nargin<3
    ndes=PSFParam.n; % refractive index of the immersion medium in which the objective lens is designed for. The maximum opening aperture which defined the amount of light that can be collected by the lens is defined of this 
end
% PSFParam.n: ri in which the light is propagating
% ndes: ri in which the light is designed to propagate
% PSFParam.n can be different from ndes
lambda=(PSFParam.lambdaEm./PSFParam.n)./ImageParam.Sampling(1:2); % the wavelength is defined in the medium where the light is propagating
na=min(PSFParam.NA./ndes, PSFParam.NA./PSFParam.n); % the opening aperture is defined by the minimum of how it was designed by the manufacturer or the available given by the imaging condition
AbbeLimit=lambda/na;   % coherent Abbe limit, central illumination, not incoherent
ftradius=ImageParam.Size(1:2)./AbbeLimit;
myscales=ftradius./ImageParam.Size(1:2); % [100 100]/(488/0. 3);
res=jinc(ImageParam.Size(1:2),myscales);  % Airy disc
