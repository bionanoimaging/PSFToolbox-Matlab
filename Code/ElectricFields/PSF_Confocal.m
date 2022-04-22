% [h,otf]=PSF_Confocal(ImageParam,PSFParam,Pinholesize) : simulates a confocal PSF based on kSimPSF.
% ImageParam: structure with image parameters. See GenericPSFSim for details
% PSFParam: structure with PSF parameters. See GenericPSFSim for details
% Pinholesize: size of the pinhole in Airy units. Default: 1.0 AU
% 
% Last modified: 17/09/2020
% Rainer Heintzmann
% Copyright (c) 2020
% Nano-imaging Research group, 
% Friedrich Schiller University and  Leibniz Institute of Photonic Technology, Jena, Germany
% 
% EXAMPLE:
% ImageParam=struct('Sampling',[10 10 50],'Size',[128 128 32]);
% PSFParam=struct('NA',1.2,'n',1.33,'lambdaEm',520,'lambdaEx',488);
% [h_conf10, OTF_conf10] = PSF_Confocal(ImageParam,PSFParam,1.0,'SlicePropagation');
%_________________________________________________________________________
%
function [h,otf]=PSF_Confocal(ImageParam,PSFParam,Pinholesize,Method, nocheck)
if nargin < 3 || isempty(Pinholesize)
    if isfield(PSFParam,'PinholeSize')
        Pinholesize=PSFParam.PinholeSize;  
    else
        Pinholesize=1.0;  
    end
end
if nargin < 4
    Method='RichardsWolff';
end
if nargin < 5 
    CheckPSFParams(ImageParam,PSFParam, 0); % make some sanity checks
else
    CheckPSFParams(ImageParam,PSFParam, nocheck); % make some sanity checks
end

if strcmp(Method,'RichardsWolff') || strcmp(Method,'RichardsWolffInt')
    h =kSimPSF({'scaleX',ImageParam.Sampling(1);'scaleY',ImageParam.Sampling(2);'scaleZ',ImageParam.Sampling(3);'sX',ImageParam.Size(1);'sY',ImageParam.Size(2);'sZ',ImageParam.Size(3);'na',PSFParam.NA;'ri',PSFParam.n;'lambdaEm',PSFParam.lambdaEm;'lambdaEx',PSFParam.lambdaEx;'circPol',1;'confocal',1;'pinhole',Pinholesize;'nonorm',1});
else
    PSFParamEx=PSFParam;
    PSFParamEx.Aplanatic=1; % excitation
    PSFParamEx.lambdaEm=PSFParam.lambdaEx; % excitation
    PSFParamEx.Mode = 'Widefield';
    h_ex=squeeze(GenericPSFSim(ImageParam,PSFParamEx,Method,1));
    PSFParamEm=PSFParam;
    PSFParamEm.Aplanatic=-1; % emission
    PSFParamEm.lambdaEx=PSFParam.lambdaEm; % excitation
    PSFParamEm.Mode = 'Widefield';
    h_em=squeeze(GenericPSFSim(ImageParam,PSFParamEm,Method,1));
    otf_em=ft(h_em)*sqrt(prod(ImageParam.Size(1:2)));
    AU = 1.22/2.0 * PSFParam.lambdaEm/PSFParam.NA./ImageParam.Sampling(1:2);
    ftradius=Pinholesize*AU;  % pixels in Fourier space (or real space, if starting in Fourier space)
    myscales=ftradius./ImageParam.Size(1:2); 
    ftpinhole=jinc(ImageParam.Size(1:2),myscales);
    pinholeArea=pi*prod(ftradius);
    ftpinhole=ftpinhole/ftpinhole(MidPosX(ftpinhole),MidPosY(ftpinhole))*pinholeArea
    % pinhole=real(ift(ftpinhole))/sqrt(prod(mysize));
    h_em_pin=ift(otf_em .* ftpinhole)/sqrt(prod(ImageParam.Size(1:2)));  % convolve with the pinhole
    h=h_ex .* h_em_pin;
end

otf = ft(h) * sqrt(prod(size(h)));
