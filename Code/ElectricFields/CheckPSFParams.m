% success=CheckPSFParams(ImageParams,PSFParams,optionalIgnore)
% This function checks whether all the parameters are in order and the Nyquist sampling is met or not. 
% success=1 means all is ok.
%  success=1 does not necessarily mean than that the Nyquist sampling is met at the given arguments if optionalIgnore=1. It could mean that the data is undersampled but still carry on with the calculation
% 
% EXAMPLE:
% ImageParam=struct('Sampling',[500 200 250],'Size',[128 128 32]);
% PSFParam=struct('NA',1.2,'n',1.33,'lambdaEm',520);
% success=CheckPSFParams(ImageParam,PSFParam,0)
% 
% Copyright (C) 2020, Rainer Heintzmann,  All rights reserved.
%_________________________________________________________________________
function success=CheckPSFParams(ImageParam,PSFParam,optionalIgnore)
success = 1;
if nargin < 3
    optionalIgnore = 0;
end
if PSFParam.NA>PSFParam.n
    if optionalIgnore
        msgbox('Error: NA is bigger than refractive index.')
    end
    error('Error while checking PSF. The NA has to be lower than the refractive index! Change PSFParam.n and/or PSFParam.NA.')
end

[IncoherentRadius,PupilRadius,OTF_KZExtent]=AbbeMaxRadiusFromPSF(PSFParam,ImageParam);

if isfield(PSFParam,'Mode') && strcmp(PSFParam.Mode,'Confocal') && (isfield(PSFParam,'PinholeSize')) && (PSFParam.PinholeSize < 2.0)
    IncoherentRadius = IncoherentRadius * (1+PSFParam.lambdaEm/PSFParam.lambdaEx);    
end
if isfield(PSFParam,'Mode') && strcmp(PSFParam.Mode,'2-Photon') 
    IncoherentRadius = IncoherentRadius * 2;
end
if isfield(PSFParam,'Mode') && strcmp(PSFParam.Mode,'Lightsheet') 
    PupilRadius = ImageParam.Size(3) .* (2*PSFparam.NA/PSFparam.lambdaEx).*ImageParam.Sampling(3)/2.0; %aperture radius (along X)
    IncoherentRadius = 2*PupilRadius;
    OTF_KZExtent = OTF_KZExtent + IncoherentRadius;
end

relSamplingX = IncoherentRadius(1)/floor(ImageParam.Size(1)/2);
relSamplingY = IncoherentRadius(2)/floor(ImageParam.Size(2)/2);
relSamplingZ = OTF_KZExtent/floor(ImageParam.Size(3)/2);

if relSamplingX > 1.0 || relSamplingY > 1.0 || relSamplingZ > 1.0
    if isfield(PSFParam,'Mode') 
        txt = sprintf('%s PSF Caculation:\n',PSFParam.Mode);
    else
        txt = '';
    end
    if relSamplingX > 1.0
        txt = sprintf('%sSampling X: %g µm, by %2.2g too large!\n',txt,ImageParam.Sampling(1),relSamplingX);
    end
    if relSamplingY > 1.0
        txt = sprintf('%sSampling Y: %g µm, by %2.2g too large!\n',txt,ImageParam.Sampling(2),relSamplingY);
    end
    if relSamplingZ > 1.0
        txt = sprintf('%sSampling Z: %g µm, by %2.2g too large!\n',txt,ImageParam.Sampling(3),relSamplingZ);
    end
    if optionalIgnore > 1
        txt = sprintf('%sUNDERSAMPLED! Calculate anyway?\n',txt);
        answer = questdlg(txt);
        if strcmp(answer,'Yes') ~= 1
            success = 0;
            return;
            % error('Error while checking PSF. The data is undersampled along X! reduce the sample spacing or reduce the NA of the objective.')
        end
    else
        fprintf('%s',txt);
        if optionalIgnore == 0
            error('Error while checking PSF. The data is undersampled! reduce the sample spacing or reduce the NA of the objective.')
        end
    end
end


