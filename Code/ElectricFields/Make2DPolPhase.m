% InputAmp=Make2DPolPhase(ImageParam,PSFParam) : makes a 2D complex valued amplitudes according the the PSFParam.polarization and PSFParam.phasemask settings
% ImageParam: Parameter set for the image
% PSFParam: Parameter set for the PSF. Choices are: 
%     PSFParam.polarization='circular'
%     PSFParam.polarization='linearX'
%     PSFParam.polarization='linearY'
%     PSFParam.polarization='radial'
%     PSFParam.polarization='azimuthal'
%     PSFParam.phasemask='flat'
%     PSFParam.phasemask='spiral'   : 1 turn of a phase spiral
%     PSFParam.phasemask='spiral2'  : 2 turns of a phase spiral
%
% default: PSFParam.polarization='circular', PSFParam.phasemask='flat'
%
% ImageParam=struct('Sampling',[20 20 50],'Size',[128 128 32]);
% PSFParam=struct('NA',1.2,'n',1.33,'lambdaEm',520,'Aplanatic',-1);  
% InputAmp=Make2DPolPhase(ImageParam,PSFParam)
% 
% Copyright (C) 2020, Rainer Heintzmann,  All rights reserved.
%_________________________________________________________________________
function InputAmp=Make2DPolPhase(ImageParam,PSFParam)
if ~isfield(PSFParam,'polarization')
    PSFParam.polarization = 'circular';
end
if ~isfield(PSFParam,'phasemask')
    PSFParam.phasemask = 'flat';
end

switch PSFParam.polarization
    case 'circular'
        InputAmp=newim([ImageParam.Size(1:2),1,2],'scomplex');
        InputAmp(:,:,:,0)=sqrt(0.5);
        InputAmp(:,:,:,1)=i*sqrt(0.5);   % Makes it circular polarisation (or also == random)
    case 'linearX'
        InputAmp=newim([ImageParam.Size(1:2),1,2],'scomplex');
        InputAmp(:,:,:,0)=1;
    case 'linearY'
        InputAmp=newim([ImageParam.Size(1:2),1,2],'scomplex');
        InputAmp(:,:,:,1)=1;
    case 'radial'
        InputAmp=newim([ImageParam.Size(1:2),1,2],'scomplex');
        InputAmp(:,:,:,0)=cos(phiphi(ImageParam.Size(1:2)));
        InputAmp(:,:,:,1)=sin(phiphi(ImageParam.Size(1:2)));
    case 'azimuthal'
        InputAmp=newim([ImageParam.Size(1:2),1,2],'scomplex');
        InputAmp(:,:,:,0)=sin(phiphi(ImageParam.Size(1:2)));
        InputAmp(:,:,:,1)=cos(phiphi(ImageParam.Size(1:2)));
    otherwise
        error('unknown polarization mode. Available are circular, linearX, linearY, radial, azimuthal, spiral')
end

switch PSFParam.phasemask
    case 'flat'
        % do nothing
    case 'spiral'
        InputAmp = InputAmp * exp(1i*phiphi(phiphi(ImageParam.Size(1:2))));
    case 'spiral2'
        InputAmp = InputAmp * exp(1i*phiphi(2*phiphi(ImageParam.Size(1:2))));
     otherwise
        error('unknown polarization mode. Available are flat, spiral, spiral2.')
end
