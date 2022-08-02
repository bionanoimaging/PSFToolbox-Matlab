% [amp3d,h]=Proj3DSincR(ImageParam,PSFParam,FPlane,BorderRegion) : Uses the sinc(abs(rr)) trick and maps the BFP onto it
% ImageParam : structure of image parameters (see GenericPSFSim)
% PSFParam: structure of PSF parameters (see GenericPSFSim)
% FPlane: three electric field pupils to be projected
% BorderRegion: If active (i.e. a value is attributed to it), the region is used for DampEdge (default: 0.1)
% 
% Copyright (C) 2020, Rainer Heintzmann,  All rights reserved.
%_________________________________________________________________________
function [amp3d,h]=Proj3DSincR(ImageParam,PSFParam,FPlane,BorderRegion)
    if nargin < 4
        BorderRegion=0.2;
    end
% detect whether the k-sphere goes out of range
    KZExtent = ImageParam.Size(3) .* (2*PSFParam.n/PSFParam.lambdaEm).*ImageParam.Sampling(3)/2.0;
    if KZExtent/2 > ImageParam.Size(3)/2-1
         fprintf('k-Sphere kZ extent: %g, Sampling Limit Z: %g\n',KZExtent,floor(ImageParam.Size(3)/2-1));
         error('K-sphere undersampled')
    end
    
    myzz=zz(ImageParam.Size);
    myrr=rrscale(ImageParam.Size,ImageParam.Sampling *2*pi/(PSFParam.lambdaEm/PSFParam.n)); % This is in real space
    amp3d=sinc(myrr); % sinc is sin(pi x) /(pi x)
    % Edit Dina 17.06.22 -----------------------------
%     Radius=(ImageParam.Size(1)./(1+BorderRegion))./2;
%     Disc = disk(ImageParam.Size(1:2),[0 0],Radius);
%     amp3d=amp3d.*Disc; % limit the sinc by a disc of radius equal to (2/3)*bigsize (real space)
    % end edit ------------------------------------------

    % Edit Dina 26.05.22 -----------------------------
%     Radius=1/2*ImageParam.Size(1);
%     Radius=ImageParam.Size(1:2)./(1+BorderRegion);%BorderRegion*ImageParam.Size(1:2);
%     Rec=extract(newim(round(Radius))+1,ImageParam.Size(1:2));
%     amp3d=amp3d.*Rec; % limit the sinc by a disc of radius equal to (2/3)*bigsize (real space)
    % end edit ------------------------------------------
    
    % The below approaches destoy the uniform brightness property and should NOT be used!
    % amp3d=DampEdge(amp3d,0.1,2,0,1);
    % amp3d=DampEdge(amp3d,BorderRegion,3,1,3);
    % amp3d=DampEdge(amp3d,BorderRegion,3,1,1);   % dimming down to zero seems a good choice here and creates less artefacts.
    % START
    ftAmp=ft3d(amp3d);  % 3D shell
%     ftAmp=
    % END
    if KZExtent > ImageParam.Size(3)/2-1
        fprintf('k-Sphere kZ extent: %g, Sampling Limit Z: %g\n',KZExtent,floor(ImageParam.Size(3)/2-1));
        fprintf('Warning!! kZ sphere is undersampled, but the final intensity OTF will still be OK\n');
        LocalKXY2 = (abssqr(xx(ImageParam.Size(1),'freq')./ImageParam.Sampling(1))+abssqr(yy([1 ImageParam.Size(2)],'freq')./ImageParam.Sampling(2))) ./ (PSFParam.n/PSFParam.lambdaEm).^2;
        LocalKZ = ImageParam.Size(3) * sqrt(1-LocalKXY2) *ImageParam.Sampling(3)/(PSFParam.lambdaEm/PSFParam.n);
        MZ=MidPosZ(ftAmp);
        res=newim(ftAmp,'complex');
        if MZ*2==size(ftAmp,3)
            AliasMask = LocalKZ >  MZ-1;
            res(:,:,0:MZ-1)=ftAmp(:,:,MZ:end)*~AliasMask;
            res(:,:,MZ:end)=ftAmp(:,:,0:MZ-1)*AliasMask;
            ftAmp=res;
            ftAmp(:,:,MZ)=res(:,:,MZ)/2;
            ftAmp(:,:,0)=res(:,:,0)*2;
        else
            AliasMask = LocalKZ >  MZ;
            res(:,:,0:MZ)=ftAmp(:,:,MZ:end)*~AliasMask;
            res(:,:,MZ+1:end)=ftAmp(:,:,0:MZ-1)*AliasMask;
            ftAmp=res;
        end
        clear res;
    else
        ftAmp(myzz==0)=ftAmp(myzz==0)/2;  % keep the middle plane but reduce by a factor of two
        ftAmp(myzz>0)=0;  % keep the middle plane 
        fPlaneSum=sum(ftAmp,[],3);
    end
    ftAmp=ftAmp .* FPlane;  % this used to be with the aplanatic factor +2 from the given one
    %ftAmp=ftAmp .* (FPlane.*conj(fPlaneSum) ./ (abssqr(fPlaneSum)+1e-10));
    amp3d=ift3d(ftAmp);
    if KZExtent > ImageParam.Size(3)/2-1
        amp3d = amp3d * exp(pi*1i*zz([1 1 ImageParam.Size(3)])); % shift it back (with aliasing)
    end
    h=squeeze(sum(abssqr(amp3d),[],4));
end
