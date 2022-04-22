% [indexList2D,fullIndex3D,factorList,amask]=FillProjSpherePrepare(size3d,lambda,pixelsize,NA,imatrix,n,allowwraparound): Calculates preperatory data for the FillPrjSphere function
% size3d: size of the 3d array to later fill in the data
% lambda: vakuum wavelength
% pixelsize : size of a pixel in real space (same units as lambda)
% NA=sin(alpha) : numerical aperture for refractive index 1.0. Else use NA/n
% imatrix : the interpolation matrix as a dipimage ranging from 0 to 1 subpixel difference and having the interpolation kernel in each row.
% n: refractive index of the medium
% allowwraparound: If true, interpolation can wrap around along Z
%
% Use "IterateCoefficients" to generate this matrix.
% The output arguments (indexList2D,fullIndex3D,factorList) are all one-dimensional index lists into the pixels to access
% indexList2D : indexes the pixels in the 2d projection (which is supposed to have the same XY size as the 3D array)
% indexList 3D : indexes the position where the write the interpolation results
% factorList : contains the interpolation weights for the appropriate pixel
% amask : the 2D pupil mask used
%
% Copyright (C) 2020, Rainer Heintzmann,  All rights reserved.
%_________________________________________________________________________
%
function [indexList2D,fullIndex3D,factorList,amask]=FillProjSpherePrepare(size3d,lambda,pixelsizes,NA,imatrix,n,allowwraparound)
if nargin < 6
    n=1.0;
end
if nargin < 7
    allowwraparound=1.0;
end
mid3d=floor(size3d/2);
size2d=size3d(1:2);
k0=1/(lambda/n);
kxymax=k0*(NA/n); % Note [06.07.21]: The n here correspond to the ri in immersion medium region. kxymax here is emitter dependent i.e. this is what the emitter but what is collected by the objective lens is soomething else. The maximum angular aperture is defined by the ri in design condition. The jinc aperture will take care of the maximum allowed light into the system.
myr2=(xx(size2d,'freq')/pixelsizes(1)).^2+(yy(size2d,'freq')/pixelsizes(2)).^2;
kz=newim(myr2);
tmask=k0*k0 > myr2;  % generate the fukll k-circle in 2D out of the existing k-space sphere
if length(size3d) < 3
    size3d(3)=1;
end
kz(tmask)=-sqrt(k0*k0-myr2(tmask))*pixelsizes(3)*size3d(3);  % calculate the hight map for this cap of the sphere
% changed kz to be negative to correspond to an upward propagating wave.

% amask=myr2 < kxymax*kxymax;  % limit to the aperture chosen by the user
% amask=myr2 < (kxymax*kxymax + k0*k0)/2;  % limit to the aperture chosen by the user, but made intentionally larger to account for the edges of the pupil
amask=sqrt(myr2) < (1*kxymax + 3*k0)/4;  % limit to the aperture chosen by the user, but made intentionally larger to account for the edges of the pupil

if (prod(size(imatrix))==1) && (imatrix==1) && (length(size3d) < 3 || size3d(3)<=1)
    indexarray=newim(size2d);
    indexarray(:)=[0:prod(size(indexarray))-1];
    indexList2D=double(indexarray(amask));  % To be used for DipImage indexing
    fullIndex3D=indexList2D;
    factorList=1;
    return;
end
if (mid3d+k0+(size(imatrix,1)+1)/2) > size3d(3)
    error('Datavolume is not large enough along Z to accomodate all interpolation values');
end

% Write the start z-index into the 2D mask

kzarray=floor(kz);  % get the reference pixel position in Z at each XY position
matrixindex = round((kz - kzarray)*(size(imatrix,2)-1));  % this indicates the subpixel indexing. It denotes the row in the imatrix to use for interpolation

indexarray=newim(size2d);
indexarray(:)=[0:prod(size(indexarray))-1];
%indexarray=reshape([1:prod(size2d)],size2d);  % just to extract the 2d mask indices
indexList2D=double(indexarray(amask));  % To be used for DipImage indexing
%% The code below creates an index list into the 3d array
kzoffset=mid3d(3)-(size(imatrix,1)-1)/2;
if max(kzarray+kzoffset) >= size3d(3) || min(kzarray+kzoffset)<0
    if (0)
        error('Pupil kz larger than 3D data volume allows. PSF is undersampled')
    else
        if allowwraparound
            fprintf('WARNING: Pupil kz larger than 3D data volume allows. Allowing for interpolation to wrap around and cause potential aliasing by %d Z pixels.\n',min(kzarray+kzoffset));
            kzarray=mod(kzarray+kzoffset,size3d(3));
            kzoffset=0; % as it is already accounted for
        else
            kzoffset=-min(kzarray(amask));
            fprintf('WARNING: Pupil kz larger than 3D data volume allows. PSF is undersampled. Shifting by applying an offset of %d.\n',kzoffset)
            if max(kzarray(amask))+kzoffset + size(imatrix,1) >= size3d(3)
                error('Pupil kz larger than 3D data volume allows. The Problem can be undersampling of the PSF or a too small Z data volume, which also needs to have space for the interpolators to fit completely.')
                % fprintf('Warning: Pupil kz larger than 3D data volume allows when considering the interpolation range. This will generate wrap-around, but have hopefully no effect in the intensity.\n')
            end
        end
    end
end
indexlist3d=double((kzarray(amask)+kzoffset) * (size3d(1)*size3d(2)) + indexList2D);   % This is the start index of each line to write in 3D
fullIndex3D=transpose(repmat(indexlist3d,[size(imatrix,1) 1])) + (transpose(indexlist3d)*0+1) * [0:size(imatrix,1)-1]*(size3d(1)*size3d(2));
fullIndex3D=transpose(fullIndex3D);
fullIndex3D=mod(fullIndex3D(:),prod(size3d)); % To be used as DipImage indexing
%mask3d=newim(img,'bin');
%for z=1:size(imatrix,1)
%    mask3d(floor(indexlist3d))=1;
%end
%indexarray3d=reshape([1:prod(size(mask3d))],size(mask3d));  % just to extract the 2d mask indices
%fullindex3d=indexarray3d(mask3d);
%%

%See: 
%q=[1 2 3;4 5 6;7 8 9]
%qq=q([1,3,2,3,1,1],:)
% transpose([1 2 3 4]) * [10 1 10 10]  % Outer product
matrixindexlist=double(matrixindex(amask));
factorList=permute(imatrix(:,matrixindexlist),[2 1]);  % This returns the factors to be multiplied by the amplitude values in aproj

% The normalizations below make the projections of the 3D OTF more smooth, but also slow down convergence of the deconvolution
% factorList=factorList ./ sum(factorList,[],2);  % This normalization seems necessary to avoid a frequency-dependend non-uniformity of about 4.5 %
% factorList=factorList ./ sqrt(sum(factorList,[],2));  % This normalization seems necessary to avoid a frequency-dependend non-uniformity of about 4.5 %
% Now call

% FillProjSpherePrepare(img,indexList2D,fullIndex3D,factorList)

