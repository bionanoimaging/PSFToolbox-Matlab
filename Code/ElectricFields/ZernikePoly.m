% function res = ZernikePoly(order,ImageParam,PSFParam) generates the
% Zernike polynomial of radial order n = order(1) and azimuthal order m =
% order(2). Zernike polynomial is in range of [-1,1].
% n is a nonnegative integer; m can be positive and negative but satisfies n-m=2k where k is a natural number i.e.  n-m is even and n>=m
% [n m]=[0 0] piston or bias 
%       [1 -1] and [1 1] tilt along x and y respectively
%       [2 -2]; [2 0]; [2 2] astigmatism, defocus
%       [3 -3]; [3 -1]; [3 1]; [3 3] coma trefoil
%       [4 -4]; [4 -2]; [4 0]; [4 2]; [4 4] tetrafoil, 2nd astigmatism, spherical, higher aberrations start from here
% order={[4 0 1 90];[1 -1 1 60]}; % a cell array where each input
%       corresponds to a given order of aberration [n m aberration_coeff azimuth_angle]; 
%       aberration_coeff (power of the aberration) is 1 by default and azimuth_angle in degrees is 0 (rotate the phase around z-axis from the vertical axis in anti-clockwis by this angle)
% 
% Source: Lakshminarayanan V, Fleck A. Zernike polynomials: a guide. Journal of Modern Optics. 2011 Apr 10;58(7):545-61.
% 
% ToDO: These Zernikes will not be very orthogonal for smaller array sizes. To fix this, one should supersample them first and then downsample by binning.
%
% example 1:
% ImageParam=struct('Sampling',[40 40 100],'Size',[256 256 16]);
% PSFParam=struct('NA',1.2,'n',1.33,'lambdaEm',520);
% order={[3,1];[1 -1 1 60]};
% res = ZernikePoly(order,ImageParam,PSFParam);
% 
% example 2:
% order={[3,1];[1 -1 1 60]};
% rho=rr(100,100); 
% mask=rho<max(rho)/4;
% res=ZernikePoly(order,rho,mask);


function [phse,rho,phi] = ZernikePoly(orderlist,ImageParam,PSFParam,Option)
if nargin<4 
    Option=0; % this option mean we want to have the resulted summation of Zernike. If its value is different than 0 then have each polynomial computed but do not summed them up. This last option is needed when decomposing a phase into Zernike polynomials
end
    
if ~iscell(orderlist)
    orderlist={[orderlist]};
end
if ~isstruct(ImageParam)
    mask=dip_image(PSFParam);
    rho=dip_image(ImageParam).*mask; % if ImageParam is not a structure then it should be the 2D rho 
    rho=rho/max(rho);
    origSZ=size(rho);
    midPos=(origSZ-mod(origSZ,2))/2;
    PupilRadius(1)=find(mask(midPos(1),midPos(2):end)~=1, 1, 'first');
    PupilRadius(2)=find(mask(midPos(1):end,midPos(2))~=1, 1, 'first');
    rmsk=mask;
else
    origSZ=ImageParam.Size(1:2);
    [IncoherentRadius,PupilRadius,EwaltRadii]=AbbeMaxRadiusFromPSF(PSFParam,ImageParam);

    rho=rrscale(ImageParam.Size(1:2), 1.0 ./ PupilRadius);
    rmsk=rho>1.0;
    rho(rmsk)=0.0;  % limit the radius to o ne. THIS IS PHASE! Not amplitude.
end
myxx = xx(origSZ)./PupilRadius(1);  % to get a valid coordinate system in Fourier space
myyy = yy(origSZ)./PupilRadius(2);
phi = atan2(myyy,myxx);   % only then are the pupil angles correct   

if Option==0
    phse = newim(origSZ); 
else
    phse=[];
end
for abnb=1:length(orderlist)%size(orderlist,1)  % iterate through the list of aberrations
    order=orderlist{abnb};
 
    if length(order)<4
        rotateAngle=0; 
    else 
        rotateAngle=order(4);
    end
    if rotateAngle ~= 0
        phi = phi + rotateAngle*pi/180;  % This will automatically rotate all the anglular orders
    end
    if length(order)<3
        a=1; % aberration coefficient
    else
        a=order(3);  % was a specific orientation provided?
    end
    n=order(1);   % get the radial order for this aberration
    m=abs(order(2)); % get the angular order
    p=n-m;
    if p<0
        error('The azumuthal order m should be less than the radial order n');
    end
    if mod(p,2)==1 % i.e. n-m = 2k+1 odd, k in N (natural integer)
        res=0;
    else
        l=0:1:p/2;   % to be summed over
        l=expanddim(dip_image(l),3); l=double(permute(l,[3 2 1])); 
        numerator=(-1).^l.*factorial(n-l).*dip_image(double(rho).^(n-2.*l));   % The cast is necessary since DipImge 2.9 has trouble with broadcasting in the .^ function
        denominator=factorial(l).*factorial((1/2).*(n+m)-l).*factorial((1/2).*(n-m)-l);
        R=numerator./denominator;
        if length(size(R))==3
            R=sum(R,[],3);
        end
        R=squeeze(R);
%         res=R*cos(m*phi);  % The angles is taken care of via the offset angle in phi. This means the odd j are not needed! For arbitrary angles you always have to add multiple zernikes.
        % add today +uncomment the above line
        if order(2)<0
            res=R*sin(m*phi);
        else
            res=R*cos(m*phi);
        end
        % end add
        res=a*res;%*rmsk;
    end
    if Option==0
        phse = phse + res; 
    else
        phse{abnb}=res;%cat(3,phse,res);
    end
end

% % phase outside the mask should be zero not non-zero
if Option==0
    phse=phse.*(~rmsk); 
else
    for k=1:length(phse)
        phse{k}=phse{k}.*(~rmsk);
    end
end
