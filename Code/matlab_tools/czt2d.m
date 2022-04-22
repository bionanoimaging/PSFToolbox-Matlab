% xout=czt2d(xin,scalx , scaly , idir ,method) : makes 2d chirped z transforms for 3d or 4d data
% xin: input image to make serial 2d transforms of
% scalx:  Zoom factor along X 
% scaly: Zoom factor along y
% idir: 0: x and y, 1: only x, 2: only y
% method: 'forward', 'inverse'
%
function xout=czt2d(xin,scalx , scaly , idir ,method)
if nargin < 4
    idir=0;
end
if nargin < 5
    method='forward';
end
sz=size(xin);
xout=[];

if ndims(xin)==4
 for jj=1:sz(4)
    for ii=1:sz(3)
        xout(:,:,ii,jj)=czt_dip(squeeze(xin(:,:,ii-1,jj-1)),scalx , scaly , idir, method);
    end
 end
elseif ndims(xin)==3
    for ii=1:sz(3)
        xout(:,:,ii)=czt_dip(squeeze(xin(:,:,ii-1)),scalx , scaly , idir, method);
    end
elseif ndims(xin)==2
        xout=czt_dip(xin,scalx , scaly , idir, method);
else
    error('czt2d can only be applied to 3d and 4d data')
end
xout=dip_image(xout);