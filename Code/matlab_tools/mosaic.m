% img=mosaic(imgs,layout) : assembles a series of images along a mosaic, and array of tiles.
% imgs: input images to assemble into mosaic
% layout: How many tiles should be used along X and Y direction
%
% Example:
% a=readim,b=readim('orka');
% mosaic(cat(3,a,b,b,a,a,b),[3 2])
%
function img=mosaic(imgs,layout)
sz=size(imgs);
syx=sz(1:2);
if nargin < 2
    layout=[1 1];
    layout(1)=ceil(sqrt(sz(end)));
    layout(2)=ceil(sz(end)/layout(1));
end

img=newim(sz(1)*layout(1),sz(2)*layout(2));
pz=0;
for py=1:layout(2)
    for px=1:layout(1)
        pxs=(px-1)*sz(1);
        pys=(py-1)*sz(2);
        if pz < sz(3)
            img(pxs:pxs+sz(1)-1,pys:pys+sz(2)-1)=imgs(:,:,pz);
        end
        pz=pz+1;
    end
end
