% res=rrscale(mysize,myscales): a scaled version of rr with sampling being different along the coordinates
% mysize: size of image to generate
% myscales: pixel size vector
%
% Example:
% rrscale([20 20],[1 2])

function res=rrscale(mysize,myscales)
if nargin < 1
    mysize=[256 256];
end
if nargin < 2
    myscales=mysize*0+1;
end
if numel(myscales) < numel(mysize)
    myscales(numel(myscales)+1:numel(mysize))=0;
end

res = (ramp(mysize,1).*myscales(1))^2;
for d=2:numel(mysize)
    res = res + (ramp(mysize,d).*myscales(d))^2;
end
res = sqrt(res);
