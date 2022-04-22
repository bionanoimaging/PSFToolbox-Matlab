% res=MidPos(img,myDim) : returns the integer center position
% img : dataset to subslice
% myDim : dimension along which to subslice
%
function res=MidPos(img,myDim)
if nargin < 2
    myDim=[];
end
res=floor(size(img)/2);
if ~isempty(myDim)
    res=res(myDim);
end
