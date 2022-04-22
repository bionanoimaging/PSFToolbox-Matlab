% FillProjSphere(img,aproj,indexList2D,fullIndex3D,factorList): Fills a 2D amplitude array into 3D Fourier space using an interpolation matrix

function img=FillProjSphere(img,aproj,indexList2D,fullIndex3D,factorList)

vallist=aproj(indexList2D);
vallist=reshape(vallist,[size(vallist,1) 1]);
if isa(vallist,'cuda')
    vallist=repmat(vallist,[1 size(factorList,2)]);
end
vallist=factorList .* vallist; % outer product

img(fullIndex3D)=img(fullIndex3D)+vallist(:);

