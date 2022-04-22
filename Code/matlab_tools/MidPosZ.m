% res=MidPosZ(anImg)  : determines the mid Z position of an image in the dipimage sense = floor(size/2)

function res=MidPosZ(anImg)
res=floor(size(anImg,3)/2);
