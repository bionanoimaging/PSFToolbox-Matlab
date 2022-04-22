% res=MidPosX(anImg)  : determines the mid X position of an image in the dipimage sense = floor(size/2)

function res=MidPosX(anImg)
res=floor(size(anImg,1)/2);
