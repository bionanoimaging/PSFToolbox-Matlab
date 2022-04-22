% res=MidPosY(anImg)  : determines the mid Y position of an image in the dipimage sense = floor(size/2)

function res=MidPosY(anImg)
res=floor(size(anImg,2)/2);
