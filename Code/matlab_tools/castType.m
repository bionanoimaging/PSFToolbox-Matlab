% function res = castType(data,type)
% change the type or class of the input_data into output_type and results
% 'res' is the same data as the input_data but with the output_type
%  This function is needed because Matlab cannot converst directly double
%  to DipImage
% 21.04.21

function res = castType(input_data,output_type)
switch output_type
    case 'dip_image'
        res=dip_image(input_data);
    case 'double'
        res=double(input_data);
    otherwise
        res=cast(input_data,output_type);
end