% function xout = iczt(xin, scale, dims) compute the inverse czt of input image xin
%                          along the dimension stored in vector dims and with the scale given in the
%                          vetor scale
% see also: czt_1d, czt
% 
% Chirp z transform along a single direction d of an ND array `xin` into the ND array 'xout'.
% Note that xin and xout can be the same array for inplace operations.
% This code is based on a 2D Matlab version of the CZT, written by H. Gross et al.
%     
% References: Rabiner, Schafer, Rader, The Cirp z-Transform Algorithm, IEEE Trans AU 17(1969) p. 86
% Copied and translated by Dina R from Julia code by Rainer Heintzmann - Aug 2021

function xout = iczt( xin , scale, dims)
    sz = size(xin);
    if nargin<3
        dims = 1:length(sz);
    end
    xout = conj(czt(conj(xin), scale, dims)) / prod(size(xin));
end
