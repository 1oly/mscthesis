function [f,g,r] = nnlsqfun(PSF,b,x,varargin)
% NONNEGATIVE LEAST SQUARES OBJECTIVE
% 
% Compute objective, gradient and residual of 
% Nonnegative least squares problem (NNLS): 
%   minimize_x      0.5||Ax - b||^2 
%   subject to         x >= 0
%
% Usage:    
%       f = nnlsqfun(PSF,b,x)
%       [f,g] = nnlsqfun(PSF,b,x)
%       [f,g,r] = nnlsqfun(PSF,b,x)
%
% Input:
%   PSF: Point-spread function (matrix)
%   b : Beamformer map (matrix)
%   x: Source distribution (matrix)
%   varargin: Can take precomputed PSF
%
% Output: 
%   f: Objective function value (scalar)
%   g: Gradient (matrix)
%   r: Residual (matrix)
%
%
% Author: Oliver Lylloff
% Date: 7/2/14
% Latest revision: 5/2/14
%
% http://github.com/1oly/mscthesis


if nargin==5
    Fps = varargin{1};
    FpsT = varargin{2};
else
    Fps = fft2(PSF);
    FpsT = fft2(rot90(PSF,2));
end

r = fftshift(ifft2(fft2(x).*Fps)) - b;
f = 0.5*norm(r,'fro')^2;
if (nargout > 1)
    g = fftshift(ifft2(fft2(r).*FpsT));
end

end