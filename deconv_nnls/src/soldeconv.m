function [x,info] = soldeconv(fun,PSF,b,x0,opt)
%SOLDECONV Solve deconvolution problem b = q * psf
%
% Usage:  [x,info] = soldeconv(@fun,PSF,b,x0,opt)
%
% Solve deconvolution problem by following methods:
%  GPL: Gradient projection with exact line search
%  GPBB: Gradient projection with bb-steps
%  FGP: Fast gradient projections
%
% Input:
%   @fun: Function handle with objective function and gradient
%   PSF: Point-spread function
%   b : Beamformer map
%   x0: Starting vector
%   opt: Struct with options - see deconv_opt.m for default settings
%
%
% Output: 
%   x: Source distribution for all iterations
%   info: Struct with info about
%       - info.obj
%       - info.grad_norm
%       - info.time
%       - info.iter
%       - info.sumx
%       - info.L
%       - info.x_warm
%
% Author: Oliver Lylloff
% Date: 5/2/14
% Latest revision: 5/2/14
%
% http://github.com/1oly/mscthesis

if opt.zeropad
    N = size(b,1);
    npad = 2*N;
    b = zeropad(b,npad,npad,'pad');
    PSF = zeropad(PSF,npad,npad,'pad');
    x0 = zeropad(x0,npad,npad,'pad');
end

switch opt.algo
    
    case 'GPL'
        [x,info] = gpl(fun, PSF, b, x0, opt);
        if opt.zeropad
            x = zeropad(x,N,N,'unpad');
        end
        
    case 'FGP'
        [x,info] = fgp(fun, PSF, b, x0, opt);
        if opt.zeropad
            x = zeropad(x,N,N,'unpad');
        end
        
    case 'GPBB'
        [x,info] = gpbb(fun, PSF, b, x0, opt);
        if opt.zeropad
            x = zeropad(x,N,N,'unpad');
        end
        
    otherwise
        error('Unknown method')
end

end