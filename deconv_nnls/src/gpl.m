function [x,info] = gpl(fun, PSF, b, x0, opt)
% GPL GRADIENT PROJECTION ALGORITHM WITH EXACT LINE SEARCH
%
% Usage:  [x,info] = gpbb(@fun, PSF, b, x0, opt)
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
% 
% Reference: 
% Ehrenfried, Klaus, and Lars Koop. 2008. 
% 'A Comparison of Iterative Deconvolution Algorithms for the Mapping of Acoustic Sources.' 
% American Institute of Aeronautics and Astronautics.

start_time = tic;

if opt.verbose
   fprintf('GPL ') 
end

x = x0;

% Precompute fft of PSF
Fps = fft2(PSF);
FpsT = fft2(rot90(PSF,2));

fgx = @(x) fun(PSF,b,x,Fps,FpsT);      
[~,grad] = fgx(x);
w = grad;                 
w(x == 0 & w > 0) = 0;  
ngrad = norm(w,'fro');
ngrad1 = ngrad;
info.grad_norm(1) = ngrad1;
info.sumx(1) = sum(x(:));
n = 0;
while (n < opt.maxit) && (opt.tol < ngrad/ngrad1)
    
    if opt.verbose
        if (mod(n, 30) == 0)
            fprintf('.');
        end
    end
    
    n = n+1;
    [f,grad,r] = fgx(x);        
    w = grad;                   
    w(x == 0 & w > 0) = 0;      
    
    g = fftshift(ifft2(fft2(w).*Fps));         
    alpha = dot(g(:),r(:))/dot(g(:),g(:));    
    
    x = x - alpha*w;           
    x(x < 0) = 0;              
    
    info.obj(n) = f;
    ngrad = norm(w,'fro');
    info.grad_norm(n+1) = ngrad;
    info.time(n) = toc(start_time);
end
if opt.verbose
    fprintf('\n\t'); disp([num2str(n) ' iterations in ' num2str(info.time(n)) ' seconds'])
end
info.iter = n;
end