function [x,info] = gpbb(fun,PSF,b,x0,opt)
% GPBB GRADIENT PROJECTION ALGORITHM WITH BB-STEPS
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
%   Kim, Dongmin, Suvrit Sra, and Inderjit S Dhillon. 2013. 
%   'A Non-Monotonic Method for Large-Scale Non-Negative Least Squares.' 
%   Optimization Methods and Software 28 (5) (October): 1012-1039.
%

start_time = tic;

if opt.verbose
   fprintf('GPBB ')
end

x = x0;

% Precompute fft of PSF
Fps = fft2(PSF);
FpsT = fft2(rot90(PSF,2));

fgx = @(x) fun(PSF,b,x,Fps,FpsT);      
[f,grad] = fgx(x);

% Diminishing scalars
xref = x;
fref = f;
gradref = grad;
beta = 1;           % no scaling at the beginning,

% Projected gradient
grad(x == 0 & grad > 0) = 0;
gradold = grad;
ngrad = norm(grad,'fro');
ngrad1 = ngrad;
info.grad_norm(1) = ngrad1;
info.obj(1) = f;

n = 0;
while (n < opt.maxit) && (opt.tol < ngrad/ngrad1)

   if opt.verbose
       if (mod(n, 30) == 0)
           fprintf('.');
       end
   end
   
   n = n+1;
   
   if mod(n,opt.Mds) == 0
       [f,grad] = fgx(x);
       if (fref - f) < opt.sigma*(dot(gradref(:),xref(:)-x(:)))
           beta = beta*opt.eta;
           disp('Diminishing scalars...')
       else
           xref = x;
           fref = f;
           gradref = grad;
       end
   end
   
   x(x == 0 & grad > 0) = 0;
   grad(x == 0 & grad > 0) = 0;
   gradold(x == 0 & grad > 0) = 0;
   
   % Compute BB-steps
   g = fftshift(ifft2(fft2(gradold).*Fps));  % Auxillary vector
   
   if (mod(n, 2) == 0)
       alpha = dot(gradold(:),gradold(:))/dot(g(:),g(:));
   else
       numer = dot(g(:),g(:));
       g = fftshift(ifft2(fft2(g).*FpsT));
       g(x == 0 & grad > 0) = 0;
       alpha = numer/dot(g(:),g(:));
   end
   
    x = x - beta*alpha * grad;      % Step along search path
    x(x < 0) = 0;                  
    
    [f,grad] = fgx(x);         
    
    grad(x == 0 & grad > 0) = 0;
    
    gradold = grad;
    
    info.obj(n+1) = f;
    ngrad = norm(grad,'fro');;
    info.grad_norm(n+1) = ngrad;
    info.time(n) = toc(start_time);
end
if opt.verbose
    fprintf('\n\t'); disp([num2str(n) ' iterations in ' num2str(info.time(n)) ' seconds'])
end
info.iter = n;
end