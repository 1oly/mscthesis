function [x,info] = fgp(fun, PSF, b, x0, opt)
% FGP FAST GRADIENT PROJECTION ALGORITHM
% 
% Usage:  [x,info] = fgp(@fun, PSF, b, x0, opt)
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
% Beck, Amir, and Marc Teboulle. 2009. 
% 'A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems.' 
% SIAM Journal on Imaging Sciences 2 (1) (January): 183-202.
%
start_time = tic;

if opt.verbose
   fprintf('FGP ') 
end

% Initialize variables
x = x0;
xold = x;
y = x;
t = 1; 

% Precompute fft of PSF
Fps = fft2(PSF);
FpsT = fft2(rot90(PSF,2));

% Evaluate f(x0)
fgx = @(x) fun(PSF,b,x,Fps,FpsT);      

if opt.backtrack
    L = 1;
    beta = 1.5;
else
    L = lipschitz(PSF,Fps);
end

% For n = 1
[fx,gradx] = fgx(x);
grady = gradx;
fy = fx;
gradx(x == 0 & gradx > 0) = 0;
ngradx = norm(gradx,'fro');

% Get info
info.grad_norm(1) = ngradx;
ngrad1 = ngradx;
info.obj(1) = fx;
info.sumx(1) = sum(x(:));

if opt.stop_cond
    cond = ngradx/ngrad1;
else
    cond = inf;
end

% Start iteration
n = 0;
while (n < opt.maxit) && (opt.tol < cond)
    
    if opt.verbose
        if (mod(n, 30) == 0)
            fprintf('.');
        end
    end
    
    n = n+1;
    
    x = max(0,y - (1/L)*grady);
    
    % Backtracking
    if opt.backtrack
        fx = fgx(x);
        m = 1;
        while fx > fy + dot(x(:)-y(:),grady(:))+(L/2)*norm(x-y,'fro')^2;
            L = m^(beta)*L;
            x = max(0,y - (1/L)*grady);
            fx = fgx(x);
            m = m + 1;
        end    
    end
    
    if opt.descent
        if fgx(x)>fgx(xold)
            x = xold;
        end
    end
    
    tnew = (1+sqrt(1+4*t*t))/2;
    y = x + ((t-1)/tnew)*(x-xold);
    
    [fy,grady] = fgx(y);
    xold = x;
    t = tnew;
    
    if opt.stop_cond && (mod(n, opt.gradn) == 0)
        [fx,gradx] = fgx(x);
        gradx(x == 0 & gradx > 0) = 0;
        ngradx = norm(gradx,'fro');
        cond = ngradx/ngrad1;
        
        info.obj(n+1) = fx;     
        info.grad_norm(n+1) = ngradx;
    elseif opt.stop_cond==0
        cond = inf;
    else
        info.obj(n+1) = fx;     
        info.grad_norm(n+1) = ngradx;
    end
    
    if opt.backtrack
        info.L(n) = L;
    else
        info.L = L;
    end
    
    info.time(n) = toc(start_time);
end

if opt.verbose
    fprintf('\n\t'); disp([num2str(n) ' iterations in ' num2str(info.time(n)) ' seconds'])
end

info.iter = n;

end


function L = lipschitz(PSF,Fps)
% Estimate Lipschitz constant by power iteration
x = rand(size(PSF));
for k = 1:10
    x = fftshift(ifft2(fft2(x).*Fps))/norm(x,'fro');
end
    L = norm(x,'fro')^2;    % lambda(A'A) Assuming a symmetric matric A
end
