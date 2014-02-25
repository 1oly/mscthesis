function options = deconv_opt
% DEFAULT OPTIONS STRUCT FOR DECONVOLUTION ALGORITHMS

% General settings:
options.maxit = 300;        % Max number of iterations
options.tol = 1e-4;         % Projected gradient norm
options.verbose = 1;        % Prints out progress to Command Window
options.zeropad = 1;        % Zero-padding

% FGP:
options.stop_cond = 1;      % 1: Use tol as stopping criterion else 0: just maxit
options.gradn = 1;         % Compute gradient of x every k iterations default = 20
options.backtrack = 0;      % Toggle backtracking on/off
options.descent = 0;        % Descent condition

% GPBB
options.sigma = 0.298;      % Diminishing scalars
options.eta = 0.99;
options.Mds = options.maxit/10;

end