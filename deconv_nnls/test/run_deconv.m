% Run deconvolution algorithms

plot_beamformer;    % Get beamformer map

opt = deconv_opt;   % Get options

% Gradient projection with exact line search:
opt.algo = 'GPL';
[gpl_x,gpl_info] = soldeconv(@nnlsqfun,PSF,real(b),zeros(size(b)),opt);

% Gradient projection with BB-steps:
opt.algo = 'GPBB';
[gpbb_x,gpbb_info] = soldeconv(@nnlsqfun,PSF,real(b),zeros(size(b)),opt);

% Fast gradient projection:
opt.algo = 'FGP';
[fgp_x,fgp_info] = soldeconv(@nnlsqfun,PSF,real(b),zeros(size(b)),opt);

