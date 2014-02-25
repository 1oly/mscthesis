% Get beamformer map
clear all
close all
load('../data/pizza_array4.mat');   % Load microphone positions
rn = pizza_array;
N = 70;
z0 = 1;
phi = 30;
f = 1500;

source = [27 35; 43 35];    % x,y position of sources

[b,PSF,X,Y,x,y] = psf(N,z0,f,phi,rn,source);

%figure
%pcolor(X,Y,real(b))