% Plot deconvolution results

% Deconvolved maps:
figure
pcolor(X,Y,gpl_x)
figure
pcolor(X,Y,gpbb_x)
figure
pcolor(X,Y,fgp_x)

% Norm of projected gradient:
figure
loglog(gpl_info.grad_norm/gpl_info.grad_norm(1),'b'), hold on 
plot(gpbb_info.grad_norm/gpbb_info.grad_norm(1),'r')
plot(fgp_info.grad_norm/fgp_info.grad_norm(1),'k')
legend('GPL','GPBB','FGP')  

% Objective function value:
gplp = (gpl_info.obj-gpl_info.obj(end))./norm(zeros(N,N)-gpl_x,'fro')^2;
gpbbp = (gpbb_info.obj-gpbb_info.obj(end))./norm(zeros(N,N)-gpbb_x,'fro')^2;
fgpp = (fgp_info.obj-fgp_info.obj(end))./norm(zeros(N,N)-fgp_x,'fro')^2;

figure
loglog(abs(gplp),'b'); hold on,
plot(abs(gpbbp),'r');
plot(abs(fgpp),'k');
L =  fgp_info.L;
plot(2*L./((1:opt.maxit)+1).^2,'k--')
axis([0 opt.maxit 1e-4 1e6])
legend('GPL','GPBB','FGP')  