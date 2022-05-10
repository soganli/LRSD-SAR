function xc=addphase(x);

% adds random uncorrelated uniform [-pi,pi] phase to matrix x

[M,N]=size(x);
phi=rand(M,N);
phi=pi*(2*phi-1);
phase=exp(sqrt(-1)*phi);

xc=x.*phase;


%diag(phase) would be a matrix such that xc=diag(phase)*x;