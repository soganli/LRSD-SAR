function x=HH(y,parameters)

%no normalization i.e. produces x=F'*y, where F is DFT.
N1 = parameters.N;
N2 = parameters.N;
N3 = parameters.Nlim;
N4 = parameters.Nlim;

y=reshape(y,N3,N4);

ytemp=zeros(N1,N2);

ytemp(1:N3/2,1:N4/2)=y(1:N3/2,1:N4/2);
ytemp(1:N3/2,N2-N4/2+1:N2)=y(1:N3/2,N4/2+1:N4);
ytemp(N1-N3/2+1:N1,1:N4/2)=y(N3/2+1:N3,1:N4/2); 
ytemp(N1-N3/2+1:N1,N2-N4/2+1:N2)=y(N3/2+1:N3,N4/2+1:N4);

x=ifft2(ytemp)*N1*N2;
x=x(:) / N1;
