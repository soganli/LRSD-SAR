function y=H(x,parameters)

N1 = parameters.N;
N2 = parameters.N;
N3 = parameters.Nlim;
N4 = parameters.Nlim;


x=reshape(x,N1,N2);

ytemp=fft2(x,N1,N2);
y=[ytemp(1:N3/2,1:N4/2) ytemp(1:N3/2,N2-N4/2+1:N2);ytemp(N1-N3/2+1:N1,1:N4/2) ytemp(N1-N3/2+1:N1,N2-N4/2+1:N2)];  

y=y(:) / N1;
