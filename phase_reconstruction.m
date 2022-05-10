function [ Ph ] = phase_reconstruction(H,HH,Ph,X,y,M, betaks, parameters)
%PHASE_RECONSTRUCTÝON Summary of this function goes here
%   Detailed explanation goes here

%lambda = parameters.phaselambda;
gamma = parameters.cgtol;
err=1e5;
iter = 0;
%e = H(Rt(X).*Ph) - y;
%cost(iter) = sum(abs(e(:))).^2 + lambda*sum((abs(Ph)-1).^2)
%initPh = Ph;
while err>gamma,

    iter = iter + 1;
    
    fkp1 = cg_phase(H,HH,Ph,X,y,M, betaks, parameters);
    
    err=norm(fkp1-Ph,2)/norm(Ph,2)    ;
    Ph=fkp1;
    
    
    %e = H(Rt(X).*Ph) - y;
    %cost = sum(abs(e(:))).^2 + parameters.phaselambda*sum((abs(Ph)-1).^2)
    
    betaks =exp(sqrt(-1)*(angle(fkp1)));
    
    
end

