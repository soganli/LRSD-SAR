function [cost, X,Sk, Bk] = lowRankSparseSar_Decomposition(R,Rt,X,Sk,Bk,y,f,Ph,H,HH,f_conv,parameters)

lambda = parameters.lambda;
beta = parameters.beta;

c_sparse = parameters.c_sparse;
c_rank = parameters.c_rank;

maxiter = parameters.maxiter;
N = parameters.N;
truth = parameters.x;
cost = zeros(maxiter,1);
p = parameters.p;
error = 1;
iter = 1;
while error > parameters.cgtol  || iter < 20
    iter = iter + 1;
    MSE2_Learned = (1/(N*N)) *abs(sum((abs((X))-abs(truth(:))).^2));
    SNR_learned = 10*log10(var(abs(truth(:))) / MSE2_Learned);

    
    
    %% X - subproblem
    XOLD = X;
    [X, ~, ~, ~] = cg_solver_X(Rt,H,HH,Ph,X,Sk,Bk,y, speye(N*N,N*N), parameters);
    X = abs(X);
    
    
    %% S - subproblem
    Y = X - Rt(Bk) + ( lambda / beta );
    Sk = sign(Y) .* max(abs(Y) - c_sparse/beta, 0);
    Sk(Sk < 0) = 0;
    
    
    %% B - subprolbme
    Y = R(X) - R(Sk) +  ( R(lambda) / beta );
    [U,D,VT] = svd(Y); VT = VT';
    D = diag(D);
    ind = find(D > (c_rank/beta).*D.^(p-1));
    Dold =D;
    D = diag(D(ind) - (c_rank/beta).*D(ind).^(p-1));
    Bk = U(:,ind) * D * VT(ind,:); 
    %% P - subproblem
    
    
    [Phk] = phase_reconstruction(H,HH,Ph,X,y, speye(N*N,N*N), Ph, parameters);
    Ph = exp(sqrt(-1)*(angle(Phk)));
    X = abs(Ph.*X);
    
    %% lambda - subproblem
    
    
    lambda = lambda + beta*(X - Rt(Bk) - Sk);
    parameters.lambda = lambda;
    
    parameters.beta = parameters.beta * parameters.betaIncrement;
    beta = parameters.beta;
    
    
    e = y - H(Ph.*X);
    
    
    indexes = find(abs(f) > 0.9);
    figure(1);
    
    if strcmp(parameters.type ,'Synthetic')
        figure(1);
        drawnow;
        colormap(gray);
        subplot(241)
        imagesc(abs(reshape(f,N,N)));  axis image; title('Truth');
        subplot(242)
        imagesc(reshape((abs(X)),N,N));  axis image;    title('S + L');
        subplot(243)
        imagesc(reshape(((abs(Sk))),N,N));  axis image; title('Sparse(S)');
        subplot(244)
        imagesc(reshape(((abs(Rt(Bk)))),N,N));  axis image; title('Low-rank (L)'); 
        subplot(245)
        imagesc(reshape(abs(f_conv),32,32));      axis image; title('Conventional');        
        subplot(246)
        plot(parameters.trueSingular); title('True Singular');
        subplot(247)
        plot(Dold); hold on;   plot((c_rank/beta)*ones(length(Dold)),'r'); hold off;    title('Sing Val (L). and Thresh.')
        subplot(248)        
        plot(angle(f(indexes))); hold on; title('CorrectPhase(blue) vs Estimated (red)')
        plot(angle(Phk(indexes)),'r');  hold off; axis([0 15 -4 4]);
    else
        figure(1);
        drawnow;
        colormap(sar_cmap);
        subplot(241)
        imagesc(20*log10(abs(reshape(f,N,N))));    axis image; title('Truth');
        subplot(242)
        imagesc(20*log10(reshape((abs(X)),N,N)));   axis image;    title('S + L');
        subplot(243)
        imagesc(20*log10(reshape(((abs(Sk))),N,N)));   axis image; title('Sparse(S)');
        subplot(244)
        imagesc(20*log10(reshape(((abs(Rt(Bk)))),N,N)));    axis image; title('Low-rank (L)');
        subplot(245)
        imagesc(reshape(abs(f_conv),32,32));    axis image; title('Conventional');
        subplot(246)
        plot(parameters.trueSingular); title('True Singular');
        subplot(247)
        plot(Dold); hold on;   plot((c_rank/beta)*ones(length(Dold)),'r'); hold off;    title('Sing Val (L). and Thresh.')
        subplot(248)        
        plot(angle(f(indexes))); hold on;  title('CorrectPhase(blue) vs Estimated (red)')
        plot(angle(Phk(indexes)),'r');  hold off; axis([0 15 -4 4]);
    end

     error = norm(X -XOLD ,2) / norm(XOLD,2);
     fprintf('Iteration:%d\t mse :%f\t snr:%f\t error:%f\t beta:%f\t \n' , iter, MSE2_Learned,SNR_learned,error,beta);   
     
    %waitbar(iter/maxiter);
    
end

end