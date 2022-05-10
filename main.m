clear all; close all;  tic

load testImage23;   f = b + s; 
f = f / max(f(:)); 
parameters.N = size(f,1);  %% 32
parameters.Nlim = 28;    %% Band-limitide Fourier geniþliði   (28/32)^2 = %75 veri var.
N = parameters.N;
Nlim = parameters.Nlim;
a = 1;
parameters.beta = 1;          parameters.betaIncrement = 1.02;   %% Beta parametresinin deðeri. Her iterasyonda yavaþ yavaþ arttýrýlýyor.
parameters.patch_sizex = 8;   parameters.patch_sizey = 8;       parameters.sliding_distance = 1;  % Patch deðerleri
parameters.cgiteration = 300; parameters.cgtol = 1e-5;    
parameteres.mu = 1;
parameters.c_sparse = 0.08;   parameters.c_rank = 1;   %% Sparse ve low-rank terimlerinin katsayýsý
parameters.maxiter = 300;
parameters.phaselambda = 0.05;   %% Faz kestiriminin katsayýsý (Sýfýr olunca daha iyi çalýþýyor nedense)
parameters.p = 1;  
parameters.trueSingular = svd(b/2);
f = addphase(f);

parameters.type = 'Synthetic';
parameters.x = f(:);



[X,indexes] = my_im2col(ones(N,N),parameters);  % We obtain the patches for dictionary process.  bb is 8: For 32x32 image there will be 625 patches.
[rows,cols] = ind2sub([N,N]-[parameters.patch_sizex, parameters.patch_sizey]+1,indexes);
count = 1;
Weight = zeros(N,N);
for j  = 1:length(cols)
    col = cols(j); row = rows(j);        
    Weight(row:row+parameters.patch_sizex-1,col:col+parameters.patch_sizey-1)=Weight(row:row+parameters.patch_sizex-1,col:col+parameters.patch_sizey-1)+ones(parameters.patch_sizex , parameters.patch_sizey);
    count = count + 1;
end;

parameters.patchWeight = Weight;
parameters.patchIndexes = indexes;
parameters.sizeofPatchMatrix = size(X);

parameters.lambda_rank = zeros(parameters.sizeofPatchMatrix);
parameters.lambda_sparse = zeros(N*N,1);
parameters.lambda = zeros(N*N,1);

R = @(z)R(z,parameters); % The forward Patch Extraction Operator
Rt=@(z)Rt(z,parameters); % The backward Image reconstruction (Default is taking mean of overlapping pixels)

H = @(z)H(z,parameters); % The forward Fourier sampling operator
HH=@(z)HH(z,parameters); % The backward Fourier sampling operator

sigma_noise = 0.01;
y = H(f) + sigma_noise*(randn(Nlim*Nlim,1)+sqrt(-1)*randn(Nlim*Nlim,1));
f_conv = HH(y);
MSE_conv = (1/(N*N))*abs(sum((abs(f_conv(:))-abs(f(:))).^2));
SNR_CONV = 10*log10(var(abs(f(:))) / MSE_conv)

Ph = exp(sqrt(-1)*(angle(f_conv(:))));               

X = abs((f_conv(:)));  %% Initial deðerler.
Sk = 0*abs((f_conv(:)));
Bk = 0*R(abs((f_conv(:))));

[cost,X,Sk, Bk] = lowRankSparseSar_Decomposition(R,Rt,X,Sk,Bk,y,f,Ph,H,HH,f_conv,parameters);


toc