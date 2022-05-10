function [X] = R(x,parameters)

if size(x,2) == 1
N1 = sqrt(length(x)); 
N2 = N1;
else
N1 = size(x,1);
N2 = size(x,2);
end
    
[X,~] = my_im2col(reshape(x,N1,N2),parameters);  % We obtain the patches for dictionary process.  bb is 8: For 32x32 image there will be 625 patches.


end