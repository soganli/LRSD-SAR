function x = Rt(X,parameters)

patch_size_x = parameters.patch_sizex;
patch_size_y = parameters.patch_sizey;

Weight = parameters.patchWeight;
indexes = parameters.patchIndexes;
N = parameters.N;

IMout = zeros(N,N);
[rows,cols] = ind2sub([N N]-[patch_size_x patch_size_y]+1,indexes);
count = 1;
for j  = 1:length(cols)
    col = cols(j); row = rows(j);        
    IMout(row:row+patch_size_x-1,col:col+patch_size_y-1)=IMout(row:row+patch_size_x-1,col:col+patch_size_y-1)+reshape(X(:,count),[patch_size_x,patch_size_y]);
    count = count+1;
end;

x_image = IMout ./ Weight;

x = x_image(:);

end