%% Accuracy test
% This code is applying random manual shifts to the ref volume and 
% calculating the mean error within a fixed number of iterations to test
% the accuracy while using 1.5, and 0.3737 neighbor_window size.

tic

iters = 10;
AverageError = 0;
Xerror = 0;
Yerror = 0;
Zerror = 0; 

for i=1:iters
    dx= round((rand(1) * 20),1);
    dy= round((rand(1) * 10),1);
    dz= ceil(rand(1) * 10);
    translatedVolume = imtranslate(ref, [dx,dy,dz]);
    [output] = dftregistration3D(fftn(ref(:,:,:)),fftn(translatedVolume(:,:,:)), 1, 0, 0, 0.3737);
    Xerror = Xerror + abs(abs(output(4))-dx);       
    Yerror = Yerror + abs(abs(output(3))-dy);
    Zerror = Zerror + abs(abs(output(5))-dz);
end
Xerror = Xerror / iters
Yerror = Yerror / iters
Zerror = Zerror / iters

AverageError = Xerror + Yerror + Zerror


toc