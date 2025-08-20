translatedVolume1 = imtranslate(b_scans, [20.5,10.3,5]);
% translatedVolume2 = imtranslate(b_scans, [22.4,6.3,0.0]);

%% Auto-tune the best neighbor_window size for upsfac > 2
% tic
% previous_error = 1;
% best_neighborWindow = 1.500000; 
% iterations = 50;
% for i= 0:iterations
%     neighborWindow = rand(1) * 2;
%     [output] = dftregistration3D(fftn(b_scans(:,:,:)),fftn(translatedVolume1(:,:,:)),10,0, neighborWindow);
%     error = abs(abs(abs(output(3))-16.3) + abs(abs(output(4))-20.4) + abs(abs(output(5))-17.7));
%     if error < previous_error
%         best_neighborWindow = neighborWindow;
%         previous_error = error;
%     end
% end
% toc

%% The best after 10 iterations and gives almost 0.2 error at most is neighbor_window size = 0.3737
[output] = dftregistration3D(fftn(b_scans(:,:,:)),fftn(translatedVolume1(:,:,:)), 1, 0, 0, 0.3737);


%%
% [outputFrame] = dftregistration(fft2(b_scans(:,:,50)),fft2(translatedVolume2(:,:,50)),10);

