%% Sub-Volume Registration based on Sub-Pixel Cross Correlation V1  
%% Image preprocessing
tic

% Set subpixel value: For better registration accuracies, it should be >=2
usfac = 1;

% Set the neighbor_window size if using the usfac > 2
neighbor_window_size = 0.3737;

% Binarize the scans to remove background information
binarize = 1; % 1 = Yes; 0 = No
Gaussian_Filter = 0;  % 1 = Yes; 0 = No  (gives bad results)

if binarize == 1
    AMPbin_Ref = imbinarize(ref1,0.15);
    AMPbin_Mov = imbinarize(mov1,0.15);

    RefBin = ref1.*AMPbin_Ref;
    MovBin = mov1.*AMPbin_Mov;
else
    RefBin = ref1;
    MovBin = mov1;
end

if Gaussian_Filter
    RefBin = imgaussfilt3(RefBin, 0.5);
    MovBin = imgaussfilt3(MovBin, 0.5);
end 

toc

%% Remove reflection of volume
RefBin_NoRefl = ReflectionRemover(RefBin);
MovBin_NoRefl = ReflectionRemover(MovBin);

%% cut out the black region in the Full-Volumes
RefBin = RefBin(1:400, :, 2:end);
MovBin = MovBin(1:400, :, 2:end);
RefBin_NoRefl = RefBin_NoRefl(1:400, :, 2:end);
MovBin_NoRefl = MovBin_NoRefl(1:400, :, 2:end);

%% MeanFree Volumes 
RefMean = mean(mean(mean(RefBin_NoRefl)));
MovMean = mean(mean(mean(MovBin_NoRefl)));
RefBin_NoRefl = RefBin_NoRefl - RefMean;
MovBin_NoRefl = MovBin_NoRefl - MovMean;

%% Cut off volumes Initialization 
% The minimum Cut-Off should be at least 50, 200, 200 in Y, X, and Z for
% better features retaining in  the registered Sub-Volumes.
X = 512;       % pixels in X-axis / number of A-scans
Y = 200;       % pixels in Y-axis / A-scan length
Z = 35 -1;     % pixels in Z-axis / number of B-scans 
               % with excluding the first frame from the registration process 

%% Registration
Registered_Volume = zeros(size(RefBin_NoRefl));
Registered_Volume = SubVolume_Registration_Bb_padding(RefBin_NoRefl, MovBin_NoRefl, MovBin, Y, X, Z, usfac, neighbor_window_size);      % New version with black boders padding                                  
% Registered_Volume = SubVolume_Registration(RefBin_NoRefl, MovBin_NoRefl, MovBin, Y, X, Z, usfac);               % Old version without padding 

    
%% Registered Volume Visualization 
figure;
for i = 1:size(RefBin,3) 
    imshowpair(RefBin(:,:,i), Registered_Volume(:,:,i));

    % Saving the ref, mov, registered frames for documentation
    % imwrite(Registered_Volume(:,:,i), ['Y:\Mitarbeiter\Aly\00_Lab_Book\Registration Results\2048X35_SBFingertip_100mm_10upsfac_0.3737\Registered_Volume_(200_512_35)_' num2str(i) '.png']);
    % imwrite(RefBin_NoRefl(:,:,i), ['Y:\Mitarbeiter\Aly\00_Lab_Book\Registration Results\2048X35_SBFingertip_100mm_10upsfac_0.3737\RefBin_'  num2str(i)  '.png']);
    % imwrite(MovBin_NoRefl(:,:,i), ['Y:\Mitarbeiter\Aly\00_Lab_Book\Registration Results\2048X35_SBFingertip_100mm_10upsfac_0.3737\MovBin_' num2str(i) '.png']);
    % fusedpair = imfuse(RefBin_NoRefl(:,:,i), Registered_Volume(:,:,i));
    % imwrite(fusedpair, ['Y:\Mitarbeiter\Aly\00_Lab_Book\Registration Results\2048X35_SBFingertip_100mm_10upsfac_0.3737\impair_' num2str(i) '.png']);

    colormap('gray')
    xlabel(i);
    pause(0.01)    
    % hold on;
    % title(['frame: ', num2str(i)])    
end

%% Save Video
% figure('Units','normalized','Position',[0 0 0.5 0.5],'visible', 'on')
% set(gcf, 'Color','w')
% 
% % open video file
% myVideo = VideoWriter([OCTPath,'\',ID,'_registered_video_48subvolumes_nothreshold']);
% % myVideo = VideoWriter([OCTPath,'\',ID,'_reference_video']);
% % myVideo = VideoWriter([OCTPath,'\',ID,'_moved_video']);
% % myVideo = VideoWriter([OCTPath,'\',ID,'_noregistration']);
% 
% myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
% open(myVideo)
% 
% for i = 1:(size(ref1,3)-1)
%     imshowpair(RefBin(:,:,i), Registered_Volume(:,:,i));
%     % imagesc(mov1(:,:,i))
%     % imshowpair(RefBin(:,:,i), MovBin(:,:,i));
%     axis off
%     colormap('gray')
%     pause(0.01)
% 
%     frame = getframe(gcf); %get frame
%     writeVideo(myVideo, frame);
% end
% 
% close(myVideo)





