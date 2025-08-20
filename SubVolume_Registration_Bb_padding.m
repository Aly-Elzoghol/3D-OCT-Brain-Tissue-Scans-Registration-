function Com_Registered = SubVolume_Registration_Bb_padding(RefBin_NoRefl, MovBin_NoRefl, MovBin, Y, X, Z, usfac, neighbor_window_size)
%% Cut off volumes Initialization 
% Initialize the starting and ending points in each direction X, Y, and Z
X_s = 1; 
X_f = X;
Y_s = 1;
Y_f = Y;
Z_s = 1;
Z_f = Z;


%% Initialization
tic

% Initialize the the Sub-Volume cell size
[y, x, z] = size(RefBin_NoRefl);  % The whole volume size before Cutting off
C =  fix(x / X);    % no.of Sub-Volumes in X
R = fix(z / Z);     % no.of Sub-Volumes in Y
D = fix(y / Y);     % no.of Sub-Volumes in Z

% Initialization of the Sub-Volumes cells
output = cell(D, C, R);             % Shift values for each Sub-Volumes
% Greg = cell(D, C, R);             % Registered Sub-Volumes cell before ifft process
% Moved_Image = cell(D, C, R);      % Moved Sub-Volumes
% Ref_Image = cell(D, C, R);        % Ref Sub-Volumes
Img_Registered = cell(D, C, R);     % Final Img_Registered Sub-Volume cell 
CCabs = cell(D, C, R);              % Sub-Volumes cross correlation cell 

for r = 1: R
    for c = 1: C
        for d = 1: D
            output{d, c, r} = zeros(5, 1);
            % Greg{d, c, r}(:,:,:) = zeros(Y, X, Z);
            % Moved_Image{d, c, r}(:,:,:) = zeros(Y, X, Z);
            % Ref_Image{d, c, r}(:,:,:) = zeros(Y, X, Z);
            Img_Registered{d, c, r}(:,:,:) = zeros(Y, X, Z);
            if usfac == 1
                CCabs{d, c, r}(:,:,:) = zeros(Y, X, Z);
            end 
            if usfac == 2
                CCabs{d, c, r}(:,:,:) = zeros(2*Y, 2*X, 2*Z);
            end  
        end
    end
end

% Sub-Volumes Initializations with X, Y, Z:
Ref_XY = zeros(Y, X, Z);
Move_XY= zeros(Y, X, Z);

% Padding of the moved volume to form the black borders to avoid the
% negative shifts outside the volume borders.
% The padding can be done using higher Max_Bb depending on the
% requirements.
Max_Bb = max(X/2,Y/2);
MovBin = padarray(MovBin, [Max_Bb Max_Bb Max_Bb], 0, 'both');

toc

%% Putting the whole volumes on the GPU for faster calculations (it saves about 10 sec, but look at the problem of resizeing)
% RefBin_NoRefl = gpuArray(RefBin_NoRefl);
% MovBin_NoRefl = gpuArray(MovBin_NoRefl);

%% Registration
tic

for r = 1: R
    for c = 1: C
        for d = 1: D
            % Cut out the REF Sub-Volume based on the starting and ending points each time there is new Sub-Volume
            % Cut out the MOV Sub-Volume based on the starting and ending points each time there is new Sub-Volume
            Ref_XY = fftn(RefBin_NoRefl(Y_s:Y_f, X_s:X_f, Z_s:Z_f));     % Getting the fftn for the REF and MOV ones (Note that only fftn gives the correct results)
            Move_XY = fftn(MovBin_NoRefl(Y_s:Y_f, X_s:X_f, Z_s:Z_f));
            % Moved_Image{d, c, r} = MOV;  % Assigning the REF and MOV ones in the Sub-Volumes cell 
            % Ref_Image{d, c, r} = REF;

            [output{d, c, r}, ~ , CCabs{d, c, r}] = dftregistration3D(Ref_XY, Move_XY, usfac, 0, 0, neighbor_window_size);    % Perform the 3D registration for the Cut-Off Sub-Volumes and get the amount of 
                                                                                                        % shifts in  X, Y, and Z
            
            % Additional step here might be useful (neglect or replace them with the ones from the previous frame, the unconsistent values from the output calculated shifts)  
            %
            %
            %

            
            % Calculate the new starting and ending points using the amount of shifts from the registration process
            % Note that using usfac = 2 lead to float values of the shifts
            % Calculate the new starting and ending points using the amount of shifts from the registration process
            % Note that using usfac = 2 lead to float values of the shifts
            Y_s_new = floor(Y_s - output{d, c, r}(3)) + Max_Bb;
            Y_f_new = floor(Y_f - output{d, c, r}(3)) + Max_Bb;
            X_s_new = floor(X_s - output{d, c, r}(4)) + Max_Bb;
            X_f_new = floor(X_f - output{d, c, r}(4)) + Max_Bb;
            Z_s_new = floor(Z_s - output{d, c, r}(5)) + Max_Bb;
            Z_f_new = floor(Z_f - output{d, c, r}(5)) + Max_Bb;

            % Check if the new starting and ending points exceed the padded MovBin limits.
            if Y_s_new < 1;                     Y_s_new = 1;                    end
            if Y_f_new > size(MovBin, 1);       Y_f_new = size(MovBin, 1);      end

            if X_s_new < 1;                     X_s_new = 1;                    end
            if X_f_new > size(MovBin, 2);       X_f_new = size(MovBin, 2);      end 

            if Z_s_new < 1;                     Z_s_new = 1;                    end
            if Z_f_new > size(MovBin, 3);       Z_f_new = size(MovBin, 3);      end 

            % Getting the registered back Sub-Volumes from the Bb Moved Full-Volume
            % based on the new calculated starting and ending points after
            % updating the amounts of shifts from the main starting and
            % ending points.
            Img_Registered{d, c, r} = MovBin(Y_s_new:Y_f_new, X_s_new:X_f_new, Z_s_new:Z_f_new);        

            % Updating the Sub-Volumes starting and ending points each time
            % there is a new Sub-Volume, and also check if the new points
            % exceed the main Volume limits.
            % Update Y and check its limits
            Y_s = Y_s + Y - 1;
            Y_f =  Y_f + Y - 1;    
            if Y_s < 1;             Y_s = 1;            end
            if Y_s > (y - Y + 1);   Y_s = y - Y + 1;    end            
            if Y_f > y;             Y_f = y;            end
        end
        % Update X and check its limits conidering intitializing Y again
        Y_s = 1;
        Y_f = Y;
        X_s = X_s + X - 1;
        X_f =  X_f + X - 1;
        if X_s < 1;             X_s = 1;            end
        if X_s > (x - X + 1);   X_s = x - X + 1;    end        
        if X_f > x;             X_f = x;            end
    end
    % Update Z and check its limits conidering intitializing Y and X again
    Y_s = 1;
    Y_f = Y;
    X_s = 1;
    X_f = X;
    Z_s = Z_s + Z - 1;
    Z_f =  Z_f + Z - 1;    
    if Z_s < 1;                 Z_s = 1;            end
    if Z_s > z - Z + 1;         Z_s = z - Z + 1;    end    
    if Z_f > z;                 Z_f = z;            end
end

toc


%% Stitching the SubVolumes in the moved image into the new Registered Image 
% Combining all the Sub-Volumes into one big Volume with the size of the
% main Volume
tic
   Com_Registered = zeros(size(RefBin_NoRefl)); 
cat1 = cell(C, R);
cat2 = cell(1, R);
for j = 1:R
    for i = 1:C
        cat1{i, j}= cat(1, Img_Registered{:, i, j});
    end
    cat2{1, j} = cat(2, cat1{:,j});
end
Com_Registered = cat(3, cat2{1, :});
toc

return
