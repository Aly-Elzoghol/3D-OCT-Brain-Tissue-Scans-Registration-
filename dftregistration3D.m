function [output, Greg, CCabs] = dftregistration3D(buf1ft,buf2ft,usfac,shiftedoutput, Subpixel_Zshift, neighborWindow)
% dftregistration3D Function for 3D-OCT scans registration based on the cross correlation matrix between the volumes.
% Inputs:
%        buf1ft: fft of the first volume
%        buf2ft: fft of the second volume
%        usfac: the upper sampling factor in case of the subpixel reg.
%        shiftedoutput: a selector for the user to choose if out registered
%                       image is required or not.
%        Subpixel_Zshift: a selector for the user to choose if Zshift will
%                         be updated in the case of the upsfac > 2 or not (normally user 
%                         should use integer shifts in the frame-wise dimension because frames in between will be blurry)
%        neighborWindow: The neighborhood window sizein the ccase of upsfac>2 
%                         in which the new CC matrix between the subpixel volumes will be calculated 

% Outputs:
%        output: A vector that contains
%               1. phase shift between the two correlated volumes.
%               2. error value in the calculted shift values.
%               3. shift values in X, Y, and Z
%        Greg: Registerd image using the esmated shift values.
%        CCabs: Cross correlation matrix between the ref and moved volumes

% Default usfac to 1
if exist('usfac')~=1, usfac=1; end

% Compute error for no pixel shift
if usfac == 0
    CCmax = sum(sum(sum(buf1ft.*conj(buf2ft))));
    rfzero = sum(abs(buf1ft(:)).^2);
    rgzero = sum(abs(buf2ft(:)).^2);
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero*rfzero); 
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));
    row_shift = 0;
    col_shift = 0;
    z_shift = 0;
    output=[error,diffphase,row_shift,col_shift,z_shift];


% Whole-pixel shift - Compute crosscorrelation by an IFFT and locate the
% peak
elseif usfac == 1
    [m,n,z]=size(buf1ft);
    
    % Single pixel registration
    CC = ifftn(buf1ft .* conj(buf2ft));
    % CC = xcorr3(ifftn(buf1ft), ifftn(buf2ft));

    CCabs = abs(CC);     
    [rloc, cloc, zloc] = ind2sub(size(CC), find(CC == max(CC(:))));

    % If subvolume can not get registered an error should be issued
    % Otherwise rloc, cloc and zloc is larger than the size of CC
    if size(rloc,1) > 1 || size(rloc,1) > size(CC,1)
        rloc = size(CC, 1);
        % warning('')
    end

    if size(cloc,1) > 1 || size(cloc,1) > size(CC,2)
        cloc = size(CC, 2);
    end
    if size(zloc,1) > 1 || size(zloc,1) > size(CC,3)
        zloc = size(CC, 3);
    end
    
    % Calculating the max based on one single peak
    CCmax = CC(rloc, cloc, zloc);

    % Calculate the max based on the max of the averaged xy, yz, xz cc
    % matrices as a refined max point (not worth!!)
    % max_xy = max(mean(CCabs(:,:,fix(size(CCabs,3)/2))));    
    % max_xz = max(mean(squeeze(CCabs(fix(size(CCabs,1)/2),:,:)),1));
    % max_yz = max(mean(squeeze(CCabs(:,fix(size(CCabs,2)/2),:)))); 
    % CCmax = max([max_xy max_xz max_yz]);

    rfzero = sum(abs(buf1ft(:)).^2)/(m*n*z);   
    rgzero = sum(abs(buf2ft(:)).^2)/(m*n*z);    
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero(1,1,1)*rfzero(1,1,1));
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));

    md2 = fix(m/2); 
    nd2 = fix(n/2);
    zd2 = fix(z/2);

    if rloc > md2
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end

    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end

    if zloc > zd2
        z_shift = zloc - z - 1;
    else
        z_shift = zloc - 1;
    end
    output=[error,diffphase,row_shift,col_shift, z_shift];
    
% Partial-pixel shift
else
    % First upsample by a factor of 2 to obtain initial estimate
    % Embed Fourier data in a 2x larger array
    [m,n,z]=size(buf1ft);
    mlarge=m*2;
    nlarge=n*2;
    zlarge=z*2;
    CC=zeros(mlarge,nlarge,zlarge);
    CC(m+1-fix(m/2):m+1+fix((m-1)/2),n+1-fix(n/2):n+1+fix((n-1)/2),z+1-fix(z/2):z+1+fix((z-1)/2)) = ...
        fftshift(buf1ft).*conj(fftshift(buf2ft));

    % Compute crosscorrelation and locate the peak 
    CC = ifftn(ifftshift(CC)); % Calculate cross-correlation  
    CCabs = abs(CC);
    [rloc, cloc, zloc] = ind2sub(size(CC), find(CC == max(CC(:))));
    
    % If subvolume can not get registered an error should be issued
    % Otherwise rloc, cloc and zloc is larger than the size of CC
    if size(rloc,1) > 1 || size(rloc,1) > size(CC,1)
        rloc = size(CC, 1);
        % warning('')
    end
    if size(cloc,1) > 1 || size(cloc,1) > size(CC,2)
        cloc = size(CC, 2);
    end
    if size(zloc,1) > 1 || size(zloc,1) > size(CC,3)
        zloc = size(CC, 3);
    end
    
    CCmax = CC(rloc, cloc, zloc);
    
    % Obtain shift in original pixel grid from the position of the
    % crosscorrelation peak 
    [m,n,z] = size(CC); md2 = fix(m/2); nd2 = fix(n/2); zd2 = fix(z/2);
    if rloc > md2 
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end
    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    if zloc > zd2
        z_shift = zloc - z - 1;
    else
        z_shift = zloc - 1;
    end    
    row_shift= row_shift/2;
    col_shift= col_shift/2;
    z_shift= z_shift/2;

    % If upsampling > 2, then refine estimate with matrix multiply DFT
    if usfac > 2
        % Initial shift estimate in upsampled grid
        row_shift = round(row_shift*usfac)/usfac; 
        col_shift = round(col_shift*usfac)/usfac;

        % in case of the sup-pixel with upsfac more than 2, Z axis(frame-wise is excluded)
        if Subpixel_Zshift
            z_shift = round(z_shift*usfac)/usfac;
        end

        % neighborWindow = 0.3737;  % The best value because in 2d we had
                                    % only one maxima within the 15 pixels window size which means
                                    % neighborWindow = 1.5, but here in 3d, we got multible maxima
                                    % within this window, then we need to make it much smaller to be 2
                                    % which means neighborWindow = 0.3737
        dftshift = fix(ceil(usfac*neighborWindow)/2); %% Center of output array at dftshift+1        

        % Matrix multiply DFT around the current shift estimate  
        CC = conj(dftups(buf2ft.*conj(buf1ft),ceil(usfac*neighborWindow),ceil(usfac*neighborWindow),ceil(usfac*neighborWindow),usfac,...
            dftshift-row_shift*usfac,dftshift-col_shift*usfac,dftshift-z_shift*usfac))/(md2*nd2*zd2*usfac^2);
        CCabs = abs(CC);

        % Locate maximum and map back to original pixel grid 
        [rloc, cloc, zloc] = ind2sub(size(CC), find(CC == max(CC(:))));

        % If subvolume can not get registered an error should be issued
        % Otherwise rloc, cloc and zloc is larger than the size of CC
        if size(rloc,1) > 1 || size(rloc,1) > size(CC,1)
            rloc = size(CC, 1);
            % warning('')
        end
        if size(cloc,1) > 1 || size(cloc,1) > size(CC,2)
            cloc = size(CC, 2);
        end
        if size(zloc,1) > 1 || size(zloc,1) > size(CC,3)
            zloc = size(CC, 3);
        end

        CCmax = CC(rloc, cloc, zloc);
        rg00 = dftups(buf1ft.*conj(buf1ft),1,1,1,usfac)/(md2*nd2*zd2*usfac^2);
        rf00 = dftups(buf2ft.*conj(buf2ft),1,1,1,usfac)/(md2*nd2*zd2*usfac^2);  
        rloc = rloc - dftshift - 1;
        cloc = cloc - dftshift - 1;
        zloc = zloc - dftshift - 1;

        row_shift = row_shift + rloc/usfac;
        col_shift = col_shift + cloc/usfac;    
        if Subpixel_Zshift
            z_shift = z_shift + zloc/usfac;            
        end

    % If upsampling = 2, no additional pixel shift refinement
    else    
        rg00 = sum(sum(sum( buf1ft.*conj(buf1ft) )))/m/n/z;
        rf00 = sum(sum(sum( buf2ft.*conj(buf2ft) )))/m/n/z;
    end
    error = 1.0 - CCmax.*conj(CCmax)/(rg00.*rf00);
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));

    % If its only one row or column the shift along that dimension has no
    % effect. We set to zero.
    if md2 == 1
        row_shift = 0;
    end
    if nd2 == 1
        col_shift = 0;
    end
    if zd2 == 1
        z_shift = 0;
    end    
    output=[error,diffphase,row_shift,col_shift,z_shift];
end  
if shiftedoutput == 1  % Check if the user wants the new registered image or not.
    % Compute registered version of buf2ft
    if (nargout > 1)&&(usfac > 0)
        [nr,nc,nz]=size(buf2ft);
        Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
        Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
        Nz = ifftshift(-fix(nz/2):ceil(nz/2)-1);
        [Nc,Nr,Nz] = meshgrid(Nc,Nr,Nz);
        Greg = buf2ft.*exp(1i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc-z_shift*Nz/nz));
        Greg = Greg*exp(1i*diffphase);

    elseif (nargout > 1)&&(usfac == 0)
        Greg = buf2ft*exp(1i*diffphase);
    end

else 
    Greg = zeros(size(buf1ft));
end
return

% Another implementation of the cross correlation matrix
function CC = xcorr3(buf1ft, buf2ft)
% CC = zeros(size(buf1ft));
tic
for i=1:size(buf1ft, 3)
    CC(:,:,i) = xcorr2_fft(buf1ft(:,:,i), buf2ft(:,:,i));
end
toc
return

function out = dftups(in, nor, noc, noz, usfac, roff, coff, zoff)
% function out=dftups(in,nor,noc,noz,usfac,roff,coff,zoff);
% Upsampled DFT by matrix multiplies, can compute an upsampled DFT in just
% a small region.
% usfac         Upsampling factor (default usfac = 1)
% [nor, noc, noz]     Number of pixels in the output upsampled DFT, in
%               units of upsampled pixels (default = size(in))
% roff, coff, zoff    Row, column, and depth offsets, allow to shift the output array to
%               a region of interest on the DFT (default = 0)
% Recieves DC in upper left corner, image center must be in (1,1,1) 
% Manuel Guizar - Dec 13, 2007
% Modified from dftus, by J.R. Fienup 7/31/06
% Modified from dftups_2D to 3D, by Aly Elzoghol 31/1/2024

% This code is intended to provide the same result as if the following
% operations were performed
%   - Embed the array "in" in an array that is usfac times larger in each
%     dimension. ifftshift to bring the center of the image to (1,1,1).
%   - Take the FFT of the larger array
%   - Extract an [nor, noc, noz] region of the result. Starting with the 
%     [roff+1 coff+1 coff+1] element.

% It achieves this result by computing the DFT in the output array without
% the need to zeropad. Much faster and memory efficient than the
% zero-padded FFT approach if [nor noc] are much smaller than [nr*usfac nc*usfac]

[nr,nc, nz]=size(in);
% Set defaults
if exist('roff')~=1, roff=0; end
if exist('coff')~=1, coff=0; end
if exist('zoff')~=1, zoff=0; end
if exist('usfac')~=1, usfac=1; end
if exist('noc')~=1, noc=nc; end
if exist('nor')~=1, nor=nr; end
if exist('noz')~=1, noz=nz; end

% % Compute kernels and obtain DFT by matrix products
kernc=exp((-1i*2*pi/(nc*usfac))*( ifftshift(0:nc-1).' - floor(nc/2) )*( (0:noc-1) - coff ));
%kernc=exp((-1i*2*pi/(nc*usfac))*( (0:noc-1).' - coff )*( ifftshift(0:nc-1) - floor(nc/2)));
kernr=exp((-1i*2*pi/(nr*usfac))*( (0:nor-1).' - roff )*( ifftshift(0:nr-1) - floor(nr/2)));
kernz=exp((-1i*2*pi/(nz*usfac))*( (0:noz-1).' - zoff )*( ifftshift(0:nz-1) - floor(nz/2)));

X = pagemtimes(kernr,in);
Y = pagemtimes(X,kernc);
di = size(Y);
Y = reshape(Y, di(2),di(3),di(1));
out = pagemtimes(Y,kernz');
% size(out)

return
