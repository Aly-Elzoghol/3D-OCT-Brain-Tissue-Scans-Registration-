# 3D-OCT-Brain-Tissue-Scans-Registration-
dftregistration3D Function for 3D-OCT scans registration based on the cross correlation matrix between the volumes.

# Approach
3D registration based on Sub-Pixel DFT matrix multiplication using black borders padding Sub-Volumes stiching function.
1. SubPixel DFT matrix multiplication: Use upsampling factor to enable half pixel shift values for upsfac=2, or more accurate tiny shift values using DFT matrix mutiplication within a predefined neighboring window size using the first shift estimate from usfac=2 step. 
2. Black borders Padding: By padding the moved volume in x, y, and z with a padding size equal to half of the sub-volume size to avoid the negative shift values and maintain the registration process being as  accurate as possible.
3. Sub-Volumes stching: After the registration is done, the subvolumes are needed to be stiched together to form the main registered volume with the size of the original reference volume.

# Oct-scans charachterisitics
## PVAPhantom
    - 1312CWL 
    - 101BW
    - 100nblank
    - 40phase
    - 50offset
    - Nopumpedtable
    - 1000x75
    - shaking_100µm
    - XYZ Movement

 ## Fingertip
    - 1312CWL
    - 101BW
    - 200nblank
    - 30phase
    - 90offset
    - Nopumpedtable
    - 2048x35


# Function Parameters
## Usfac: For more accurate registration results, upsampling factor need to be greater than 2 for smaller Sub-Volumes, however in the case that the hardware system is doing upsampling. Then software upsampling is no longer required.
## Cut-off sub-volume: The minimum Cut-Off should be at least 50, 200, 200 in Y, X, and Z for better features retaining in the    registered Sub-Volumes.
    - X = 200; % pixels in X-axis / number of A-scans
    - Y = 200; % pixels in Y-axis / A-scan length
    - Z = 75;   % pixels in Z-axis / number of B-scans  (Higher frame size leads to better accuracies)
## neighbor_window_size: The size of the window beyond which the first estimate shift values are updated in case of DFT matrix multiplication.
    - 1.5 in 2D , 15 pixels window
    - 0.3737 in 3D , 4 pixels window which is much smaller than in 2D to avoid the multible optimal points in larger window sizes.
## Subpixel_Zshift: Supixels in Z axis lead to to blurry images.


# Results
## XYZ Movements 1000 frames 
<img width="1563" height="625" alt="image" src="https://github.com/user-attachments/assets/224a24dd-0b45-4cdb-bb71-920a9e43822c" />


## Fingertip
<img width="1697" height="331" alt="image" src="https://github.com/user-attachments/assets/df406909-7832-401e-8e28-29fcaa616ba0" />


# Software Accuracy
Applying random manual shifts to the ref volume and calculating the mean error within iterations number = 10 for four different scans:
<img width="3557" height="1067" alt="image" src="https://github.com/user-attachments/assets/a7147e9e-e0d0-4a5f-86f7-7b0d69b94b65" />






