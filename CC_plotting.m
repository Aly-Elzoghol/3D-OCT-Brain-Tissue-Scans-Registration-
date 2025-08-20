
%% CC plotting Function

% In case of the subvolume CC list
% load('CCabs.mat');

% X = zeros(size(CCabs{1,3}, 2));
% Y = zeros(1,size(CCabs{1,3}, 1));
% 
% for i=1:size(CCabs{1,3}, 3)
%     CC = CCabs{1,4}(:,:,i);
% %     im = mat2gray(CC);
% %     im = ind2rgb(im);
% %     imwrite(im, ['Y:\Mitarbeiter\Aly\00_Lab_Book\Registration Results\Cross Correlation Sub-Volumes\CCabs_(200_200)_' num2str(i) '.png']);    
%     Y(i) = max(CC(:));
%     X(i) = i;
% %     [X(i), Y(i)] = ind2sub(size(CC), find(CC == max(CC(:))));
% end
% 
% figure;
% plot(Y,'*-')

    
% In case of the only one CC matrix for the whole volume
figure(5)
for i=1:size(CCabs,3)
    surf(CCabs(:,:,i));
    hold on
    grid on
end
