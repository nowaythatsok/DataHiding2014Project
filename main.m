%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Course:	Data Hiding 2014
% Project: 	Digital Image Forgery Detection Using JPEG Features and Local Noise Discrepancies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;clc;

start_time = cputime;

%we need JPEG image I suppose
file_name = 'DataHiding2014Project/img/lena.tif';
img       = imread(file_name);
[imH imW] = size(img);

%% Exclude edges. For preserving BAGs we only excluded edges within certain range.
% we detect edges using Sobel operator. edge = img * S
edges_img = edge(img,'sobel');
%imshow(edges_img);

% Then we defined whether a pixel is excluded
% First I need the edges gradient
[Gmag,Gdir] = imgradient(edges_img);
%figure; imshowpair(Gmag, Gdir, 'montage');
% Then we defined whether a pixel is excluded R=1 means excluded
R = zeros(imH, imW);
theta = 0.1;
% Define R using formula (1) of the paper
for i = 1:imH*imW
    if ((Gmag(i) > 0 && Gmag(i) < theta) || ((Gmag(i) > (pi/2 - theta)) && (Gmag(i) < (pi/2 + theta))) || ((Gmag(i) > (pi - theta)) && (Gmag(i) < pi + theta)))
        R(i) = 0;
    else
        R(i) = 1;
    end
end

%% Extract BAGs
d = zeros(imH, imW);
% absolute second-order difference formula (2), I've adjust index to the
% bound limit
for i = 2:imH-1
    for j = 1:imW
        d(i,j) = abs(2*img(i,j) - img(i-1,j) - img(i+1,j));
    end
end

%% print the images
%subplot(1,2,1); imshow(img);
%title(sprintf('Original %s', file_name));
%% Time evaluation
stop_time = cputime;
fprintf('Execution time = %0.5f sec\n',abs( start_time - stop_time));