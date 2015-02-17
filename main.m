%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Course:	Data Hiding 2014
% Project: 	Digital Image Forgery Detection Using JPEG Features and Local Noise Discrepancies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;clc;

start_time = cputime;

%we need JPEG image I suppose
file_name = 'DataHiding2014Project/img/tulipano10.jpg';
x       = imread(file_name);
img = rgb2gray(x);

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
theta = 1;
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

a = zeros(imH, imW);
ar = zeros(imH, imW);
b = zeros(imH, imW);

for i = 1:imH
    for j = 17:imW
        for k = 1:16
            if(R(i,j) == 0 && d(i,j) < 0.1)
                a(i,j) = a(i,j) + d(i,k);
            end
        end
    end
end

med = medfilt2(a, [16 1]);

for i = 1:imH
    for j = 1:imW
        ar(i,j) = a(i,j) - med(i,j);
    end
end

%imshow(ar);

for i = 17:imH-17
    for j = 17:imW-17
        temph = medfilt1(ar([i-16 i-8 i i+8 i+16],j),5);
        tempv = medfilt1(ar(i,[j-16 j-8 j j+8 j+16]),5);
        b(i,j) = max(temph) + max(tempv);
    end
end

%imshow(b);

%% Noise estimation
H = zeros(imH, imW);
for i = 1:imH*imW
    if(edges_img(i) ~= 0)
        H(i) = 1;
    end
end

V = ones(imH, imW);
Z = xor(H, V);

%% print the images
%subplot(1,2,1); imshow(img);
%title(sprintf('Original %s', file_name));
%% Time evaluation
stop_time = cputime;
fprintf('Execution time = %0.5f sec\n',abs( start_time - stop_time));