clear all; close all;clc;
start_time = cputime;

file_name = 'tulipano10';
G = imread(file_name,'jpg');
[imH imW] = size(G);

%% Noise estimation

%clustering with SLIC algo
% 1000 superpixels using a weighting factor of 1.5 
number_of_segment = 1000;
[l, Am, C] = slic(G, number_of_segment, 1.5, 1, 'median');
%imshow(drawregionboundaries(l, G, [255 255 255]))

%rgb2gray converts RGB values to grayscale values by forming 
%a weighted sum of the R, G, and B components:
%0.2989 * R + 0.5870 * G + 0.1140 * B 
I = rgb2gray(G);
%sobel operation to get edges
%G=int32(G);
E2=edge(I,'sobel');
[imH imW] = size(I);
%[E,A]=Sobel(G);

H = zeros(imH, imW);
for i = 1:imH*imW
    if(E2(i) ~= 0)
        H(i) = 1;
    end
end

V = ones(imH, imW);
Z = xor(H, V);

%% creation of filters
im_red = G(:, :, 1);
im_green = G(:, :, 2);
im_blue = G(:, :, 3);

Hr1 = medfilt2(im_red);
Hg1 = medfilt2(im_green);
Hb1 = medfilt2(im_blue);
gaussian_filter = fspecial('gaussian',[3 3], 0.5);
Hr2 = imfilter(im_red, gaussian_filter, 'replicate');
Hg2 = imfilter(im_green, gaussian_filter, 'replicate');
Hb2 = imfilter(im_blue, gaussian_filter, 'replicate');
avarage_filter = fspecial('average');
Hr3 = imfilter(im_red, avarage_filter, 'replicate');
Hg3 = imfilter(im_green, avarage_filter, 'replicate');
Hb3 = imfilter(im_blue, avarage_filter, 'replicate');
Hr4 = wiener2(im_red,[3 3]);
Hg4 = wiener2(im_green,[3 3]);
Hb4 = wiener2(im_blue,[3 3]);
Hr5 = wiener2(im_red,[5 5]);
Hg5 = wiener2(im_green,[5 5]);
Hb5 = wiener2(im_blue,[5 5]);

%% extraction of feature for every combination color/filter
% color_filter(color, filter)
col_fil = zeros(3, 5, imH, imW);
col_fil(1,1,:,:) = im_red - Hr1;
col_fil(2,1,:,:) = im_green - Hg1;
col_fil(3,1,:,:) = im_blue - Hb1;
col_fil(1,2,:,:) = im_red - Hr2;
col_fil(2,2,:,:) = im_green - Hg2;
col_fil(3,2,:,:) = im_blue - Hb2;
col_fil(1,3,:,:) = im_red - Hr3;
col_fil(2,3,:,:) = im_green - Hg3;
col_fil(3,3,:,:) = im_blue - Hb3;
col_fil(1,4,:,:) = im_red - Hr4;
col_fil(2,4,:,:) = im_green - Hg4;
col_fil(3,4,:,:) = im_blue - Hb4;
col_fil(1,5,:,:) = im_red - Hr5;
col_fil(2,5,:,:) = im_green - Hg5;
col_fil(3,5,:,:) = im_blue - Hb5;

% PARAMETERS:
% 1) number_of_segment: number of segment in the image
% 2) red=1, green=2, blue=3
% 3) median=1, gaussian=2, average=3, weiner3=4, weiner5=5
% 4) mean=1, standard deviation=2
F = zeros(number_of_segment, 3, 5, 2);
counter = zeros(number_of_segment, 3, 5);

%fill the multidim array with the mean
    for i=1:imH
        for j=1:imW
            if(Z(i,j))
                for k=1:3 %color
                    for c=1:5 %filter
                        F(l(i,j), k, c, 1) = F(l(i,j), k, c, 1) + col_fil(k,c,i,j);
                        counter(l(i,j), k, c) = counter(l(i,j), k, c) + 1;
                    end
                end
            end
        end
    end
    
    for i=1:number_of_segment
        for k=1:3 %color
           for c=1:5 %filter
              if(counter(i, k, c) ~= 0)
                F(i, k, c, 1) = F(i, k, c, 1)/counter(i, k, c);
              end
           end
        end
    end

    %fill the array with the standard deviation
    for i=1:imH
        for j=1:imW
            if(Z(i,j))
                for k=1:3 %color
                    for c=1:5 %filter
                        F(l(i,j), k, c, 2) = F(l(i,j), k, c, 2) + (col_fil(k,c,i,j) - F(l(i,j), k, c, 1))^2;
                    end
                end
            end
        end
    end
    
    for i=1:number_of_segment
        for k=1:3 %color
           for c=1:5 %filter
              if(counter(i, k, c) ~= 0)
                F(i, k, c, 2) = sqrt(F(i, k, c, 2))/counter(i, k, c);
              end
           end
        end
    end
    
% figure;
% subplot(2,2,1); imshow(Fr5);
% subplot(2,2,2); imshow(Fg5);
% subplot(2,2,3); imshow(Fb5);
% subplot(2,2,4); imshow(Hg5);


%% print the images
%subplot(1,2,1); imshow(img);
%title(sprintf('Original %s', file_name));
%% Time evaluation
stop_time = cputime;
fprintf('Execution time = %0.5f sec\n',abs( start_time - stop_time));