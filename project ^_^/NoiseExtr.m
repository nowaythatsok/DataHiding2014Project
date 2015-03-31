clear all; close all;clc;
start_time = cputime;

 file_name = 'tulipano90_t';
 G = imread(file_name,'jpg');
%G = imread('tulipano10.jpg');
%% Noise estimation

%clustering with SLIC algo
% 1000 superpixels using a weighting factor of 1.5 
number_of_segment = 200;
[l, Am, C] = slic(G, number_of_segment, 1.5, 1, 'median');

[imH imW] = size(l);

%compute matrix with neighbours
adj_matrix = zeros(number_of_segment, number_of_segment);
for i=2:imH-1
    for j=2:imW-1
        for k=-1:1
            for x=-1:1
                if(l(i,j) ~= l(i+k,j+x) && adj_matrix(l(i,j), l(i+k,j+x)) ~= 0)
                    adj_matrix(l(i,j), l(i+k,j+x)) = 1;
                end
            end
        end
     end
end

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
            if(Z(i,j) == 0)
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

    counter = zeros(number_of_segment, 3, 5);
    
    %fill the array with the standard deviation
    for i=1:imH
        for j=1:imW
            if(Z(i,j) == 0)
                for k=1:3 %color
                    for c=1:5 %filter
                        F(l(i,j), k, c, 2) = F(l(i,j), k, c, 2) + (col_fil(k,c,i,j) - F(l(i,j), k, c, 1))^2;
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
                F(i, k, c, 2) = sqrt(F(i, k, c, 2))/counter(i, k, c);
              end
           end
        end
    end
    
%calculate F_mean that contains mean values for each feature calculate before 
F_mean = zeros(3, 5, 2);

for i=1:number_of_segment
   for k=1:3 %color
      for c=1:5 %filter
          F_mean(k, c, 1) = F_mean(k, c, 1) + F(i, k, c, 1);
          F_mean(k, c, 2) = F_mean(k, c, 2) + F(i, k, c, 2);
      end
   end
end

for k=1:3 %color
   for c=1:5 %filter
        F_mean(k, c, 1) = F_mean(k, c, 1)/number_of_segment;
        F_mean(k, c, 2) = F_mean(k, c, 2)/number_of_segment;
   end
end

%now I need the most distant vector from mean in the set (using euclidian distance)
max_distance = -1;
F_max=zeros(3, 5, 2);


for i=1:number_of_segment
    distance = my_euclidian_distance(squeeze(F(i,:,:,:)), F_mean);
    if(distance > max_distance)
        max_distance = distance;
        F_max = squeeze(F(i,:,:,:));
    end
end

% now we evaluate the weight factor of everi segmente calculated by SLIC
w = zeros(number_of_segment, 2); 
% 1=weigth alpha(original), 
% 2=weight beta(forged), 
for i=1:number_of_segment
    w(i,1) = my_euclidian_distance(squeeze(F(i,:,:,:)), F_mean);
    w(i,2) = my_euclidian_distance(squeeze(F(i,:,:,:)), F_max);
end


K = 0.1; %constant K
GW = zeros(number_of_segment+2, number_of_segment);
C = zeros(number_of_segment+2, number_of_segment);
L = zeros(imH, imW);

%init label  
for i=1:number_of_segment
  if(w(i,1) < w(i,2)) %if original is greater then tampered init L to 0
  %if(w(i,1) < (w(i,2)-128))
          C(number_of_segment+1, i) = 1;
  else
          C(number_of_segment+2, i) = 1;
  end   
end

%C([number_of_segment+1, number_of_segment+2], :)



E=1;
while(E ~= 0) %if the function return 0 there is no possible improvement of E
    [GW, C, E] = mincut(GW, C, adj_matrix, number_of_segment, w);
end

% create label matrix L
for i=1:imH
    for j=1:imW
        L(i,j) = C(number_of_segment+2, l(i,j)); %if tbeta = 1 label 1 else label 0
    end
end

% create block indicator A
A = zeros(imH, imW);
for i=1:8:imH-8
    for j=1:8:imW-8
        temp = 0;
        for k=0:7
            for z=0:7
                temp = temp + L(i+k,j+z);
            end
        end
        temp = temp / 64;
        for k=0:7
            for z=0:7
                A(i+k,j+z) = temp;
            end
        end
    end
end

for i=1:imH
    for j=1:imW
        if(A(i,j) >= 0.5)
        %if(L(i,j) >= 0.5)
            G(i,j,1) = 255;
            G(i,j,2) = 255;
            G(i,j,3) = 255;
        end
    end
end

imshow(G);

%% print the images
%subplot(1,2,1); imshow(img);
%title(sprintf('Original %s', file_name));
%% Time evaluation
stop_time = cputime;
fprintf('Execution time = %0.5f sec\n',abs( start_time - stop_time));