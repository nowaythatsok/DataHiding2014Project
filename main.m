%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Course:	Data Hiding 2014
% Project: 	Digital Image Forgery Detection Using JPEG Features and Local Noise Discrepancies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;clc;

start_time = cputime;

%we need JPEG image I suppose
file_name = 'DataHiding2014Project/img/tulipano90.jpg';
x       = imread(file_name);
img = rgb2gray(x);

quality_score = jpeg_quality_score(x);

alpha = 0.0;

if(quality_score < 2)
    alpha = 1;
elseif(quality_score < 6.9)
    alpha = (-0.0213 * quality_score) + 1.0469;
elseif(quality_score < 8.9)
    alpha = (-0.2980 * quality_score^2) + (4.2584 * quality_score) - 14.2952;
else
    alpha = 0;
end

alpha

% evaluate alpha*B + (1- alpha)*A


%% Time evaluation
stop_time = cputime;
fprintf('Execution time = %0.5f sec\n',abs( start_time - stop_time));