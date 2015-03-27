function [distance] = my_euclidian_distance(x, y)
%EUCLIDIAN_DISTANCE in a 3x5x2 array
%   Detailed explanation goes here
    dist = x - y;
    square_dist = dist.^2;
    distance = sqrt(sum(sum(sum(square_dist))));
end

