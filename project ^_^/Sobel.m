function [S,A]=Sobel(G)
    %define directional operators
    Oy=[-1 -2 -1; 0 0 0; 1 2 1];
    Ox=Oy';
    
    %use them
    Gx = conv2(double(G),Ox, 'same');
    Gy = conv2(double(G),Oy,'same');
    S = uint8(sqrt(Gx.*Gx + Gy.*Gy));
    A=abs(angle(Gx+Gy*i));
end
    