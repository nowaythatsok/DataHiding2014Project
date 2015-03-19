clear all; close all;clc;

%parameters:
theta=pi()/180*5;
debug_on=0;
plot_switch=1;

if debug_on == 0 
    %file_name = 'lena.tif';
    %img  = imread(file_name);
    %imwrite(uint8(img),'lena.jpg','jpg','Quality',10);
    
    
    %file_name = 'lena.jpg'
    %G = imread(file_name,'jpg');
    
    file_name = 'tulipano10.jpg';
    %file_name = 'tampered.jpg';
    img = imread(file_name,'jpg');
    G   = rgb2gray(img);
    %G=G(700:1300,700:1300);
   
    
else

    %call test image - grayscale
    %G=bagTester1(12);
    G=bagTester2(60);
    %G=bagTester3(600);
end

G=int32(G);

%get dimensions
[Xdim,Ydim]=size(G);

%use sobel operator:
%E2=edge(G,'sobel');
[E,A]=Sobel(G);
if plot_switch==1
    subplot(2,2,1);imshow(uint8(G));
    subplot(2,2,2);imshow(E);
    %subplot(2,2,3);imshow(E2);
end

%create exclusion function, now 0 is excluded
R=zeros(Xdim, Ydim);
for i=1:Xdim
    for j=1:Ydim
        if ( (A(i,j)<theta) || (A(i,j)> (pi()-theta)) || (A(i,j)>(pi()/2-theta) && A(i,j)<(pi()/2+theta)))
            R(i,j)=1;
        end 
    end
end

%show likely result
if plot_switch==1
    E3=E.*uint8(R);
    subplot(2,2,3);imshow(E3);
end

%extract BAGs
%second ordex difference:   d_x in the X direction (1st coordinate) 
%                           d_y
%of range 0-255
d_x=zeros(Xdim,Ydim);
d_y=zeros(Xdim,Ydim);

for i=2:Xdim-1
    for j=2:Ydim-1
        if R(i,j)==0
            d_x(i,j)=0;         %actual exclusion
            d_y(i,j)=0;
        else
            d_x(i,j)=abs(2*G(i,j)-G(i-1,j)-G(i+1,j));
            d_y(i,j)=abs(2*G(i,j)-G(i,j-1)-G(i,j+1));
            if d_x(i,j)>25
                d_x(i,j)=0;     %here too
            end
            if d_y(i,j)>25
                d_y(i,j)=0;
            end
        end
    end
end


%accumulation
%of range 0-255*33
a_x=zeros(Xdim,Ydim);
a_y=zeros(Xdim,Ydim);
for i=18:Xdim-17
    for j=18:Ydim-17
        for k=-16:16
            a_x(i,j)=a_x(i,j)+d_x(i,j+k);
            a_y(i,j)=a_y(i,j)+d_y(i+k,j);
        end
    end
end


%median reduction filtering
%of range -255-255*32
a_rx=zeros(Xdim,Ydim);
a_ry=zeros(Xdim,Ydim);
for i=18:Xdim-17
    for j=18:Ydim-17
        a_rx(i,j)=a_x(i,j)-median(a_x((i-16:i+16),j));
        a_ry(i,j)=a_y(i,j)-median(a_y(i,(j-16:j+16)));
    end
end

%median filtering
%of range same
b=zeros(Xdim,Ydim);
for i=18:Xdim-17
    for j=18:Ydim-17
        b(i,j)=median([a_rx(i-16,j) a_rx(i-8,j) a_rx(i,j) a_rx(i+8,j) a_rx(i+16,j)])...
            +  median([a_ry(i,j-16) a_ry(i,j-8) a_ry(i,j) a_ry(i,j+8) a_ry(i,j+16)]); 
    end
end

if plot_switch==1
    subplot(2,2,4); imshow(uint8(b*255/max(max(b,[],1),[],2)));
end

if plot_switch==1
    figure;
    subplot(2,2,1); imshow(E);
    subplot(2,2,2); imshow(d_x);
    subplot(2,2,3); imshow(a_x);
    subplot(2,2,4); imshow(a_rx);
end


    figure;
    imshow(uint8(G)+uint8(b));


%b=uint8(b);

%get starting point for feature extracion
%collapse cols and rows
Cx=zeros(Xdim,1);%collapse for 1st coordinate
for i=1:Ydim
    Cx=Cx+b(1:Xdim,i);
end
cx8=zeros(8,1); %1st coord to 8
for i=1:(floor(Xdim/8)-1)
    cx8=cx8+Cx(8*i:8*i+7);
end

%we start from the maximal gradient
gx=zeros(8,1);
for i=1:8
    gx(i)=cx8(i)-cx8(mod(i,8)+1);
end
[~,X]=max(gx);

if plot_switch==1
    figure; 
    subplot(2,2,1);plot(cx8);
    subplot(2,2,3);plot(gx);
end
    
Cy=zeros(1,Ydim);
for i=1:Xdim
    Cy=Cy+b(i, 1:Ydim);
end
cy8=zeros(1,8);
for i=1:(floor(Ydim/8)-1)
    cy8=cy8+Cy(8*i:8*i+7);
end
gy=zeros(1,8);
for i=1:8
    gy(i)=cy8(i)-cy8(mod(i,8)+1);
end

if plot_switch==1
    subplot(2,2,2);plot(cy8);
    subplot(2,2,4);plot(gy);
end
[~,Y]=max(gy);

translate_x=X-1;
translate_y=Y-1;

%plot check MTRX
Ch=ones(Xdim, Ydim);
for i=1:(floor(Xdim/8)-2)
    for j=1:(floor(Ydim/8)-2)
        Ch(8*i+translate_x,8*j+translate_y)=0;
    end
end

if plot_switch==1
    figure;
    imshow(Ch);
end

%BAG feature EXTRACTION
B=zeros(floor(Xdim/8),floor(Ydim/8));
line_sum=zeros(1,8);
column_sum=zeros(8,1);
for i=1:(floor(Xdim/8)-2)
    for j=1:(floor(Ydim/8)-2)
        %Accumulate lines and columns
        line_sum=zeros(1,8);
        column_sum=zeros(8,1);
        for k=1:6
            
            line_sum=line_sum+b(8*i+k+translate_x,((8*j:8*j+7)+translate_y));
            
            column_sum=column_sum+b(((8*i:8*i+7)+translate_x),8*j+k+translate_y);
        end
        B(i,j)=max(line_sum(2:7))+max(column_sum(2:7))...
            -min([line_sum(1), line_sum(8)])-min([column_sum(1), column_sum(8)]);
    end
end
Bu=uint8((255*int64(B))/max(max(B,[],1),[],2));
Ba=uint8((255*int64(abs(B)))/max(max(abs(B),[],1),[],2));
%B(1,1)=20*255*
%B
figure;
subplot(2,2,2);imshow(Bu);
subplot(2,2,4);imshow(Ba);
subplot(2,2,1);imshow(B);





        