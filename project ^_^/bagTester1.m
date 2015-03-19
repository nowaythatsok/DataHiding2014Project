function [B]=bagTester1(length)

    A1=zeros(length,floor(length/2));
    A2=255*ones(length,floor(length/2));
    B=[A1,A2];
    %imshow(B);
end 



