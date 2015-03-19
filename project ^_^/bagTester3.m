function [B]=bagTester3(size)

    B=zeros(size);
    for i=1:size
        for j=1:i
            B(i,j)=255;
        end 
    end        
    %imshow(B);
end 