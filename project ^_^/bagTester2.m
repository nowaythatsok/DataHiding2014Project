function [B]=bagTester2(piece)

    s=8;
    A1=zeros(s);
    A2=255*ones(s);
    C1=A1;
    C2=A2;
    for i=1:(piece)
        C1=[C1,A2,A1];
        C2=[C2,A1,A2];
    end 
    B=C1;
    for i=1:(piece)
        B=[B;C2;C1];
    end 
%imshow(B);
end 