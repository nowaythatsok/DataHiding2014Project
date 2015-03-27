function [newG, newC, E] = mincut(GW, C, adj_matrix, number_of_segment, w)
%mincut
%   Detailed explanation goes here
    %build the graph
for i=1:number_of_segment
   for j=1:number_of_segment
       if(adj_matrix(i,j) == 1 && GW(l(i,j), l(i+k,j+x)) == 0 && C(l(i,j)) ~= C(l(i+k,j+x)))
           GW(l(i,j), l(i+k,j+x)) = K;
           C(l(i,j), l(i+k,j+x)) = 1;
       end
   end
end

for i=number_of_segment+1:number_of_segment+2
  for j=1:number_of_segment
      %evaluate neighbours values
      temp = 0;
      for k=1:number_of_segment
          temp = temp + GW(j, k);
      end
      GW(i,j) = w(j, i-number_of_segment) + temp; 
  end  
end
    

%evaluate and update the cut
    oldE = sum(sum(xor(GW, C)));
    min = 10^99;
    [h, w] = size(GW);
    for i=1:h-2
      for j=1:w
         if(adj_matrix(i,j) == 1)
             all_zero = GW(i, h-1) + GW(j, h-1); 
             all_one = GW(i, h) + GW(j, h);
             zero_one = GW(i, h-1) + GW(j, h) + GW(i,j);
             one_zero = GW(i, h) + GW(j, h-1) + GW(j,i);
             if(all_zero == min(all_zero, all_one, zero_one, one_zero))
                C(h-1, i) = 1;
                C(h-1, j) = 1;
                C(h, i) = 0;
                C(h, j) = 0;
                C(i, j) = 0;
                C(j, i) = 0;
             elseif(all_one == min(all_zero, all_one, zero_one, one_zero))
                C(h-1, i) = 0;
                C(h-1, j) = 0;
                C(h, i) = 1;
                C(h, j) = 1;
                C(i, j) = 0;
                C(j, i) = 0;
             elseif(zero_one == min(all_zero, all_one, zero_one, one_zero))
                C(h-1, i) = 1;
                C(h-1, j) = 0;
                C(h, i) = 0;
                C(h, j) = 1;
                C(i, j) = 1;
                C(j, i) = 1;
             elseif(one_zero == min(all_zero, all_one, zero_one, one_zero))
                C(h-1, i) = 0;
                C(h-1, j) = 1;
                C(h, i) = 0;
                C(h, j) = 1;
                C(i, j) = 1;
                C(j, i) = 1;
             end
         end
      end
    end
    newG = GW;
    newC = C;
    newE = sum(sum(xor(GW, C)));
    if(newE < oldE)
        E = newE
    else
        E = 0; % E = 0 means that we have already reach the min
    end
end