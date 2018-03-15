function C = allcomb_cell(A,B)
% A is a m by 1 cell array
% B is a n by 1 cell array
% C is mn by 2 cell array

m = size(A,1);
n = size(B,1);

C=cell(m*n,2);
for i=1:m
    for j = 1:n
        C((i-1)*n+j,:) = {A{i,1} B{j,1}};
    end
end
