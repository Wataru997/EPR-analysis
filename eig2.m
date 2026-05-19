function [V,Ds] = eig2(A)
[V,D] = eig(A);
[~,I] = sort(diag(real(D)),'descend');
V = V(:, I);
Ds = D(I, I);