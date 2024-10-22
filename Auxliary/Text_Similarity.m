function A=Text_Similarity(X)

% % Important: evry column is normalized before any calculation

% X: inputed data set, each row is a sample
n = size(X,1);% number of samples
nn = 7; knn = floor(log2(n)) + 1; %default value

% normalized every row of X to be unit
d = sqrt(sum(X.^2, 2));
d=1./d;
X=bsxfun(@times, X, d);
% calculate distant matrix D
tempx = full(sum(X.^2, 2));
tempc = full(sum(X.^2, 2)');
D= tempx(:, ones(1,n)) + ...
    tempc(ones(1,n), :) - ...
    2.*(X*(X'));

Weigh=X*X';

% calculate the similarity matrix
[sorted, idx] = sort(D);% sort every column of D
% ls = sorted(nn+1, :); % find the first 7 nearest samples
% ls = sqrt(ls)';% calculate sigma_j
% get the first knn distant indexes
i = idx(1:knn+1, :);
j = kron(ones(knn+1,1),[1:n]);

i = i(:);
j = j(:);
I = find(i ~= j);% A(i,i)=0
i = i(I);
j = j(I);
index = [i, j; j, i];


% s = s(:);
% s = s(I);
% s = [s; s];
% % A_s = exp( -s ./ (ls(index(:,1)).*ls(index(:,2))) );
% A_s=s;


[index, i, j] = unique(index, 'rows');
for k=1:size(index,1)
    xind=index(k,1);
    yind=index(k,2);
    A_s(k)=X(xind,:)*X(yind,:)';
end



% A_s = A_s(i);
A = sparse(index(:,1), index(:,2), A_s, n, n);

% normalize the similarity matrix
dd = 1 ./ sum(A);
dd = sqrt(dd);
A = bsxfun(@times, A, dd);
A = A';
A = bsxfun(@times, A, dd);
A = (A + A') / 2;