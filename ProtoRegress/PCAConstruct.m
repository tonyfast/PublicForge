function [PC,Var,Base]=PCAConstruct(Data,Rank)
% A wrapper around SVD for PCA. ROWS as data points, COLUMNS as dimensions. 
K=Data;
M=repmat(mean(K,1),size(K,1),1);
Data=K-M;
[U S V] = pca(Data, Rank);
PC=U*S;
Var=S^2;
Base=S*V;