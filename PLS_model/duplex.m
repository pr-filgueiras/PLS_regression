
% -------------------------------------------------------------------------
% Function: [model,test]=duplex(X,k)
% -------------------------------------------------------------------------
% Aim:
% Subset selection with Duplex algorithm; uniform design of model and test
% sets
% -------------------------------------------------------------------------
% Input:
% X, matrix (n,p), predictor variables in columns
% k, number of objects to be selected to test set, (test set can contain 
% at most 0.5n objects. If less objest are selected to test set, than k 
% first objects of the model and test sets are designed uniformly and the
% remaining objects not selected by Duplex algorithm are included 
% into model set)
% -------------------------------------------------------------------------
% Output:
% model, vector (k+(n-k),1), list of objects selected to model set
% test, vector (k,1), list of objects selected to test set (optionally)
% -----------------------------------------------------------------------
% Example: 
% [model,test]=duplex(X,10)
% -----------------------------------------------------------------------
% Reference:
% R.D. Snee, Technometrics 19 (1977) 415-428

% Written by Michal Daszykowski
% Department of Chemometrics, Institute of Chemistry, 
% The University of Silesia
% December 2004

function [model,test]=duplex1(X,k)

[m,n]=size(X);
ma=floor(0.5*m);

if k>ma
    k=ma;
end

x=[[1:size(X,1)]' X];
n=size(x,2);

% Fist two most distant points to model set
p=tril(fastdist(x(:,2:n),x(:,2:n)));
[i1 i2]=find(p==max(max(p)));
model=x([i1 i2],1);
x([i1 i2],:)=[];

% Another two most distant points to test set
p=tril(fastdist(x(:,2:n),x(:,2:n)));
[i1 i2]=find(p==max(max(p)));
test=x([i1 i2],1);
x([i1 i2],:)=[];


h=waitbar(0,'Please wait ...'); 
h=waitbar(0/k,h);
iter=2;

while length(model)<k
    [ii,ww]=max(min(fastdist(x(:,2:n),X(model,:))));
    model=[model;x(ww,1)];
    x(ww,:)=[]; 
    [ii,ww]=max(min(fastdist(x(:,2:n),X(test,:))));
    test=[test;x(ww,1)];
    x(ww,:)=[];
    iter=iter+1;
    h=waitbar(iter/k,h);
end

if ~isempty(x);
    model=[model;x(:,1)];
end
    
close(h);


function D=fastdist(x,y)

% Calculated Euclideam distances between two sets of objetcs

D=((sum(y'.^2))'*ones(1,size(x,1)))+(ones(size(y,1),1)*(sum(x'.^2)))-2*(y*x');