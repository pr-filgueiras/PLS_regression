function [objetos,Xcal,Xtest,ycal,ytest]=caltestda(X,y,ncal,alg,rep,method)
% algoritmo para selecionar amostras de calibração e teste para problemas
% de classificação. 
%      X : matriz de dados espectrais; 
%      y : vetor contendo as classes; (ordenado)
%   ncal : percentagem de amostras de treinamento; default: 70
%    alg : tipo de algoritmo a ser utilizado para repartição dos dados:
%          kenston('k'); duplex('d') ou segmentado (s); default: 'k'
%    rep : vetor indicando as réplicas espectrais; default: [] sem réplicas;
% method : método de pre-processamento dos dados antes da separação dos 
%          conjuntos de calibração e teste; default: 'none'.
% 
% exemplos:
% [objetos,Xcal,Xtest,ycal,ytest]=caltestda(X,y,70,'k',rep,{'none'});
% [objetos,Xcal,Xtest,ycal,ytest]=caltestda(X,y,70,'k',rep);
% [objetos,Xcal,Xtest,ycal,ytest]=caltestda(X,y,70,'d',rep,);
% [objetos,Xcal,Xtest,ycal,ytest]=caltestda(X,y,70,'s',rep,{'none'});
% [objetos,Xcal,Xtest,ycal,ytest]=caltestda(X,y,70,'k');
%
% Paulo R. Filgueiras  - 13/08/2014
% 

if nargin<3, ncal=70; end
if nargin<4; alg='k'; end
if nargin<5; rep=[]; end
if nargin<6;  method={'none'}; end


% pre-processamento dos espectros.
X2=pretrat(X,[],method);

% verificação das réplicas a partir do vetor rep
if isvector(rep)
    y2.sample=[];y2.sample=[y2.sample;1];
    y2.numero=[];
    for ki=2:size(X2,1)
        if sum(rep(ki,:)==rep(ki-1,:))==size(rep,2)
            y2.numero=ki;
        else
            y2.sample=[y2.sample;ki];
        end
    end
    aa2=y2.sample(2:end)-1;aa2=[aa2;y2.numero];
    y2.sample=[y2.sample,aa2,aa2-y2.sample+1];  % vetor com as posições iniciais e finais para as réplicas.
    clear aa2 ki
    % Calculando a média espectral 
    y2.medioX=[];y2.medioy=[];y2.class=[];
    for ki=1:size(y2.sample,1)
        posicao_medio = y2.sample(ki,1):y2.sample(ki,2);
        matriz_medio = mean(X2(posicao_medio,:));
        y2.medioX=[y2.medioX;matriz_medio];         % matriz com os espectros médios.
        y2.medioy=[y2.medioy;rep(y2.sample(ki,1))]; 
        y2.class=[y2.class;y(y2.sample(ki),:)];     % vetor com as classes.
    end
else
    y2.medioX=X2;    % matriz com os espectros médios.
    y2.medioy=ones(size(X2,1),1);   
    y2.class=y;      % vetor com as classes.
end

% determinando o número de amostras em cada classe a ser selecionada para
% treinamento
classes.classes=unique(y,'rows');
classes.amostras{1,1}='Classes';
classes.amostras{1,2}='Amostras';
classes.amostras{1,3}='Treinamento';
for ki=1:length(classes.classes)
    amostras = find(y2.class==classes.classes(ki));
    classes.amostras{ki+1,1}=classes.classes(ki);
    classes.amostras{ki+1,2}=amostras;
    classes.amostras{ki+1,3}=round(ncal*length(amostras)/100);
end

% Selecionando as amostras de treinamento pelo algoritmo
%  alg : tipo de algoritmo a ser utilizado para separação dos dados:
classes.train=[];
for ki=1:length(classes.classes)
    X3=y2.medioX(classes.amostras{ki+1,2},:);
    y3=y2.medioy(classes.amostras{ki+1,2},:);
    if strcmp(alg,'k')
        classes.train{ki,1}=classes.amostras{ki+1,1};
        classes.train{ki,2}=kenston(X3,classes.amostras{ki+1,3},1,0,y3); %não gera gráfico de saída
        classes.train{ki,3}=setxor(1:length(y3),classes.train{ki,2});
        classes.train{ki,2}=classes.amostras{ki+1,2}(classes.train{ki,2});
        classes.train{ki,2}=sort(classes.train{ki,2}); %classes.train{ki,2}=classes.train{ki,2}';  
        classes.train{ki,3}=classes.amostras{ki+1,2}(classes.train{ki,3});
        classes.train{ki,3}=sort(classes.train{ki,3});
    elseif strcmp(alg,'d')
        classes.train{ki,1}=classes.amostras{ki+1,1};
        classes.train{ki,2}=duplex(X3,length(y3)-classes.amostras{ki+1,3});
        classes.train{ki,3}=setxor(1:length(y3),classes.train{ki,2});
        classes.train{ki,2}=classes.amostras{ki+1,2}(classes.train{ki,2});
        classes.train{ki,2}=sort(classes.train{ki,2}); %classes.train{ki,2}=classes.train{ki,2}';  
        classes.train{ki,3}=classes.amostras{ki+1,2}(classes.train{ki,3});
        classes.train{ki,3}=sort(classes.train{ki,3});
    elseif strcmp(alg,'s')
        classes.train{ki,1}=classes.amostras{ki+1,1};
        aa2=repmat((1:ncal)',length(y3),1);  aa2=aa2(1:length(y3));
        classes.train{ki,2}=find(aa2~=1);
        classes.train{ki,3}=find(aa2==1);
        classes.train{ki,2}=classes.amostras{ki+1,2}(classes.train{ki,2});
        classes.train{ki,2}=sort(classes.train{ki,2}); %classes.train{ki,2}=classes.train{ki,2}';  
        classes.train{ki,3}=classes.amostras{ki+1,2}(classes.train{ki,3});
        classes.train{ki,3}=sort(classes.train{ki,3});
    end
end

% Organizando as amostras de treinamento e teste.
classes.train2=[];
classes.teste2=[];
for ki=1:length(classes.classes)
    classes.train2=[classes.train2;classes.train{ki,2}];
    classes.teste2=[classes.teste2;classes.train{ki,3}];
end
classes.train2=sort(classes.train2);   
classes.teste2=sort(classes.teste2);

% Separando as matrizes de treinamento e teste.
if isvector(rep)
    objetos.treinamento=[];
    for ki=1:length(classes.train2);
        aa1=classes.train2(ki);
        aa2=y2.sample(aa1,1):y2.sample(aa1,2);
        objetos.treinamento=[objetos.treinamento,aa2];
    end
    objetos.treinamento=objetos.treinamento';
    objetos.teste=setxor(1:size(X,1),objetos.treinamento);
else
    objetos.treinamento=classes.train2;
    objetos.teste=classes.teste2;
end

objetos.classes=classes;
objetos.classes2=y2;

Xcal=X(objetos.treinamento,:);
Xtest=X(objetos.teste,:);
ycal=y(objetos.treinamento,:);
ytest=y(objetos.teste,:);

function [object,xm,ym,xt,yt]=kenston(x,no_p,men,pl,y)
%#									
%#  function [object,xm,ym,xt,yt]=kenston(x,no_p,men,pl,y)				
%#									
%#  AIM: 	Kennard-Stone design (e.g. for subset selection).	
%#									
%#  PRINCIPLE:  Based on Computer Aided Design of Experiments. 		
%#		The first point can be the mean or the furthest from 	
%#		the mean.						
%# 		REF :   R. W. Kennard and L. A. Stone			
%# 			Technometrics Vol. 11, No. 1, 1969 		
%# 									
%#  INPUT:	x :  (n x m) absorbent matrix with n spectra		
%#			 and m variables				
%#		no_p : number of objects to be selected			
%#		men: position of the first point selected		
%#		     (1 closest to mean; 0 furthest from mean)		
%#		pl: 1 plot; 0 no plot					
%#		y: matrix of responses							
%#  OUTPUT:	object : (1 x no_p) the vector of designed objects	
%#																	
%#  AUTHOR: 	Wen Wu 							
%#	    	Copyright(c) 1997 for ChemoAC				
%#          	FABI, Vrije Universiteit Brussel            		
%#          	Laarbeeklaan 103 1090 Jette				
%#    	    								
%# VERSION: 1.1 (28/02/1998)						
%#									
%#  TEST:   	Roy De Maesschalck					
% Exemplo
%#	[object,Xcal,ycal,Xtest,ytest]=kenston(xsnv,30,1,1,y);								

[n,m]=size(x);	
t=x;

% Kennard and Stone method to select objects
  % starting (1st) point is the closest (or furthest) from the centroid
	meant=mean(t);
	t1=t-ones(n,1)*meant;
        for i=1:n
	    a(i)=t1(i,:)*t1(i,:)'; 
	end %i
	if men==1,
		[b,c]=min(a);
	else
		[b,c]=max(a);
	end
	object(1)=c;
	clear a b c t1

   % 2nd points is the furthest point from 1st point
	t1=t-ones(n,1)*t(object(1),:);
	for i=1:n
	    a(i)=t1(i,:)*t1(i,:)';
	end %i
	[b,c]=max(a);
	object(2)=c;
	clear a b c t1
	
  % k+1 point
	for pi=3:no_p
	    list=1:n;
	    tt=t;
	    k=length(object);
	    list(object)=[];
	    tt(object,:)=[];
	    nl=length(list);
	    for j=1:nl
		for i=1:k
	 	    t1=tt(j,:)-t(object(i),:);
		    a(i)=t1*t1';
		end % i
	    	[b,c]=min(a);
		dmin(j)=b;
	    end %j
	    [b,c]=max(dmin');
	    object(pi)=list(c);
	    clear dmin a b c list
	end %pi
object=object(1:no_p);

% plot
if pl==1
	[a,b]=size(t);
   if b>1
	plot(t(:,1),t(:,2),'.')
	for i=1:n
	    text(t(i,1),t(i,2),int2str(i))
	end
	hold on
	plot(t(object,1),t(object,2),'r*')
	hold off
	xlabel('Variable 1')
	ylabel('Variable 2')
   end 
   if b==1
	plot(t(:,1),t(:,1),'.')
	for i=1:n
	    text(t(i,1),t(i,1),int2str(i))
	end
	hold on
	plot(t(object,1),t(object,1),'r*')
	hold off
	xlabel('Variable 1')
	ylabel('Variable 1')
   end 
end                    
xm=x(object,:);
ym=y(object,:);
ind=[1:n]';
ind(object)=[];
xt=x(ind,:);
yt=y(ind,:);

function [model,test]=duplex(X,k)
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
