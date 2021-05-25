function [Sample,Fator] = outlier_detec2(X,y,nvl,graf)
% Rotina para determinar influência de cada amostra em modelos PLS.
% Quanto maior a influência da amostra maior sua probabilidade de ser um 
% outlier.
% Rotina utilizada junta com a libPLS
% input:
%  X : matriz com os dados espectrais.
%  y : vetor com a propriedade de interesse.
% nvl: número de variáveis latentes a ser testada.
%graf: saida gráfica (0 = não, 1 sim [default])
%
% output
% Samples : vetor contendo a ordem de influência. A primeira amostra 
% representa a mais influente.
%   Fator : influência de cada amostra.
%
% Exemplo
% [Samples,Fator] = outlier_detec(X,y,nvl,1);
% [Samples,Fator] = outlier_detec(X,y,4,0);
%
% Paulo R. Filgueiras  - 19/01/2016

if nargin==3; graf=1; end


Amostras.samples=1:length(y);



X2 = X(Amostras.samples,:);
y2 = y(Amostras.samples);


cal2=setxor(1:length(Sample),Sample1);
X2 = X(cal2,:);
y2 = y(cal2);

Fator=[]; ypred2=[];
for ki=1:size(X2,1)
    Xcal=X2(setxor(1:size(X2,1),ki),:);
    ycal=y2(setxor(1:size(X2,1),ki),:);
    Xtest=X2(ki,:);
    ytest=y2(ki,:);

    PLS=pls(Xcal,ycal,nvl);
    ypred=plsval(PLS,Xtest,ytest);
    RMSEC=rmse(ycal,PLS.y_fit,nvl);

    Fator=[Fator;RMSEC];
    ypred2=[ypred2;ypred];
end
RMSECV_loo=rmse(y2,ypred2)
[~,Sample]=sort(Fator,'descend'); % Amostras mais influentes

Sample1 = Amostras.samples(Sample(1));





if graf==1 
figure(1),plot(Fator,'k.'), hline(1,'r:')
for ki=1:size(X,1)
    text(ki,Fator(ki),num2str(ki));
end
xlabel('Sample','FontSize',12)
ylabel('RMSEC(full)/RMSEC','FontSize',12)

figure(2),plot(ypred2,Fator,'k.')
for ki=1:size(X,1)
    text(ypred2(ki),Fator(ki),num2str(ki));
end
xlabel('Error^2','FontSize',12)
ylabel('RMSEC(full)/RMSEC','FontSize',12)

figure(3),plot(y,ypred2,'k.'), vline(mean(y),'r:')
for ki=1:size(X,1)
    text(y(ki),ypred2(ki),num2str(ki));
end
ylabel('Error^2','FontSize',12)
xlabel('Reference','FontSize',12)
end

































