function Subpred=sub_pls_pred(SubModels,Xt,yt)
% Rotina para previsão em sub-modelos PLS.
% input:
% SubModels : saida da rotina sub_pls_model.m
% exemplo: SubModels=sub_pls_model(X,y,SubGrupos,4,'center');
% X : matriz Xtest contendo os espectros das amostras de previsão;
% yt : vetor ytest contendo a propriedade de interesse das amostras de
% previsão;
%
% output:
% Subpred: estrutura contendo a previsão da propriedade de interesse de
% todas as amotras de teste (Xt).
%
%  
% exemplo:
% SubGrupos = sub_grupos(X,10,'euclidean')
% SubModels=sub_pls_model(X,y,SubGrupos,4,'center');
% Subpred=sub_pls_pred(SubModels,Xt,yt);
% 
% Paulo R. Filgueiras  - 28/01/2016
%

SubGrupos=SubModels.SubGrupos;
MTR=SubModels.SubGrupos.MTR;

%----- Determinando os espectros médios.
Xmedio=[];
for ki=1:size(SubGrupos.Samplesi,1)
    xm=mean(SubModels.X(SubGrupos.Samplesi{ki,1},:));
    Xmedio=[Xmedio;xm];
end
clear xm
%----- Calculando as distâncias e prevendo a propriedade de interesse
ypred=[];
for ki=1:size(Xt,1)
    SBP.d1 = pdist([Xt(ki,:);Xmedio],MTR);
    SBP.d2 = squareform(SBP.d1);
    SBP.dist=SBP.d2(2:end,1);
    [~,Grupo]=min(SBP.dist);

    ypred2=plsval(SubModels.models{Grupo,1},Xt(ki,:),yt(ki,1));
    ypred2=[Grupo,yt(ki,1),ypred2];
    ypred=[ypred;ypred2];
    clear ypred2
end

[fm.RMSEP,fm.R2p,fm.biasp]=rmse(ypred(:,2),ypred(:,3));


Subpred.pred=ypred;
Subpred.fm=fm;

