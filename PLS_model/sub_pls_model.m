function SubModels=sub_pls_model(X,y,SubGrupos,nvl,method)
% Rotina para construir sub-modelos PLS
% input:
% X : matriz Xcal contendo os espectros das amostras de calibração;
% y : vetor ycal contendo a propriedade de interesse das amostras de
% calibração;
% SubGrupos: saida da rotina sub_grupos.m;
% nvl : número de variáveis latentes;
% method: método de pre-processamento dos dados (default: 'center').
%
% output:
% SubModels: estrutura com todos os submodelos.
%
% exemplo:
% SubGrupos = sub_grupos(X,10,'euclidean')
% SubModels=sub_pls_model(X,y,SubGrupos,4,'center');
% Subpred=sub_pls_pred(SubModels,Xt,yt);
%
% Paulo R. Filgueiras  - 28/01/2016
%  
if nargin==5;   method='center'; end

for ki=1:size(SubGrupos.Samplesi,1)
    Nsamples=SubGrupos.Samplesi{ki,1};
    Xcal=X(Nsamples,:);
    ycal=y(Nsamples,1);

    PLS=pls(Xcal,ycal,nvl,method);  
    [fm.yc,RMSEC]=plsval(PLS,Xcal,ycal); 
    [fm.RMSEC,fm.R2c,fm.biasc]=rmse(ycal,PLS.y_fit,nvl);
    fm.yc=[ycal,fm.yc];
    SubModels.models{ki,1}=PLS;
    SubModels.models{ki,2}=fm;
end
SubModels.X=X;
SubModels.y=y;
SubModels.SubGrupos=SubGrupos;
