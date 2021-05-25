function objv=obj_gapls(variable,Xcal,ycal,nvl,method,kfold)
% objv=obj_gapls(variable,Xcal,ycal,nvl,method,kfold)
% obj_gapls: função objetiva para a modelagem GA-PLS.
%+++ input
% variable: variáveis selecionadas para cálculo do fitness.
% Xcal: Matriz contendo os dados espectrais.
% ycal: vetor contendo a variável a ser modelada.
% nvl: número de variáveis latentes a ser avaliada na modelagem.
% method: método de pré-processamento de variáveis a ser aplicado a matriz Xcal.
% kfold: número de partes "k" do método k-fold.
%+++ output
% objv: valores de aptidão de cada cromossomo, utilizando o parâmetro
% RMSECV.
%
% Exemplo: objv=obj_gapls(variable,Xcal,ycal,7,{'center'},5)
%
% Paulo Roberto Filgueiras 30/12/2020

% Pre-processamento dos dados Xcal
Xcal2 = pretrat(Xcal,[],method);

% Determinação dos valores de aptidão para cada cromossomo.
% Cálculo da aptidão utilizando o parâmetro RMSECV
objv = [];
for ki = 1:size(variable,1)
    variable2 = find(variable(ki,:));
    rmsecv_i = plscv(Xcal2(:,variable2),ycal,0,nvl,kfold);
    objv = [objv; min(rmsecv_i)];
end
