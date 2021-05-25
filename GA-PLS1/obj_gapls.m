function objv=obj_gapls(variable,Xcal,ycal,nvl,method,kfold)
% objv=obj_gapls(variable,Xcal,ycal,nvl,method,kfold)
% obj_gapls: fun��o objetiva para a modelagem GA-PLS.
%+++ input
% variable: vari�veis selecionadas para c�lculo do fitness.
% Xcal: Matriz contendo os dados espectrais.
% ycal: vetor contendo a vari�vel a ser modelada.
% nvl: n�mero de vari�veis latentes a ser avaliada na modelagem.
% method: m�todo de pr�-processamento de vari�veis a ser aplicado a matriz Xcal.
% kfold: n�mero de partes "k" do m�todo k-fold.
%+++ output
% objv: valores de aptid�o de cada cromossomo, utilizando o par�metro
% RMSECV.
%
% Exemplo: objv=obj_gapls(variable,Xcal,ycal,7,{'center'},5)
%
% Paulo Roberto Filgueiras 30/12/2020

% Pre-processamento dos dados Xcal
Xcal2 = pretrat(Xcal,[],method);

% Determina��o dos valores de aptid�o para cada cromossomo.
% C�lculo da aptid�o utilizando o par�metro RMSECV
objv = [];
for ki = 1:size(variable,1)
    variable2 = find(variable(ki,:));
    rmsecv_i = plscv(Xcal2(:,variable2),ycal,0,nvl,kfold);
    objv = [objv; min(rmsecv_i)];
end
