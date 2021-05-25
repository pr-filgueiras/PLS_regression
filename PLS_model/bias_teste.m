function teste = bias_teste(valor_real,valor_previsto,alpha)
% Teste para erros sistem�ticos - bias
% input:
%    valor_real: valor de refer�ncia
%    valor_previsto: resultado calculado (modelado)
%    alpha: valor de alpha teste bivaudal segundo a distribui��o t-studente (padr�o 0.05)
% 
% teste = bias_teste(valor_real,valor_previsto,alpha);
%
% Paulo R. Filgueiras 03/12/2012
%
if nargin==2
    alpha=0.05;
end
teste.bias=(sum(valor_real-valor_previsto))/length(valor_real);
teste.SVD=sqrt(sum((valor_real-valor_previsto-teste.bias).^2)/(length(valor_real)-1));
teste.t=(abs(teste.bias)*sqrt(length(valor_real)))/teste.SVD;
%% teste para erro sistem�tico
alpha=alpha/2;
teste.ttab=abs(tinv(alpha,(length(valor_real)-1)));
disp('  ')
if teste.t < teste.ttab
    s = sprintf('tcal = %g < ttab = %g',teste.t,teste.ttab); disp(s)
    disp('Erros sistem�ticos N�O significativos')
else
    s = sprintf('tcal = %g > ttab = %g',teste.t,teste.ttab); disp(s)
    disp('Os erros sistem�ticos s�o significativos')
end
disp('  ')
teste.pvalor = 2*(1-tcdf(teste.t,length(valor_real)-1));
