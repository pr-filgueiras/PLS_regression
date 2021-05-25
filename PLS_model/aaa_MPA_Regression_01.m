% Aplicação da ideia do MPA para regressão
% Ideia:
% Construir um grande número de sub-modelos (regressores) e aplicar cada
% regressor ao conjunto de teste.
% 1. Utilizar o bootstrap para selecionar amostras.
% 2. Selecionar um subconjunto de variáveis.
% 3. Construir um regressor. 
% 4. Aplicar o sub-modelo de regressão aos conjunto de calibração e teste.
% 5. Repetir os passos 1-4 um número N de vezes.

% Uma variante pode ser utilizar esta ideia para identificação de outliers.
% Aplicar o regressor na matriz X nas amostras que não estão no modelo

cd('G:\Meu Drive\Analise Multivariada\Modulo 03 - Calibração Multivariada')
load('dados_prova.mat')

plot(num,X)
cd('G:\Meu Drive\Analise Multivariada\Banco de dados\libPLS')
[objetos,Xcal,Xtest,ycal,ytest]=caltest(X,y,70,'k',0,{'center'});

nvl = 10;
kfold = 5;
method = 'center';

CV = plscv(Xcal,ycal,nvl,kfold,method);
plot(CV.RMSECV),hold on , plot(CV.RMSECV,'ro')
xlabel('# latent variable')
ylabel('RMSECV')

nvl = CV.optLV_1SD;

modelo=pls(Xcal,ycal,nvl,method);  
[ypred,RMSEP]=plsval(modelo,Xtest,ytest); 
[fm.RMSEC,fm.R2c,fm.biasc] = rmse(ycal,modelo.y_fit,nvl);
[fm.RMSEP,fm.R2p,fm.biasp] = rmse(ytest,ypred)

%+++ Gráfico Referência vs Previsto
plot(ycal,modelo.y_fit,'bo'); hold on
plot(ytest,ypred,'r*')
refline(1,0)
legend('Calibration','Prediction')
xlabel('Reference')
ylabel('Predicted')


subplot(2,1,1), plot(X')
subplot(2,1,2), plot(modelo.VIP)
xlabel('Wavelength nm')

%+++++ Construção do MPA-PLS +++++
%+++++ Seleção das variáveis
nvar = 50; % Porcentagem do número de variáveis a serem utilizadas na modelagem
var_sel = round((nvar/100)*size(Xcal,2)); % Variáveis selecionadas
var_sel = randperm(size(Xcal,2),var_sel); var_sel = sort(var_sel);
%+++++ Seleção das amostras
nrep = 120; % Número de regressores utilizados.
[out]=bootrsp(1:size(Xcal,1),nrep);

%+++++ Construindo os regressores.
nvl = 10;
kfold = 5;
method = {'msc','center'};
ytest_out_bag = [];
ypred_out_bag = [];
for ki=1:nrep
    in_bag = unique(out(:,ki));
    outbag = setxor(1:size(Xcal,1),in_bag);
    %+++ Separação dos sub-conjuntos
    Xcal_in_bag = Xcal(out(:,ki),:);
    ycal_in_bag = ycal(out(:,ki),1);
    Xcal_out_bag = Xcal(outbag,:);
    ycal_out_bag = ycal(outbag,1);
    % Pre-processamento dos dados espectrais e seleção das variáveis
    % Pre-processa e depois seleciona as variáveis
    [Xcal_in_bag2,Xcal_out_bag] = pretrat(Xcal_in_bag,Xcal_out_bag,method);
    [Xcal_in_bag,Xtest_out] = pretrat(Xcal_in_bag,Xtest,method);
    Xcal_in_bag = Xcal_in_bag(:, var_sel);
    Xcal_out_bag = Xcal_out_bag(:, var_sel);
    Xtest_out = Xtest_out(:, var_sel);
    %+++ Determinação do nvl do regressor PLS
    CV = plscv(Xcal_in_bag,ycal_in_bag,nvl,kfold,'center');
    %nvl = CV.optLV_1SD;
    nvl = 7;
    %+++ Construção do modelo PLS
    modelo = pls(Xcal_in_bag,ycal_in_bag,nvl,'center');  
    ypred = plsval(modelo,Xcal_out_bag,ycal_out_bag);
    ytest_pred = plsval(modelo,Xtest_out,ytest);
    %+++ Armazenando os resultados 
    ypred_out_bag = [ypred_out_bag; [outbag,ycal_out_bag,ypred]];
    ytest_out_bag = [ytest_out_bag, ytest_pred];
end

% calculando
[~,a1] = sort(ypred_out_bag(:,1),1);
ypred_out_bag = ypred_out_bag(a1,:);

plot(ypred_out_bag(:,2),ypred_out_bag(:,3),'k.')
refline(1,0)

yhat_out_bag = [];
yhat_unique = unique(ypred_out_bag(:,1));
for ki=1:length(yhat_unique)
    a1 = find(ypred_out_bag(:,1)==ki);
    yhat_sub = [mean(ypred_out_bag(a1,[2,3])) std(ypred_out_bag(a1,3))];
    yhat_out_bag = [yhat_out_bag;[ki,length(a1),yhat_sub,min(ypred_out_bag(a1,3)),max(ypred_out_bag(a1,3))]];
end

ytest_out_bag_hat = [ytest, mean(ytest_out_bag,2) std(ytest_out_bag,0,2)];

[fm2.RMSEC,fm2.R2c,fm2.biasc] = rmse(ypred_out_bag(:,2),ypred_out_bag(:,3));
[fm2.RMSEP,fm2.R2p,fm2.biasp] = rmse(ytest_out_bag_hat(:,1),ytest_out_bag_hat(:,2))

% yhat_out_bag
% (1): número da amostra;
% (2): número de vezes utilizadas em out bag;
% (3): valor y de referência;
% (4): valor y médio calculado pelo out bag;
% (5): valor de desvio padrão de y calculado pelo out bag;
% (6): valor mínimo de y calculado pelo out bag;
% (7): valor máximo de y calculado pelo out bag;

plot(ytest_out_bag_hat(:,1),ytest_out_bag_hat(:,2),'ko')
refline(1,0)
