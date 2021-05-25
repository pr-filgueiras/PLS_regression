% Trabalho de quantificação do biodiesel de gordura animal em blend com
% biodiesel de soja. O teor máximo de ordura animal varia de 0% a 69%

% Abrir os espectros
[X,y,num,amostras,FileName]=leitor_espec('ASC','1.ASC','biodiesel.xlsx','bio',2);

% G:\Unicamp\Tool MatLab\PLS_Chemometrics
% separar os conjuntos de calibração e teste 
[objetos,Xcal,Xtest,ycal,ytest]=caltest(X,y,59,'k',0,{'none'});

% opções para otimização do número de variáveis latentes do modelo PLS 
options=[]; 
options.Xpretreat      = {'msc'}; % pretratamento da matriz Xcal; ycal is always mean center;
options.vene           = 5;       % Cross validation: k from venetian blind k-fold;
options.alpha          = 0.95;    % scope level of confidence level.
options.xaxis          = num;     % variáveis do eixo-X;
options.vl             = 10;      % número de variáveis latentes.
options.precisao       = [];      % matriz para o cálculo de repetibilidade.
options.graficos       = 1;       % output graphic.end 
    
modelo = plsmodel(Xcal,ycal,options);

options.vl=3;     % número ótimo de variáveis latentes.
modelo = plsmodel(Xcal,ycal,Xtest,ytest,options)

% Excluir outlier (Amostra 91) 
am=setxor(1:100,91);
X=X(am,:);
y=y(am,:);
amostras=amostras(am,:);
% retornar a linha 9


% teste de permutação para avaliação de tendência
perm = perm_teste(ytest,modelo.yp,10000,2)

%% Como construir um modelo PLS?
% Filgueiras PR, et al. Quantification of animal fat biodiesel in soybean biodiesel and B20 diesel 
% blends using near infrared spectroscopy and synergy interval support vector regression. 
% Talanta (Oxford), v. 119, p. 582-589, 2014.
%https://doi.org/10.1016/j.talanta.2013.11.056

% http://libpls.net/
cd('G:\Meu Drive\Analise Multivariada\Modulo 03 - Calibração Multivariada\Dados\exemplo_01')

% Leitura dos dados de teor de biodiesel de sebo bovino
ysebo = xlsread('biodiesel.xlsx','bio','B2:B101');
% Leitura dos dados espectrais
NIR = [];
files = dir('*.ASC');
for ki=1:numel(files)
    FileName{ki,1}=char(files(ki,1).name);
    espectro_i=dlmread(FileName{ki,1},'\t',56,0);
    NIR = [NIR; espectro_i(:,2)'];
end
num = espectro_i(:,1)';

plot(num,NIR)

% +++++ Construção do modelo +++++
% Separar as amostras em conjunto de Calibração e Teste.
% Utilização do algoritmo Kennard-Stone
cd('G:\Meu Drive\Analise Multivariada\Banco de dados\libPLS')
[info,Xcal,Xtest,ycal,ytest]=caltest(NIR,ysebo,70,'k',0,{'none'});
plot(ycal,ycal,'b*',ytest,ytest,'ro')
legend ('ycal','ytest'), grid on
xlabel ('Referência'); ylabel ('Referência')
% Utilização do algoritmo Rank_KS
cd('G:\Meu Drive\Analise Multivariada\Banco de dados\libPLS')
[Xcal,Xtest,ycal,ytest,info]=rank_ks(NIR,ysebo,FileName,70,10,1);
% [Xcal,Xtest,ycal,ytest,info]=rank_ks(X,y,FileName,ncal,inter,prop,graf)
% Pre-processamento dos dados.
cd('G:\Meu Drive\Analise Multivariada\Banco de dados\libPLS')
[Xcalp,Xtestp]=pretrat(Xcal,Xtest,{'msc'});
subplot(2,1,1), plot(num,Xcal,'b',num,Xtest,'r')
xlabel ('Número de onda (cm^{-1})'); ylabel ('Absorbância')
subplot(2,1,2), plot(num,Xcalp,'b',num,Xtestp,'r')
xlabel ('Número de onda (cm^{-1})'); ylabel ('Absorbância')

% Determinação do número correto de variáveis latentes - PLS
cd('G:\Meu Drive\Analise Multivariada\Banco de dados\libPLS')
nvl = 10;
kfold = 5;
method = 'center';
% Método k-fold cross validation
CV = plscv(Xcalp,ycal,nvl,kfold,method);
plot(CV.RMSECV),hold on , plot(CV.RMSECV,'ro')
grid on
xlabel('# variáveis latentes')
ylabel('RMSECV')
% Método Monte Carlo cross validation
Nmc = 1000;
ratiomc = 0.8; % 80% das amostras do conjunto de calibração.
MCCV = plsmccv(Xcalp,ycal,nvl,method,Nmc,ratiomc,1);
plot(MCCV.RMSECV),hold on , plot(MCCV.RMSECV,'ro')
grid on
xlabel('# variáveis latentes')
ylabel('RMSECV')

% Construção do modelo PLS
nvl = 4;
modelo=pls(Xcalp,ycal,nvl,method);  
[ypred,RMSEP]=plsval(modelo,Xtestp,ytest); 
[fm.RMSEC,fm.R2c,fm.biasc] = rmse(ycal,modelo.y_fit,nvl);
[fm.RMSEP,fm.R2p,fm.biasp] = rmse(ytest,ypred)
%Gráfico Referência vs Previsto
plot(ycal,modelo.y_fit,'bo'); hold on
plot(ytest,ypred,'r*')
refline(1,0), grid on
legend('Calibração','Previsão')
xlabel('Referência')
ylabel('Previsto')
%Gráfico dos resíduos
residuos_cal = ycal - modelo.y_fit; % Resíduos de calibração
residuos_test = ytest - ypred; % Resíduos de previsão
plot(ycal,residuos_cal,'bo'); hold on
plot(ytest,residuos_test,'r*')
hline(0,'k'), grid on
legend('Calibração','Previsão')
xlabel('Referência')
ylabel('Resíduos')

% Avaliação do modelo
% teste para viés
teste = bias_teste(ycal,modelo.y_fit,0.05)
% teste de permutação para tendência
perm = perm_teste(ycal,modelo.y_fit,10000,2)
% Determinação das figuras de mérito
saida = fmerito(modelo,ycal,ypred,fm.RMSEP)
% Intervalo de confiança
confidence = confidence_interval(Xcalp,Xtestp,ycal,{'msc'},nvl,0.05)
errorbar(ycal,confidence.cal(:,1),confidence.cal(:,2),'bo'), hold on
errorbar(ytest,confidence.test(:,1),confidence.test(:,2),'r*')
refline(1,0), grid on
legend('Calibração','Previsão')
xlabel('Referência')
ylabel('Previsto')
% Detecção de outliers
figure(1)
plot(1:length(ycal),confidence.leverage_cal,'b.'), hold on
plot(length(ycal)+1:length(ycal)+length(ytest),confidence.leverage_test,'r.'), hold on
for ki=1:length(ycal)
    text(ki,confidence.leverage_cal(ki),num2str(ki),'Color','b')
end
for ki=length(ycal)+1:length(ycal)+length(ytest)
    m=ki-length(ycal);
    text(ki,confidence.leverage_test(m),num2str(m),'Color','r')
end
hline(confidence.leverage_limit,'k'), grid on
legend('Calibração','Previsão')
xlabel('Amostra')
ylabel('Leverage')

figure(2) % Leverage vs residuo
plot(confidence.leverage_cal,abs(residuos_cal),'b.'), hold on
plot(confidence.leverage_test,abs(residuos_test),'r.')
for ki=1:length(ycal)
    text(confidence.leverage_cal(ki),abs(residuos_cal(ki)),num2str(ki),'Color','b')
end
for ki=length(ycal)+1:length(ycal)+length(ytest)
    m=ki-length(ycal);
    text(confidence.leverage_test(m),abs(residuos_test(m)),num2str(m),'Color','r')
end
vline(confidence.leverage_limit,'k'), grid on
legend('Calibração','Previsão')
xlabel('Leverage')
ylabel('Residuo')

% Aplicação do modelo em novas amostras 
Xtest2 = Xtest(3:5,:); % Exemplo de amostras de teste.
[~,Xtest2p]=pretrat(Xcal,Xtest2,{'msc'}); % Pre-processamento dos dados
ypred2 = plsval(modelo, Xtest2p, ones(size(Xtest2p,1),1)  )
confidence = confidence_interval(Xcalp,Xtest2,ycal,{'msc'},nvl,0.05);
confidence.test
confidence.leverage_test
confidence.leverage_limit


%% Utilização da rotina regression_toolbox_1.3
% http://www.michem.unimib.it/
cd('G:\Meu Drive\Analise Multivariada\Banco de dados\regression_toolbox_1.3')
reg_gui

