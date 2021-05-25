load('..\libs\');

%% Como construir um modelo GA-PLS?
% Leitura dos dados de teor de biodiesel de sebo bovino
ysebo = xlsread('biodiesel.xlsx','bio','B2:B101');
% Leitura dos dados espectrais
NIR = [];
files = dir('ascFiles/*.ASC');
for ki=1:numel(files)
    FileName{ki,1}=char(files(ki,1).name);
    espectro_i=dlmread(FileName{ki,1},'\t',56,0);
    NIR = [NIR; espectro_i(:,2)'];
end
num = espectro_i(:,1)';

plot(num,NIR);

% +++++ Constru��o do modelo +++++
% Separar as amostras em conjunto de Calibra��o e Teste.
% Utiliza��o do algoritmo Kennard-Stone


[info,Xcal,Xtest,ycal,ytest]=caltest(NIR,ysebo,70,'k',0,{'none'});
plot(ycal,ycal,'b*',ytest,ytest,'ro')
legend ('ycal','ytest'), grid on
xlabel ('Refer�ncia'); ylabel ('Refer�ncia')
% Utiliza��o do algoritmo Rank_KS
[Xcal,Xtest,ycal,ytest,info]=rank_ks(NIR,ysebo,FileName,70,10,1);

% Pre-processamento dos dados.
[Xcalp,Xtestp]=pretrat(Xcal,Xtest,{'msc'});
subplot(2,1,1), plot(num,Xcal,'b',num,Xtest,'r')
xlabel ('N�mero de onda (cm^{-1})'); ylabel ('Absorb�ncia')
subplot(2,1,2), plot(num,Xcalp,'b',num,Xtestp,'r')
xlabel ('N�mero de onda (cm^{-1})'); ylabel ('Absorb�ncia')

% Construindo o modelo GA-PLS

options.nind = 30; % popula��o inicial
options.Nvar = size(Xcal,2);
options.var_init = 10 ; % 10 porcento da popula��o inicial
options.nvl = 8;
options.method = {'center'};
options.kfold = 5;
options.maxgen = 25;
options.ggap = 0.8;
options.Nga = 30;

modelo = ga_pls(Xcalp,ycal,options)
best_var = sum(modelo.best_var);
mm = [min(best_var) max(best_var)]
RMSECV = [];
ki = 1;
while ki <= max(best_var)
    variable = find(best_var>=ki);
    if length(variable) > options.nvl
        objv=obj_gapls(variable,Xcalp,ycal,options.nvl,options.method,options.kfold);
        RMSECV = [RMSECV; [ki length(variable) objv] ];
        ki = ki + 1;
    else 
        return
    end
end

plot(RMSECV(:,1),RMSECV(:,3),'r',RMSECV(:,1),RMSECV(:,3),'bo')
variable = find(best_var>=8);
objv=obj_gapls(variable,Xcalp,ycal,options.nvl,options.method,options.kfold);
% Agora basta construir o modelo PLS.
% Primeiro pre-processa os dados e depois seleciona as vari�veis.


nvl = 10;
kfold = 5;
method = 'center';
% M�todo Monte Carlo cross validation
Nmc = 1000;
ratiomc = 0.8; % 80% das amostras do conjunto de calibra��o.
MCCV = plsmccv(Xcalp(:,variable),ycal,nvl,method,Nmc,ratiomc,1);
plot(MCCV.RMSECV),hold on , plot(MCCV.RMSECV,'ro')
grid on
xlabel('# vari�veis latentes')
ylabel('RMSECV')

nvl = 4;
modelo=pls(Xcalp(:,variable),ycal,nvl,method);  
[ypred,RMSEP]=plsval(modelo,Xtestp(:,variable),ytest); 
[fm.RMSEC,fm.R2c,fm.biasc] = rmse(ycal,modelo.y_fit,nvl);
[fm.RMSEP,fm.R2p,fm.biasp] = rmse(ytest,ypred)
%Gr�fico Refer�ncia vs Previsto
plot(ycal,modelo.y_fit,'bo'); hold on
plot(ytest,ypred,'r*')
refline(1,0), grid on
legend('Calibra��o','Previs�o')
xlabel('Refer�ncia')
ylabel('Previsto')
%Gr�fico dos res�duos
residuos_cal = ycal - modelo.y_fit; % Res�duos de calibra��o
residuos_test = ytest - ypred; % Res�duos de previs�o
plot(ycal,residuos_cal,'bo'); hold on
plot(ytest,residuos_test,'r*')
hline(0,'k'), grid on
legend('Calibra��o','Previs�o')
xlabel('Refer�ncia')
ylabel('Res�duos')

% Avalia��o do modelo
% teste para vi�s
teste = bias_teste(ycal,modelo.y_fit,0.05)
% teste de permuta��o para tend�ncia
perm = perm_teste(ycal,modelo.y_fit,10000,2)
% Determina��o das figuras de m�rito
saida = fmerito(modelo,ycal,ypred,fm.RMSEP)
% Intervalo de confian�a
confidence = confidence_interval(Xcalp(:,variable),Xtestp(:,variable),ycal,{'msc'},nvl,0.05)
errorbar(ycal,confidence.cal(:,1),confidence.cal(:,2),'bo'), hold on
errorbar(ytest,confidence.test(:,1),confidence.test(:,2),'r*')
refline(1,0), grid on
legend('Calibra��o','Previs�o')
xlabel('Refer�ncia')
ylabel('Previsto')
% Detec��o de outliers
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
legend('Calibra��o','Previs�o')
xlabel('Amostra')
ylabel('Leverage')



%% Construindo um modelo GA-PLS com outra rotina
options.nopop    = 30;    %Population size
options.maxgen   = 25;    %Max number of generations
options.mut      = 0.01;  %Mutation Rate
options.window   = 5;     %Window width for spectral channels
options.converge = 80;    % % of pop the same at convergence
options.begfrac  = 50;    %Fraction of terms included in beginning
options.cross    = 2;     %Double or single cross over, 1 = single
options.maxlv    = 10;    %No. LVs, only needed with reg = 1
options.cvopt    = 0;     %CV option, 0 = random, 1 = contiguous
options.split    = 5;     %No. subsets to divide data into
options.iter     = 1;     %No. iterations of CV at each generation
options.Nga      = 10;     % n�mero de vezes que o GA ser� rodado.

gamodel = ga_pls(Xcalp,ycal,options)

