%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  This script is used to test whether the functions in this  %%%%%%%%%%  
%%%%%%%%%%  package can run smoothly.                                  %%%%%%%%%%
%%%%%%%%%%  Suggest: every time you make some modification of codes    %%%%%%%%%%
%%%%%%%%%%  this pakcage, please run this script to debug. Notice that %%%%%%%%%%
%%%%%%%%%%  this script is only for testing so model parameters may not %%%%%%%%%
%%%%%%%%%%  be optimal.                                                %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  H.D. Li, lhdcsu@gmail.com                                  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+++ Cross validation
load corn_m51;
A=6;
K=5;
method='center';
N=500;
Nmcs=50;
CV=plscv(X,y,A,K,method)
MCCV=plsmccv(X,y,A,method,N)
DCV=plsdcv(X,y,A,K,method,Nmcs)
RDCV=plsrdcv(X,y,A,K,method,Nmcs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+++ build a model and make predictions on test set.
load corn_m51;
Rank=ks(X); %+++ Data partition using Kennard-Stone algorithm
Xcal=X(Rank(1:60),:);
ycal=y(Rank(1:60),:);
Xtest=X(Rank(61:80),:);
ytest=y(Rank(61:80),:);
PLS=pls(Xcal,ycal,10);  %+++ Build a PLS regression model using training set
[ypred,RMSEP]=plsval(PLS,Xtest,ytest); %+++ make predictions on test set
figure;
plot(ytest,ypred,'.',ytest,ytest,'r-');

% dados beer
[Xcal2,Xtest2] = pretrat(Xcal,Xtest,{'msc'});

nvl = 10;
CV=plscv(Xcal2,ycal,nvl,K,method)
plot(1:length(CV.RMSECV),CV.RMSECV,'b-',1:length(CV.RMSECV),CV.RMSECV,'ro')
nvl = 3;
PLS=pls(Xcal2,ycal,nvl);  %+++ Build a PLS regression model using training set
[ypred,RMSEP]=plsval(PLS,Xtest2,ytest); %+++ make predictions on test set
[fm.RMSEC,fm.R2c,fm.biasc]=rmse(ycal,PLS.y_fit,nvl);
[fm.RMSEP,fm.R2p,fm.biasp]=rmse(ytest,ypred)

figure(1)
plot(ycal,PLS.y_fit,'bo'); hold on ,
plot(ytest,ypred,'r*');refline(1,0)

figure(2)
subplot(211), plot(xaxis,Xcal)
subplot(212), plot(xaxis,PLS.VIP), hline(1,'g')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+++ Outlier detection
load corn_m51;
F=mcs(X,y,12,'center',1000,0.7)
figure;
plotmcs(F);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+++++++++++ CARS-PLS for variable selection  +++++++++++
load corn_m51;
A=10;
K=5; 
N=50;
method='center';
CARS=carspls(X,y,A,K,method,N,0,1);  % original version of CARS
figure;
plotcars(CARS);
CARS1=carspls(X,y,A,K,method,N,1,1);  % original version of CARS but with optimal LV selection using the Standard deviation information.
figure;
plotcars(CARS1);
CARS2=carspls(X,y,A,K,method,N,0,0);  % exactly reproducible version of CARS with random elements removed
figure;
plotcars(CARS2);
CARS3=carspls(X,y,A,K,method,N,1,0);  % exactly reproducible version of CARS but with optimal LV selection using the Standard deviation information.
figure;
plotcars(CARS3);


%+++++++++++ moving window PLS +++++++++++
[WP,RMSEF]=mwpls(X,y);
figure;
plot(WP,RMSEF);
xlabel('wavelength');
ylabel('RMSEF');

%+++ Random frog: here N is set to 1000 Only for testing.
%    N usually needs to be large, e.g. 10000 considering the huge variable
%    space.
Frog=randomfrog_pls(X,y,A,method,1000,5);
figure;
plot(Frog.probability);
xlabel('variable index');
ylabel('selection probability');
hold on
plot(X')

%+++++++++++ SPA  +++++++++++
% Subwindow Permutation Analysis for variable assessment.
F=spa(X,y,2,K,10,1000,0.75,method,0,1)
% [Q2in,Q2out]=plotspa(F,12)
nvl=6;
rmsecv_sel=[];
for ki=5:length(F.RankedVariable)
    varsel=F.RankedVariable(1:ki);
    CV=plscv(X(:,varsel),y,nvl,5,method);
    rmsecv_sel=[rmsecv_sel;[ki CV.optLV CV.RMSECV_min]];
end
plot(rmsecv_sel(:,1),rmsecv_sel(:,3),'b-',rmsecv_sel(:,1),rmsecv_sel(:,3),'ro')
xlabel('# variables'), ylabel('RMSECV')
% Escolher o RMSECV minimo
rmsecv_min=54;

VS = rmsecv_sel(find(rmsecv_sel(:,1)==rmsecv_min),:)
var_sel = F.RankedVariable(1:VS(1));
PLS=pls(X(:,var_sel),y,VS(2)); 


%+++++++++++ Verifica��o de n�o linearidade +++++++++++
load('corn_m51.mat')
Rank=ks(X); 
Xcal=X(Rank(1:60),:);
ycal=y(Rank(1:60),:);
Xtest=X(Rank(61:80),:);
ytest=y(Rank(61:80),:);
% valida��o cruzada
method='center';
nvl=10;
CV=plscv(Xcal,ycal,nvl,5,method)
plot(1:length(CV.RMSECV),CV.RMSECV,'b-',1:length(CV.RMSECV),CV.RMSECV,'ro')

% modelo PLS
nvl=4;
PLS=pls(Xcal,ycal,4)
[ypred,RMSEP]=plsval(PLS,Xtest,ytest); 
[fm.ypred,fm.RMSEP]=plsval(PLS,Xtest,ytest);
[fm.RMSEC,fm.R2c,fm.biasc]=rmse(ycal,PLS.y_fit,nvl);
[fm.RMSEP,fm.R2p,fm.biasp]=rmse(ytest,fm.ypred)
plot(ycal,PLS.y_fit,'bo'), hold on
plot(ytest,ypred,'r*');refline(1,0)
% figuras de merito
saida=fmerito(PLS,ycal,fm.ypred)
metodo={'center'}
confidence = confidence_interval(Xcal,Xtest,ycal,metodo,nvl,0.95)


% verifica��o de n�o linearidade presente no conjunto de dados.
%------------ APaRP ------------
T = PLS.X_scores;
yref = ycal;
AP = aparp(yref,T);


% Figuras de m�rito para a rotina libpls_1.95
nvl = 4;
PLS=pls(Xcal,ycal,nvl);  %+++ Build a PLS regression model using training set
[ypred,RMSEP]=plsval(PLS,Xtest,ytest); %+++ make predictions on test set
[fm.RMSEC,fm.R2c,fm.biasc]=rmse(ycal,PLS.y_fit,nvl);
[fm.RMSEP,fm.R2p,fm.biasp]=rmse(ytest,ypred)

saida=fmerito(PLS,ycal,ytest)

options=[]; 
options.Xpretreat      = {'mean'}; % pretratamento da matriz Xcal; ycal is always mean center;
options.vene           = 5;      % Cross validation: k from venetian blind k-fold;
options.alpha          = 0.95;   % scope level of confidence level.
options.xaxis          = [];     % vari�veis do eixo-X;
options.vl             = 4;     % n�mero de vari�veis latentes.
options.precisao       = [];     % matriz para o c�lculo de repetibilidade.
options.graficos       = 1;      % output graphic.end   
modelo = plsmodel(Xcal,ycal,Xtest,ytest,options)


%+++++++++++ nnls  +++++++++++
load('corn_m51.mat')
Rank=ks(X); 
Xcal=X(Rank(1:60),:);
ycal=y(Rank(1:60),:);
Xtest=X(Rank(61:80),:);
ytest=y(Rank(61:80),:);


XtX=Xcal'*Xcal;
Xty=Xcal'*ycal;

[x,w] = fnnls(XtX,Xty);

%%+++++++++++ Rotina PLS sub-Models ---------------------------------
load corn_m51;
Rank=ks(X); 
Xcal=X(Rank(1:60),:);
ycal=y(Rank(1:60),:);
Xtest=X(Rank(61:80),:);
ytest=y(Rank(61:80),:);

K=5;
method='center';
nvl = 10;
CV=plscv(Xcal,ycal,nvl,K,method)
plot(1:length(CV.RMSECV),CV.RMSECV,'b-',1:length(CV.RMSECV),CV.RMSECV,'ro')
nvl = 6;
PLS=pls(Xcal,ycal,nvl);  %+++ Build a PLS regression model using training set
[ypred,RMSEP]=plsval(PLS,Xtest,ytest); %+++ make predictions on test set
[fm.RMSEC,fm.R2c,fm.biasc]=rmse(ycal,PLS.y_fit,nvl);
[fm.RMSEP,fm.R2p,fm.biasp]=rmse(ytest,ypred)

SubGrupos = sub_grupos(Xcal,8,'euclidean');
SubModels=sub_pls_model(Xcal,ycal,SubGrupos,nvl,method);
Subpred=sub_pls_pred(SubModels,Xtest,ytest);

marc={'ko';'bo';'ro';'go';'co';'ks';'bs';'rs';'gs';'cs'};
for ki=1:size(SubModels.models,1)
    plot(SubModels.models{ki,2}.yc(:,1),SubModels.models{ki,2}.yc(:,2),marc{ki})
    hold on
end


%+++++++++++ Determina��o de outliers +++++++++++
% Determina a influencia de cada amostra na constru��o do modelo PLS
% Utilizando a rotina
load corn_m51;
method='center';
CV=plscv(X,y,10,5,method)
plot(1:length(CV.RMSECV),CV.RMSECV,'b-',1:length(CV.RMSECV),CV.RMSECV,'ro')

[Sample,Fator] = outlier_detec(X,y,6);


%+++++++++++ Sele��o de amostras +++++++++++
% Sele��o de amostras para conjuntos de calibra��o, previs�o e
% transfer�ncia
load('corn_m51.mat')
Rank=ks(X); 
Xcal=X(Rank(1:60),:);
ycal=y(Rank(1:60),:);
Xtest=X(Rank(61:80),:);
ytest=y(Rank(61:80),:);

X2=pretrat(X,[],{'msc';[9,2,1]});% preprocessar os dados
vd = dist_point(X2,mean(X2));    % determinar a dist�ncia ao espectro m�dio
[~,vdsamp]=sort(vd);             % ordemar as amostas.
samples=10;                      % N�mero de amostras a serem selecionadas
vd_ord=round(linspace(1,length(vdsamp),samples+2));
vd_ord=vd_ord(2:end-1);
vd_ord=vdsamp(vd_ord);

Xtransf=X(vd_ord,:);
ytrans=y(vd_ord,:);


%+++++++++++ Monte Carlo +++++++++++
% Utiliza��o do m�todo Monte Carlo na avalia��o da exatid�o de um modelo.
% ---------------- PLS ----------------
load('corn_m51.mat')
nvl=10;
method={'center';[7,2,1]};
N=500;
X2=pretrat(X,[],method);
MCCV=plsmccv(X2,y,nvl,'center',N)
plot(1:length(MCCV.RMSECV),MCCV.RMSECV,'b-',1:length(MCCV.RMSECV),MCCV.RMSECV,'ro')

nvl=6; % escolhido ap�s observa��o do gr�fico de RMSECV
nMC=1000;
ratio=70; % porcentagem da popula��o utilizada para treinamento 
mc.ratio=round(length(y)*ratio/100);
mc.Ypred=[];
mc.RMSEP=[];
mc.R2p=[];
ki=0;
while ki<nMC
    mc.random=randperm(length(y));
    Xcal=X(mc.random(1:mc.ratio),:);
    ycal=y(mc.random(1:mc.ratio),:);
    Xtest=X(mc.random(mc.ratio+1:end),:);
    ytest=y(mc.random(mc.ratio+1:end),:);
    mc.cal{ki+1,1}=mc.random(1:mc.ratio);
    mc.cal{ki+1,2}=mc.random(mc.ratio+1:end);
    % constru��o do modelo PLS
    [Xcal,Xtest]=pretrat(Xcal,Xtest,method);% pr�-processamento
    PLS=pls(Xcal,ycal,nvl,'center');  
    [ypred,RMSEP]=plsval(PLS,Xtest,ytest);
    [fm.RMSEP,fm.R2p]=rmse(ytest,ypred);
    mc.Ypred=[mc.Ypred;[ytest,ypred]];
    mc.RMSEP=[mc.RMSEP;fm.RMSEP];
    mc.R2p=[mc.R2p;fm.R2p];
    mc.cal{ki+1,3}=fm.RMSEP;
    mc.cal{ki+1,4}=fm.R2p;
    ki=ki+1;
end

figure(1) % frequ�ncia dos valores de RMSEP 
subplot(1,2,1), hist(mc.RMSEP,round(sqrt(nMC)))
xlabel('RMSEP','FontSize',12),ylabel('Frequency','FontSize',12)
subplot(1,2,2), hist(mc.R2p,round(sqrt(nMC)))
xlabel('R^2 prediction','FontSize',12),ylabel('Frequency','FontSize',12)

figure(2) % gr�ficos dos vaores previstos
plot(mc.Ypred(:,1),mc.Ypred(:,2),'k.')
refline(1,0)
xlabel('Reference','FontSize',12),ylabel('Prediction','FontSize',12)

faixa=[0.0352 0.0354];
samp.samp1=find(mc.RMSEP>=faixa(1));
samp.samp2=find(mc.RMSEP<=faixa(2));
samp.samp3=intersect(samp.samp1,samp.samp2);

%---------------------------------------
% ---------------- OPLS ----------------
load('corn_m51.mat')
% D:\Documentos\Filgueiras\Tool MatLab\libPLS_1.95\opls
cpnto=10;
K=5;
method={'center'};
X2=pretrat(X,[],method);
CVopls = oplscv(X2,y,K,cpnto);

cpnto=5; % escolhido ap�s observa��o do gr�fico de RMSECV
nMC=1000;
ratio=70; % porcentagem da popula��o utilizada para treinamento 
mc.ratio=round(length(y)*ratio/100);
mc.Ypred=[];
mc.RMSEP=[];
mc.R2p=[];
ki=0;
while ki<nMC
    mc.random=randperm(length(y));
    Xcal=X(mc.random(1:mc.ratio),:);
    ycal=y(mc.random(1:mc.ratio),:);
    Xtest=X(mc.random(mc.ratio+1:end),:);
    ytest=y(mc.random(mc.ratio+1:end),:);
    mc.cal{ki+1,1}=mc.random(1:mc.ratio);
    mc.cal{ki+1,2}=mc.random(mc.ratio+1:end);
    % constru��o do modelo OPLS
    [Xcal,Xtest]=pretrat(Xcal,Xtest,method);
    OPLSmodel = oplsmodel(Xcal,ycal,cpnto,Xtest,ytest,0);
    [fm.RMSEP,fm.R2p]=rmse(ytest,OPLSmodel.y_testpred);
    [RMSEP,R2p] = rmse(ytest,OPLSmodel.y_testpred);    
    mc.Ypred=[mc.Ypred;[ytest,OPLSmodel.y_testpred]];
    mc.RMSEP=[mc.RMSEP;fm.RMSEP];
    mc.R2p=[mc.R2p;fm.R2p];
    mc.cal{ki+1,3}=fm.RMSEP;
    mc.cal{ki+1,4}=fm.R2p;
    ki=ki+1;
end

figure(1) % frequ�ncia dos valores de RMSEP 
subplot(1,2,1), hist(mc.RMSEP,round(sqrt(nMC)))
xlabel('RMSEP','FontSize',12),ylabel('Frequency','FontSize',12)
subplot(1,2,2), hist(mc.R2p,round(sqrt(nMC)))
xlabel('R^2 prediction','FontSize',12),ylabel('Frequency','FontSize',12)

figure(2) % gr�ficos dos vaores previstos
plot(mc.Ypred(:,1),mc.Ypred(:,2),'k.')
refline(1,0)
xlabel('Reference','FontSize',12),ylabel('Prediction','FontSize',12)

faixa=[0.0352 0.0354];
samp.samp1=find(mc.RMSEP>=faixa(1));
samp.samp2=find(mc.RMSEP<=faixa(2));
samp.samp3=intersect(samp.samp1,samp.samp2);


% Dados de N total
load('dados.mat')

options.pretrat={'center'};
options.nvl = 8; % N�mero de vari�veis latentes.
options.nMC = 100; % n�mero de corridas doo algoritmo Monte Carlo
options.ratio=70; % Porcentagem das amostras de calibra��o utilizadas para constru��o do modelo.
options.alpha=0.95; % scoe level confidence

modelo = pls_ensemble(Xcal,ycal,Xtest,ytest,options)


