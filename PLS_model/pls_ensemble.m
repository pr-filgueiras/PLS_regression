function modelo = pls_ensemble(X,y,Xtest,ytest,options)
% modelo = pls_ensemble(X,y,Xtest,ytest,options)
% Ensemble PLS
% Rotina para construção do modelo ensemble PLS.
% input
%
% options.pretrat={'center'};
% options.nvl = 6; % Número de variáveis latentes.
% options.nMC = 100; % número de corridas doo algoritmo Monte Carlo
% options.ratio=70 % Porcentagem das amostras de calibração utilizadas para construção do modelo.
% options.alpha=0.95 % scoe level confidence
%
% Paulo R. Filgueiras  - 19/04/2016
%

% load('corn_m51.mat')
% nvl=10;
% method={'center';[7,2,1]};
% N=500;
% X2=pretrat(X,[],method);
% MCCV=plsmccv(X2,y,nvl,'center',N)
% plot(1:length(MCCV.RMSECV),MCCV.RMSECV,'b-',1:length(MCCV.RMSECV),MCCV.RMSECV,'ro')




% nvl=6; % escolhido após observação do gráfico de RMSECV
% nMC=1000;
% ratio=70; % porcentagem da população utilizada para treinamento 

mc.ratio=round(length(y)*options.ratio/100);
mc.ycalt=[];
mc.Ypred=[];
mc.RMSEP=[];
mc.R2p=[];


ki=0;
while ki<options.nMC
    mc.random=randperm(length(y));
    Xcal=X(mc.random(1:mc.ratio),:);
    ycal=y(mc.random(1:mc.ratio),:);
    Xcalt=X(mc.random(mc.ratio+1:end),:);
    ycalt=y(mc.random(mc.ratio+1:end),:);
    
    mc.cal{ki+1,1}=mc.random(1:mc.ratio);
    mc.cal{ki+1,2}=mc.random(mc.ratio+1:end);
    % construção do modelo PLS
    [~,Xcalt]=pretrat(Xcal,Xcalt,options.pretrat);% pré-processamento de calibração
    [Xcal,Xtest2]=pretrat(Xcal,Xtest,options.pretrat);% pré-processamento de teste
    PLS=pls(Xcal,ycal,options.nvl,'center');  % Construção do sub-modelo
    [ypred,RMSEP]=plsval(PLS,Xtest2,ytest);  % Aplicação nas amostras de teste
    ypredcal=plsval(PLS,Xcalt,ycalt);  % Aplicação nas amostras de calibração  
    
    mc.ycalt=[mc.ycalt;[ycalt,ypredcal]];
    [fm.RMSEP,fm.R2p]=rmse(ytest,ypred);
    mc.Ypred=[mc.Ypred,ypred];
    mc.RMSEP=[mc.RMSEP;fm.RMSEP];
    mc.R2p=[mc.R2p;fm.R2p];
    mc.cal{ki+1,3}=fm.RMSEP;
    mc.cal{ki+1,4}=fm.R2p;
    ki=ki+1;
end

modelo.data=date;
modelo.mc=mc;
modelo.predtest=[ytest,mean(mc.Ypred,2),std(mc.Ypred,0,2)];

modelo.nu=length(ytest)-options.nvl-1;
modelo.T=tinv((1-((1-options.alpha)/2)),modelo.nu);
modelo.IC=(modelo.T*(1/sqrt(length(ytest)))).*modelo.predtest(:,3);


figure(1) % frequência dos valores de RMSEP 
subplot(1,2,1), hist(mc.RMSEP)
xlabel('RMSEP','FontSize',12),ylabel('Frequency','FontSize',12)
subplot(1,2,2), hist(mc.R2p)
xlabel('R^2 prediction','FontSize',12),ylabel('Frequency','FontSize',12)

figure(2) % gráficos dos vaores previstos
plot(modelo.predtest(:,1),modelo.predtest(:,2),'bo')
refline(1,0)
xlabel('Reference','FontSize',12),ylabel('Prediction','FontSize',12)

% faixa=[0.0352 0.0354];
% samp.samp1=find(mc.RMSEP>=faixa(1));
% samp.samp2=find(mc.RMSEP<=faixa(2));
% samp.samp3=intersect(samp.samp1,samp.samp2);







