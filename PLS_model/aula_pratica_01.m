%% Aula Prática - Calibração Multivariada

% dados corn_m51
load('corn_m51.mat')
plot(X')

[objetos,Xcal,Xtest,ycal,ytest]=caltest(X,y,70,'k',0,{'center'});

nvl = 10;
kfold = 5;
method = 'center';

CV = plscv(Xcal,ycal,nvl,kfold,method);
plot(CV.RMSECV),hold on , plot(CV.RMSECV,'ro')
xlabel('# latent variable')
ylabel('RMSECV')

nvl = 4;

modelo=pls(Xcal,ycal,nvl,method);  
[ypred,RMSEP]=plsval(modelo,Xtest,ytest); 
[RMSEC,R2c,biasc] = rmse(ycal,modelo.y_fit,nvl)
[RMSEP,R2p,biasp] = rmse(ytest,ypred)

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


%% Dados beer
load('nirbeer.mat')

nvl = 10;
kfold = 5;
method = 'center';

CV = plscv(Xcal,ycal,nvl,kfold,method);
plot(CV.RMSECV),hold on , plot(CV.RMSECV,'ro')
xlabel('# latent variable')
ylabel('RMSECV')

nvl = 3;

modelo=pls(Xcal,ycal,nvl,method);  
[ypred,RMSEP]=plsval(modelo,Xprev,yprev); 
[RMSEC,R2c,biasc] = rmse(ycal,modelo.y_fit,nvl)
[RMSEP,R2p,biasp] = rmse(yprev,ypred)

%+++ Gráfico Referência vs Previsto
plot(ycal,modelo.y_fit,'bo'); hold on
plot(yprev,ypred,'r*')
refline(1,0)
legend('Calibration','Prediction')
xlabel('Reference')
ylabel('Predicted')

subplot(2,1,1), plot(Xcal')
subplot(2,1,2), plot(modelo.VIP)
xlabel('Wavelength nm')

