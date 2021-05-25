function AP = aparp(yref,T)
% Método APaRP para determinação de não linearidade presente em um
% conjunto de dados. O algoritmo também aplica o teste de 
% Durbin-Watson (test for autocorrelation in linear regression)
% input
% T : Matriz de escores do modelo PLS ou PCR. T = PLS.X_scores;
% yref: vetor com os valores de referência.   yref = ycal;
% output
% AP : respostas do método APaRP e do teste Durbin-Watson 
% 
% Paulo R. Filgueiras  - 09/12/2015 
%
% Referência
% [V. Centner, O.E. de Noord, D.L. Massart, Detection of nonlinearity in 
%         multivariate calibration, Anal. Chim. Acta 376 (1998) 153-168]

score2 = [ones(size(T,1),1), T, T(:,1).^2]; % APaRP

b_APaRP=pinv(score2)*yref;
ypred_APaRP=score2*b_APaRP;
err_APaRP=yref-ypred_APaRP;
APaRP2=err_APaRP+b_APaRP(2)*score2(:,2)+b_APaRP(end)*score2(:,end);

plot(score2(:,2),APaRP2,'k+','Linewidth',1)
title('APaRP','FontSize',16)
xlabel('PC_1','FontSize',14)
ylabel('e_A_P_a_R_P + b_1*PC_1+b_1_1*PC_1^2','FontSize',14)
p1=polyfit(score2(:,2),APaRP2,1); refline(p1(1),p1(2));

[p,dw]=dwtest(err_APaRP,score2,'Tail','right');

AP.b_APaRP=b_APaRP;
AP.ypred_APaRP=ypred_APaRP;
AP.err_APaRP=err_APaRP;
AP.APaRP2=APaRP2;
AP.p=p;
AP.dw=dw;


