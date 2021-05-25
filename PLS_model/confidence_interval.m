function confidence = confidence_interval(Xcal,Xtest,ycal,Xpretreat,vl,alpha)
% confidence = confidence_interval(Xcal,Xtest,ycal,Xpretreat,vl,alpha)
%
% Rotina para calcular o intervalo de confiança em modelos de regressão PLS.
% Os cálculos são baseados na norma ASTM E1655-05 (2012).
% NORMA: ASTM E1655-05 (2012) Standard Practices for Infrared Multivariate 
%        Quantitative Analysis.
% input
% Xcal  : Matriz X de calibração.
% Xtest : Matriz X de teste.
% ycal  : vetor y de calibração.
% Xpretreat : método de preprocessamento dos dados X. 
%             Xcal e Xtest serão pre-processados conforme construção do 
%             modelo de calibração. Ex: {'center'}
%  vl  : escalar com o número de variáveis latentes.
% alpha : nível de abrangência do intervalo de confiança. default: 0.05
%
% output
% confidence : estrutura contendo o intervalo de confiança e o leverage das
%              amostras
% saidas:
%   + cal : matriz contendo [valor previsto(PLS) intervalo de confiança] calibração
%   + test: matriz contendo [valor previsto(PLS) intervalo de confiança]  teste.
%   + leverage_cal  : vetor contendo o leverage das amostras de calibração.
%   + leverage_test : vetor contendo o leverage das amostras de teste.
%   + leverage_limit : valor limite de leverage segundo a norma ASTM e 1655-05(2012).
%
% Exemplo:
% confidence = confidence_interval(Xcal,Xtest,ycal,{'msc';'center'},3,0.95);
%
% Paulo R. Filgueiras  - 11/09/2014
%

if nargin<6, alpha=0.05; end

% -------------------- Cálculo dos valores previstos --------------------
[ycal_prev,coef] = submodel_pls(Xcal,ycal,Xcal,Xpretreat,vl,ycal); %previsão das amostras de calibração.
ytest_prev       = submodel_pls(Xcal,ycal,Xtest,Xpretreat,vl,ones(size(Xtest,1),1)); %previsão das amostras de teste.


% -------------------- Standard Error of Calibration  SEC --------------------
RMSEC = sqrt((sum((ycal-ycal_prev).^2))/(length(ycal)-vl-1));
% norma ASTM E1655-05 (2012) pag.12

% --------------------  leverage  --------------------
[~,Xtest_pretreat]=pretrat(Xcal,Xtest,Xpretreat);  

tv=Xtest_pretreat*coef.R(:,1:vl)*inv(coef.P(:,1:vl)'*coef.R(:,1:vl));

h_limit=3*(vl)/size(Xcal,1); % leverage limite, segundo a norma ASTM E1655-05 (2012) pag.15

% Leverage from calibration set
for ki=1:size(Xcal,1)
    h_cal(ki,1)=coef.T(ki,1:vl)*inv(coef.T(:,1:vl)'*coef.T(:,1:vl))*coef.T(ki,1:vl)';
end

% Leverage from test set 
for ki=1:size(Xtest,1)
    h_test(ki,1)=tv(ki,1:vl)*inv(coef.T(:,1:vl)'*coef.T(:,1:vl))*tv(ki,1:vl)';
end
% h_limit : leverage limite
%  h_cal  : leverage das amsostras de calibração.
%  h_test : leverage das amsostras de teste.
% --------------------  end leverage computation  --------------------


% --------------------  Estatística t-student  --------------------
t_stat=abs(tinv((alpha)/2,length(ycal)-1));


% --------------------  Intervalo de confiança  --------------------
% Calculado segundo a norma ASTM E1655-05 (2012) pag.14
confidence.cal  =[ycal_prev  t_stat*RMSEC*sqrt(1+h_cal)];
confidence.test =[ytest_prev t_stat*RMSEC*sqrt(1+h_test)];
confidence.leverage_cal = h_cal;
confidence.leverage_test = h_test;
confidence.leverage_limit = h_limit;

%if isvector(ytest)&& modelo.options.graficos~=0
%errorbar(confidence.cal(:,1),confidence.cal(:,2),confidence.cal(:,3),'ko'), hold on
%errorbar(confidence.test(:,1),confidence.test(:,2),confidence.test(:,3),'r*'), refline(1,0)
%xlabel('Reference','fontsize',14);ylabel('Predicted','fontsize',14);
%legend('Calibration','Prediction');
%end


function [yprev,coef] = submodel_pls(Xcal_sub,ycal_sub,Xtest_sub,Xpretreat,vl,ytest)
% sub-rotina para construção do modelo PLS.
% output:
%  yprev : valor de y previsto.
%  coef  : coeficientes de regressão do algoritmo SIMPLS
%
% pre-processamento dos dados X e y
% [Xcal_sub,Xtest_sub]=pretrat(Xcal_sub,Xtest_sub,Xpretreat);      
% [ycal_sub,para1,para2]=pretreat(ycal_sub,'center');
% determinação dos coeficientes de regressão do modelos PLS.
% [coef.B,coef.C,coef.P,coef.T,coef.U,coef.R,coef.R2X,coef.R2Y]=plssim(Xcal_sub,ycal_sub,vl);
% yprev=Xtest_sub*coef.B;
% yprev=yprev+para1; % valor de y previsto.

[Xcal_sub,Xtest_sub] = pretrat(Xcal_sub,Xtest_sub,Xpretreat);
PLS=pls(Xcal_sub,ycal_sub,vl,'center');
yprev=plsval(PLS,Xtest_sub,ytest); 
coef.B=PLS.regcoef_pretreat;
%coef.C=PLS.;
coef.P=PLS.X_loadings;
coef.T=PLS.X_scores;
%coef.U=PLS.;
coef.R=PLS.W;
coef.R2X=PLS.R2X;
coef.R2Y=PLS.R2Y;



