function saida=fmerito(modelo,ycal,yprev,RMSEP)

%------------------------------------------------------------------------------------------
%    Este programa determina o net analyte signal para uma matriz de espectros 
% por pls, tendo como algoritmo base o NIPALS. Os dados X ja devem estar
% preprocessados.
% 
% Sintaxe:  saida=fmerito(modelo,ycal)
%
% Entradas:
% modelo = modelo PLS-Chemometrics
% ycal   = vetor y de calibração
% yprev  = vetor y de teste
% RMSEP  = RMSEP do modelo
%
% Saidas:
%      
%      nas = matriz que contém os vetores nas para as I amostras.
%      
%      nnas = valor escalar do nas
%      
%      corre_H = fator de correçao para o calculo do vetor nas para amostras futuras.
%          Para um conjunto de amostras desconhecidas(Xv) 
%      
%      seletividade = vetor que conten a seletividade de cada amostra. 
%
%      sensibilidade = sensibilidade do modelo.
%      
%      sen_analitica = sensibilidade analítica do modelo.
%
%      ruido = ruido espectral
%
%      sinal_ruido = razão sinal ruído.
%
%      ld = limite de detecção do modelo.
%
%      lq = limite de quantificação.
%      
% Exemplo:
% saida=fmerito(PLS,ycal,yprev, RMSEP)
%
%---------------------------------------------------------------------------------------------
x = modelo.X_pretreat;
y = ycal;
b = modelo.regcoef_pretreat;
T = modelo.X_scores;
P = modelo.X_loadings;
vl = size(modelo.X_scores,2);


[L1,C1]=size(x);
[Rc,mx]=center(x,1);
[cc,my]=center(y,1);
cest=Rc*b;
% calculo do nnas_c
nnas_c=cest/norm(b);
Z(1:L1,1)=1;
Z(1:L1,2)=nnas_c;
ccc=cc(:,1)+ones(L1,1)*my;
cte=pinv(Z)*ccc;
constante=-cte(1,1)/cte(2,1);
nnas_cm(:,1)=nnas_c(:,1)-constante;
   
% Fator de correçao 
H=b/norm(b);

% Sensibilidade(escalar)
sen=1/norm(b);

% seletividade
for i=1:L1
    sel=nnas_cm./norm(x(i,:));
end

% Vetor nas para as I amostras de R
nas_c=nnas_cm*norm(b)*pinv(b);

% Vetor de sensibilidade
for i=1:L1
    s(i,:)=nnas_c(i,:)./cest(i,1);
end


%saidas
saida.nas=nas_c;
saida.nnas=nnas_cm;
saida.corre_H=H;
saida.seletividade=sel;
saida.sensibilidade=sen;
saida.sen_analitica = [];
saida.inv_sen_analitiva =[];

xp=T*P'; 
xpm=xp+ones(size(x,1),1)*mx;
er=x-xpm;
va=var(er);
vva=mean(va); %média da variancia do ruido
dx=sqrt(vva); % ruido espectral

%++++++++++ Razão Sinal/Ruido ++++++++++
sr=saida.nnas./dx;

%++++++++++ Sensibilidade analitica ++++++++++
sena=sen/dx; %  (SENSIBILIDADE ANALÍTICA) = 410.5

%++++++++++ Limite de detecção ++++++++++
ld=3*dx*(1/sen);
% ld=0.0073

%++++++++++ Limite de quantificação ++++++++++
lq=10*dx*(1/sen);
% lq=0.0244

saida.sen_analitica = sena;
saida.inv_sen_analitiva = 1/sena;
saida.ruido_espectral = dx;
saida.sinal_ruido = sr;
if min(sr)<0; saida.sinal_ruido_min =0; else saida.sinal_ruido_min = min(sr);end
saida.sinal_ruido_max = max(sr);
saida.LD = ld;
saida.LQ = lq;
saida.RPD = std(yprev)/RMSEP;
saida.RPIQ = iqr(yprev)/RMSEP;

ajuste = polyfit(modelo.y_fit,ycal,1);
saida.ajuste.intercepto = ajuste(2);
saida.ajuste.inclinacao = ajuste(1);
saida.ajuste.Q2 = corrcoef(modelo.y_fit,ycal);  saida.ajuste.Q2=saida.ajuste.Q2(1,2)^2; 

yc_nas = saida.nas*b;
para1 = mean(ycal);     yc_nas = yc_nas+para1;
ajuste_nas = polyfit(yc_nas,ycal,1);
saida.ajuste_nas.intercepto = ajuste_nas(2);
saida.ajuste_nas.inclinacao = ajuste_nas(1);
saida.ajuste_nas.Q2 = corrcoef(yc_nas,ycal);  saida.ajuste_nas.Q2=saida.ajuste_nas.Q2(1,2)^2; 

% if nargin>2
%     [Rp,mx]=center(modelo.datap.prev,1);
%     [cpp,my]=center(yprev,1);
%     L2=size(Rp,1);
%     cest=Rp*b;
%     nnas_p=cest/norm(b);
%     Zp(1:L2,1)=1;
%     Zp(1:L2,2)=nnas_p;
%     cccp=cpp(:,1)+ones(L2,1)*my;
%     cte=pinv(Zp)*cccp;
%     constante=-cte(1,1)/cte(2,1);
%     nnas_pm(:,1)=nnas_p(:,1)-constante;
%     saida.nas_p=nnas_pm*norm(b)*pinv(b);
%     saida.nnas_p=nnas_pm;
% end




%% Sub-function

function [cdata,me,cnewdata]=center(data,opt,newdata)

%#										
%#  function [cdata,me,ctest]=center(data,opt,newdata);			
%#										
%#  AIM: 	Centering along columns, rows or double centering		
%#										
%#  PRINCIPLE:  Removal of the column-, row- or overall mean from 		
%#              each column, row or both, respectively 		 		
%# 				 If a test data set is available it can ONLY be 
%#              column centered using the mean from the calibration
%#              data set.
%#
%#
%#  INPUT:	data: (m x n) matrix with m rows and n variables		
%#				opt: optional							
%#		     1 = column centering					
%#		     2 = row centering						
%#		     3 = double centering					
%#          newdata: (mt x n) test matrix with mt rows and n variables
%#					
%#			 							
%#  OUTPUT:	cdata: (m x n) matrix containing centered data			
%#				me: mean vector, overall mean (scalar)
%#              newdata: (mt*n) test matrix centered with the mean of data
%#	
%#										
%#  AUTHOR: 	Andrea Candolfi				 			
%#	    			Copyright(c) 1997 for ChemoAc					
%#          	FABI, Vrije Universiteit Brussel            			
%#          	Laarbeeklaan 103 1090 Jette					
%#   										
%# VERSION: 1.2 (25/02/2002)							
%#										
%#  TEST:   	I. Stanimirova	& S. Gourvénec & M. Zhang
%#										
	

[m,n]=size(data);

if nargin==1;
  opt=[4];
  while opt>3 || opt<=0 
    opt=input('column centering(1), row centering(2), double centering(3):');
  end
end


if opt==1			% column centering 
   me=mean(data);
   cdata=data-ones(m,1)*me;
end

if opt==2			% row centering
   me=mean(data')';
   cdata=data-me*ones(1,n);
end

if opt==3 	% double centering
   me=mean(mean(data));
   mej=mean(data');
   mei=mean(data);
   cdata=data-(ones(m,1)*mei)-(ones(n,1)*mej)'+(ones(m,n)*me);
end

if exist('newdata')==1			% center new data
    [mt,n]=size(newdata);
    
    if opt==1				% column centering 
        me=mean(data);
        cnewdata=newdata-ones(mt,1)*me;
    else
        error('Row centering and double centering are impossible to perform on a test set');
    end
    
end







