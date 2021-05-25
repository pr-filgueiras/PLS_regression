function perm = perm_teste(ymeasured,ypred,k,n)
% teste de permutação - avaliação de resíduos distribuidos segundo uma
% tendência quadrática. E(y) = b0 + b1*y + b2*y^2
% O teste é baseado na relação entre E(y) e y, dado de forma quadrática.
% H0 : b2 >= 0 - não há relação quadratica entre o E(y) e y.
% Ha : b2 < 0 - há evidências de relação quadrática entre E(y) e y.
%
% entradas:
% ymeasured : vector with values measured by reference method;
%  ypred    : vector with values predicted by model;  
%  k        : number of permutation for the test;
%  n        : grau of polinômio.
%
%  saída  :
%   perm  : estrutura com 
%    teste - relefente a estatistica de teste e;
%  permute - resposta do teste obtidos com valores permutados;
%  p_valor - valor de significância dos teste.
% 
% Publicação:
% Filgueiras PR, et al. Evaluation of trends in residuals of multivariate calibration models by permutation test.
% Chemometrics and Intelligent Laboratory Systems. Volume 133, 15 April 2014, Pages 33-41. 
% https://doi.org/10.1016/j.chemolab.2014.02.002
%
% Prof. Paulo Roberto Filgueiras (15/03/2014)
%

%% cálculo da estatística de teste
if nargin==3,n=2;end
b2_star=polyfit(ypred,[ymeasured-ypred],n);
H=waitbar(0,'Running...','name','Permutation test... ');
b2_perm=[];
for ki=1:k
    permut=randperm(length(ypred)); 
    y_perm=ypred(permut); 
    b2_perm=[b2_perm;polyfit(y_perm,[ymeasured-ypred],n)];
    waitbar(ki/k);
end
close(H);
if b2_star(1)>0
    pvalue=1-length(find(b2_perm(:,1)<=b2_star(1)))/k;
else
    pvalue=1-length(find(b2_perm(:,1)>=b2_star(1)))/k;
end
perm.b2_star=b2_star;
perm.b2_perm=b2_perm(:,1);
perm.pvalue=pvalue;
figure(1);
axes('FontSize',16,'FontName','Arial');
hist(b2_perm(:,1),50); vline(b2_star(1),'r');
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[.8 .8 .8],'EdgeColor','k')
ylabel('Frequência','FontSize',24,'FontName','arial');
xlabel('Polynomial Coeficiente','FontSize',24,'FontName','arial'); 



%% t
%n=length(y); % comprimento de vetor y.
%y0 = ones(n,1); y2=y.^2;
%ytest = [y0 y y2]; 
%est_teste = pinv(ytest)*erro; % modelo linear quadrático.
%perm.teste=est_teste(3,1);
%perm.p_valor = 0;
%% permutações
%for i=1:k
%    aperm=randperm(n); 
%    yperm=y(aperm); % valores de y permutado
%    y2=yperm.^2;
%    yrand = [y0 yperm y2]; 
%    rand_teste = pinv(yrand)*erro; % modelo linear quadrático.
%    perm.permute(i,1)=rand_teste(3,1);
%    
%    if perm.permute(i,1)<=perm.teste
%        perm.p_valor=perm.p_valor+1;
%    end
%end

%perm.p_valor2=k-perm.p_valor;
%perm.p_valor=perm.p_valor/k; % p-valor para b2<0
%if perm.teste>0, perm.p_valor=1-perm.p_valor, end
%% histograma dos resultados randomicos
%figure(1);
%axes('FontSize',16,'FontName','Arial');
%histfit(perm.permute,50); vline(perm.teste,'r');%vline(-perm.teste,'r');
%h = findobj(gca,'Type','patch');
%set(h,'FaceColor',[.8 .8 .8],'EdgeColor','k')
%ylabel('Frequência','FontSize',24,'FontName','arial');
%xlabel('Coeficiente (b_2)','FontSize',24,'FontName','arial'); 

