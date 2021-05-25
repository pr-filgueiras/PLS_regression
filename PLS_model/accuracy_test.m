function [pvalue,dist_tt,meandiff] = accuracy_test(y,yhatA,yhatB,teste,niter,alpha)
% Teste rand�mino para compara��o da acur�cias de dois modelos de
% calibra��o multivariada.
% input
%   y   : valores y de refer�ncia (do vetor de teste);
% yhatA : valores de y previsto por um dos modelos;
% yhatB : valores de y previsto pelo outro modelo;
% teste : tipo de teste:
%    randbi  : teste randomico bicaudal.
%    randuni : teste randomico unicaudal. yhatA > yhatB
%  tpareado  : teste-t pareado bicaudal.
% niter : n�mero de permuta��es;
% alpha : n�vel de signific�ncia adotado;
% 
% output
% pvalue : estat�stica de teste.
% Se o pvalue < alpha : a diferen�a na acur�cia dos modelos � diferente.
%
% Paulo R. Filgueiras  - 26/06/2014
% 
% Refer�ncia
% Van der Voet, H. Chemom. Intel. Lab. Syst. 25, 1994, 313-323.
%
if nargin==3 , teste='randbi' ;niter=100000; alpha=0.05; 
elseif nargin==4, niter=100000; alpha=0.05; 
elseif nargin==5 , alpha=0.05; end
    
eA=y-yhatA;
eB=y-yhatB;
diff=eA.^2-eB.^2;
meandiff=mean(diff);
n=length(diff);
%niter=199;
sum=0;
dist_tt=[];
if strcmp(teste,'randbi');   % teste ramd�mico bicaudal
for k=1:niter
    randomsign=2*round(rand(n,1))-1;
    signeddiff=randomsign.*diff;
    meansigneddiff=mean(signeddiff);
    dist_tt=[dist_tt;meansigneddiff];
    sum=sum+(abs(meansigneddiff)>=abs(meandiff));
end
pvalue=(sum+1)/(niter+1);

disp('  ')
if pvalue < alpha
    s = sprintf('tcal = %g < alpha = %g',pvalue,alpha); disp(s)
    disp('Modelos com DIFEREN�AS na acur�cia')
else
    s = sprintf('tcal = %g > alpha = %g',pvalue,alpha); disp(s)
    disp('Modelos com acur�cias iguais')
end
disp('  ')

figure(1);
axes('FontSize',16,'FontName','Arial');
hist(dist_tt,50), vline(meandiff,'r');
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[.8 .8 .8],'EdgeColor','k')
ylabel('Frequ�ncia','FontSize',24,'FontName','arial');
xlabel('Distribui��o aleat�ria','FontSize',24,'FontName','arial');  

elseif strcmp(teste,'randuni');    % teste ramd�mico unicaudal
for k=1:niter
    randomsign=2*round(rand(n,1))-1;
    signeddiff=randomsign.*diff;
    meansigneddiff=mean(signeddiff);
    dist_tt=[dist_tt;meansigneddiff];
    sum=sum+(meansigneddiff>=meandiff);
end
pvalue=(sum+1)/(niter+1);

disp('  ')
if pvalue < alpha
    s = sprintf('p-cal = %g < alpha = %g',pvalue,alpha); disp(s)
    disp('Modelos com DIFEREN�AS na acur�cia')
else
    s = sprintf('p-cal = %g > alpha = %g',pvalue,alpha); disp(s)
    disp('Modelos com acur�cias iguais')
end
disp('  ') 
    
figure(1);
axes('FontSize',16,'FontName','Arial');
hist(dist_tt,50), vline(meandiff,'r');
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[.8 .8 .8],'EdgeColor','k')
ylabel('Frequ�ncia','FontSize',24,'FontName','arial');
xlabel('Distribui��o aleat�ria','FontSize',24,'FontName','arial');  

elseif strcmp(teste,'tpareado');
% teste-t em para para m�dias (teste bicaudal)
diff1 = eA-eB;
meandiff1=mean(diff1);
stddiff1=std(diff1);
tcalc = (meandiff1*sqrt(length(diff1)))/(stddiff1);
tstat = tinv((1-alpha/2),length(diff1)-1);
pvalue = 1-tcdf(tcalc,length(diff1)-1);
disp(' ')
disp(' Teste-t em para para m�dias ')
disp(' ')
if tcalc < tstat
    s = sprintf('tcal = %g < t-tabelado = %g',tcalc,tstat); disp(s)
    disp('Modelos com acur�cias iguais')
    s = sprintf('pvalor = %g',pvalue); disp(s)
else
    s = sprintf('tcal = %g > t-tabelado = %g',tcalc,tstat); disp(s)
    disp('Modelos com DIFEREN�AS na acur�cia')
    s = sprintf('pvalor = %g',pvalue); disp(s)
end
disp('  ') 

end
