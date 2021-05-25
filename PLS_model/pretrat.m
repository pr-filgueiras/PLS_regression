function [Xp,Xtp]=pretrat(X,Xt,method)
    % rotina para preprocessar dados;
    % rotinas auxiliares: savgol; spdiags.
    % input:
    % X    : matriz de referência para preprocessamento;
    % Xt   : matriz de teste;
    % method : métodos para  preprodessar. É uma célula com os vários métodos de preprocessamento
    %        Ex: method={'center';'auto';'snv';'msc','pareto';'minmax';'minmax2';'norm';'deriv',[7,2,1]}
    %
    % output:
    % Xp   : matriz X preprocessada;
    % Xtp  : matriz Xt preprocessada;
    %
    % Exemplo: method={'msc';'deriv';[7,2,1];'minmax';}
    %          [Xp,Xtp]=pretrat(Xcal,Xprev,method);
    %
    % Paulo R. Filgueiras   -  27/02/2014
    %

[n1,n2]=size(method);
if n2>1; method=method'; end
npret=size(method,1);  % Número de preprocessamento    
for ki=1:npret
    if strcmp(method{ki,1},'center');[X,Xt]=center(X,Xt);
    elseif strcmp(method{ki,1},'auto');[X,Xt]=auto(X,Xt);    
    elseif strcmp(method{ki,1},'snv');[X,Xt]=snv(X,Xt);
    elseif strcmp(method{ki,1},'msc');[X,Xt]=msc(X,Xt);
    elseif strcmp(method{ki,1},'pareto');[X,Xt]=pareto(X,Xt);
    elseif strcmp(method{ki,1},'minmax');[X,Xt]=minmax(X,Xt);
    elseif strcmp(method{ki,1},'minmax2');[X,Xt]=minmax2(X,Xt);
    elseif strcmp(method{ki,1},'deriv');[X,Xt]=derivada(X,Xt,method{ki+1,1});   
    elseif strcmp(method{ki,1},'norm');[X,Xt]=normalizar(X,Xt); 
    elseif strcmp(method{ki,1},'emsc');[X,Xt]=emsc(X,Xt,method{ki+1,1}); 
    elseif strcmp(method{ki,1},'fir');[X,Xt]=fir(X,Xt,method{ki+1,1}); 
    elseif strcmp(method{ki,1},'osc');[X,Xt]=osc(X,Xt,method{ki+1,1},method{ki+2,1}); 
    else ; X=X;Xt=Xt;    
    end
end    
Xp=X;Xtp=Xt;  
    

function [X1,X1p] = center(Xe,Xet)
    % centrar na média
    % input: Xe : matriz X para centrar na média;
    %        Xet: matriz X de teste;
    % output: X1 : matriz Xe centrada na média;
    %         X1p: matriz Xet centrada na média;
    X1=Xe-ones(size(Xe,1),1)*mean(Xe);
    if size(Xe,2)==size(Xet,2)
        X1p=Xet-ones(size(Xet,1),1)*mean(Xe);
    else
        X1p=Xet;
    end
    
function [X1,X1p] = auto(Xe,Xet)
    % autoescalar
    % input: Xe : matriz X para autoescalar;
    %        Xet: matriz X de teste;
    % output: X1 : matriz Xe autoescalada;
    %         X1p: matriz Xet autoescalada;
    X1=(Xe-ones(size(Xe,1),1)*mean(Xe))./(ones(size(Xe,1),1)*std(Xe));
    if size(Xe,2)==size(Xet,2)
        X1p=(Xet-ones(size(Xet,1),1)*mean(Xe))./(ones(size(Xet,1),1)*std(Xe));
    else
        X1p=Xet;
    end
    
function [X1,X1p] = snv(Xe,Xet)
    % suavização normal padrão
    % input: Xe : matriz X para suavizar;
    %        Xet: matriz X de teste;
    % output: X1 : matriz Xe suavizada;
    %         X1p: matriz Xet suavizada;
    X1=(Xe-mean(Xe')'*ones(1,size(Xe,2)))./(std(Xe')'*ones(1,size(Xe,2)));
    if size(Xe,2)==size(Xet,2)
        X1p=(Xet-mean(Xet')'*ones(1,size(Xet,2)))./(std(Xet')'*ones(1,size(Xet,2)));
    else
        X1p=Xet;
    end
    
function [X1,X1p] = msc(Xe,Xet)
    % correção de espalhamento multiplicativo
    % input: Xe : matriz X para correçao;
    %        Xet: matriz X de teste;
    % output: X1 : matriz Xe corrigida;
    %         X1p: matriz Xet corrigida;
    para1=mean(Xe);
    for k2=1:size(Xe,1)
        coef=polyfit(para1,Xe(k2,:),1);
        X1(k2,:)=(Xe(k2,:)-coef(2))/coef(1);
    end
    % corrigindo a matriz de teste
    if size(Xe,2)==size(Xet,2)
    for k3=1:size(Xet,1)
        coef2=polyfit(para1,Xet(k3,:),1);
        X1p(k3,:)=(Xet(k3,:)-coef2(2))/coef2(1);
    end
    else
        X1p=Xet;
    end

function [X1,X1p] = pareto(Xe,Xet)
    % Preprocessamento pareto
    % input: Xe : matriz X para ser preprocessada;
    %        Xet: matriz X de teste;
    % output: X1 : matriz Xe preprocessada;
    %         X1p: matriz Xet preprocessada;
    X1=(Xe-ones(size(Xe,1),1)*mean(Xe))./(ones(size(Xe,1),1)*sqrt(std(Xe)));
    if size(Xe,2)==size(Xet,2)
        X1p=(Xet-ones(size(Xet,1),1)*mean(Xe))./(ones(size(Xet,1),1)*sqrt(std(Xe)));
    else
        X1p=Xet;
    end

function [X1,X1p] = minmax(Xe,Xet)
    % Preprocessamento no intervalo [0,1] de cada variável
    % input: Xe : matriz X para ser preprocessada;
    %        Xet: matriz X de teste;
    % output: X1 : matriz Xe preprocessada;
    %         X1p: matriz Xet preprocessada; 
    X1=(Xe-ones(size(Xe,1),1)*min(Xe))./(ones(size(Xe,1),1)*(max(Xe)-min(Xe)));
    if size(Xe,2)==size(Xet,2)
        X1p=(Xet-ones(size(Xet,1),1)*min(Xe))./(ones(size(Xet,1),1)*(max(Xe)-min(Xe)));
    else
        X1p=Xet;
    end
    
function [X1,X1p] = minmax2(Xe,Xet)
    % Preprocessamento no intervalo [-1, +1] de cada variável
    % input: Xe : matriz X para ser preprocessada;
    %        Xet: matriz X de teste;
    % output: X1 : matriz Xe preprocessada;
    %         X1p: matriz Xet preprocessada; 
    minv=min(Xe);maxv=max(Xe);para1=0.5*(maxv+minv);para2=0.5*(maxv-minv);
    for ki=1:size(Xe,2); X1(:,ki)=(Xe(:,ki)-para1(ki))/para2(ki); end
    if size(Xe,2)==size(Xet,2)
        for ki=1:size(Xet,2); X1p(:,ki)=(Xet(:,ki)-para1(ki))/para2(ki); end
    else
        X1p=Xet;
    end
    
function [X1,X1p] = derivada(Xe,Xet,para1)
    % Derivada com filtro de suavização Savitsky-Golay. 
    % input: Xe : matriz X para ser preprocessada;
    %        Xet: matriz X de teste;
    %      para1: vetor com [janela, grau do polinômio, ordem da derivada]
    % output: X1 : matriz Xe suavizada e diferenciada;
    %         X1p: matriz Xet suavizada e diferenciada;    
    X1=savgol(Xe,para1(1),para1(2),para1(3));   
    if size(Xe,2)==size(Xet,2)
        X1p=savgol(Xet,para1(1),para1(2),para1(3));
    else
        X1p=Xet;
    end
    
function [X1,X1p] = normalizar(Xe,Xet)
    % normalizar
    % input: Xe : matriz X para normalizar;
    %        Xet: matriz X de teste;
    % output: X1 : matriz Xe normalizada;
    %         X1p: matriz Xet normalizada;
    for k2=1:size(Xe,1);   X1(k2,:)=Xe(k2,:)/norm(Xe(k2,:));    end
    
    if size(Xe,2)==size(Xet,2)
        for k2=1:size(Xet,1);   X1p(k2,:)=Xet(k2,:)/norm(Xet(k2,:));    end
    else
        X1p=Xet;
    end
    
    
function [X1,X1p] = emsc(Xe,Xet,opt)
    % Extended Multiplicative Scatter Correction (EMSC) 
    % input: Xe : matriz X para ser preprocessada;
    %        Xet: matriz X de teste;
    %        opt: options
    % output: X1 : matriz Xe suavizada e diferenciada;
    %         X1p: matriz Xet suavizada e diferenciada;    
    para1=mean(Xe);
    if size(Xe,2)==size(Xet,2)
        X1=emscorr(Xe,mean(Xe),opt);
        X1p=emscorr(Xet,mean(Xe),opt);
    else
        X1=emscorr(Xe,mean(Xe),opt);
        X1p=Xet;
    end
    
    
function [X1,X1p] = fir(Xe,Xet,opt)
    % Standardization based on FIR modelling.
    % input: Xe : matriz X para ser preprocessada;
    %        Xet: matriz X de teste;
    %        opt: options
    % output: X1 : matriz Xe suavizada e diferenciada;
    %         X1p: matriz Xet suavizada e diferenciada;    
    if size(Xe,2)==size(Xet,2)
        X1=stdfir(Xe,mean(Xe),opt,1);
        X1p=stdfir(Xet,mean(Xe),opt,1);
    else
        X1=stdfir(Xe,mean(Xe),opt,1);
        X1p=Xet;
    end
    
    
function [X1,X1p] = osc(Xe,Xet,opt1,opt2)
    % OSC
    % input: Xe : matriz X para ser preprocessada;
    %        Xet: matriz X de teste;
    %        opt: options
    % output: X1 : matriz Xe suavizada e diferenciada;
    %         X1p: matriz Xet suavizada e diferenciada;    
    if size(Xe,2)==size(Xet,2)
        [nx,nw,np,nt] = osccalc(Xe,opt1,opt2);
        X1 = oscapp(Xe,nw,np);
        X1p = oscapp(Xet,nw,np);
    else
        [nx,nw,np,nt] = osccalc(Xe,opt1,opt2);
        X1 = oscapp(Xe,nw,np);
        X1p=Xet;
    end