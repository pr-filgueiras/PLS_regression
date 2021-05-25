function [Xcal,Xtest,ycal,ytest,info]=rank_ks(X,y,FileName,ncal,inter,graf)
% Função que aplica o método rank-ks para separar amostras mais
% significativas levando em conta a distribuição de y em intervalos,
% selecionando amostras em cada intervalo pelo algoritmo de Kennard-Stone
% Obs.: Retira as amostras de y mínimo e máximo colocando estes no conjunto
% selecionado e aplica ks nos intervalos definidos para as amostras
% restantes.
% 
% INPUT:
%        X: Matriz com amostras nas linhas e variáveis nas colunas
%        y: Vetor contendo a propriedade de interesse
%        FileName: Célula contendo os nomes das amostras
%        ncal: Porcentagem de amostras a serem selecionadas
%        prop: Nome da propriedade
%        inter: Número de intervalos em que y será dividido
%        graf: 1 - gera gráfico mostrando a separação feita;
%              0 - não gera gráfico
%
% OUTPUT:
%        Xcal: Matriz com dados mais representativos
%        Xtest: Matriz com o restante dos dados
%        ycal: Vetor com dados selecionados da propriedade de interesse
%        ytest: Vetor com o restante de dados da propriedade de interesse
%        info.y_novo: Vetor y em ordem crescente
%        info.X_novo: Matriz X na ordem de y crescente
%        info.FileName_novo: Nome das amostras na nova ordem
%        info.intervalos: Limites dos intervalos estabelecidos
%        info.nomecal: Nome das amostras de calibração
%        info.linhacal: Linhas das amostras de calibração na nova ordem
%        info.nometest: Nome das amostras de teste
%        info.linhatest: Linhas das amostras de teste na nova ordem
%
% Exemplo:
%         [Xcal,Xtest,ycal,ytest,info]=rank_ks(X,y,FileName,70,10,'API',1);
%
% Referência: 
%
% Rayza Rosa T. Rodrigues 02/09/2016
% 


tic;
replicas=rep(y);
if length(FileName(:,1))~=length(y)
    error ('FileName e y devem ter o mesmo tamanho!')
end

if size(X,1)~=length(y)
    error ('X e y devem ter o mesmo número de amostras!')
end

if replicas.tipo==0

% Retirar valores mínimo e máximo
[y2,i]=sort(y);
[~,menor]=min(y2);
[~,maior]=max(y2);
info.y_novo=y2;
y2=info.y_novo;y2([menor,maior])=[];
info.X_novo=X(i,:);
X2=info.X_novo;X2([maior,menor],:)=[];
info.FileName_novo=FileName(i);
FileName2=FileName(i); FileName2([menor,maior])=[];
clear i

%% Definição dos intervalos rank
[jj,j]=hist(y2,inter); % determina quantas amostras em cada intervalo

a=j(2)-j(1);
info.intervalos(inter+1,1)=0; % cria vetor com delimitadores de cada intervalo
info.intervalos(1)=(j(1)-a/2);
for k=2:inter+1
    info.intervalos(k)=info.intervalos(k-1)+a;
end
clear a k j
% determina quais amostras pertencem a cada intervalo (linha 1) e em quais
% linhas estão (linha 2)
am_int(2,inter)={0};
for int=1:inter
    if int==1
        [i,~]=find(y2<info.intervalos(int+1));
    elseif int==inter
        [i,~]=find(y2>=info.intervalos(int));
    else
       [i,~]=find(y2>=info.intervalos(int)&y2<info.intervalos(int+1));
    end
    am_int{1,int}=FileName2(i);
    am_int{2,int}=i;
    clear i
end

% Ranquear por X
am_rank{2,inter}={0};
for i=1:inter
    if length(am_int{2,i})==1
        a=1;
    elseif isempty(am_int{2,i})
        a=[];
    else
    a=ks(X2(am_int{2,i}));         
    end
    am_rank{1,i}=am_int{1,i}(a);
    am_rank{2,i}=am_int{2,i}(a);
end
clear i a

% Selecionar as ncal% primeiras de am_rank
selec{2,inter}={0};
for i=1:inter
    a=round(ncal/100*length(am_rank{1,i}));
    selec{1,i}=am_rank{1,i}(1:a);
    selec{2,i}=am_rank{2,i}(1:a);
    not_selec{1,i}=am_rank{1,i}(a+1:end);
    not_selec{2,i}=am_rank{2,i}(a+1:end);
end

nomecal=[];
linhacal=[];
for i=1:inter
    a=selec{1,i};
    nomecal=[nomecal;a];
    b=selec{2,i};
    linhacal=[linhacal;b];
end
info.nomecal=nomecal;
info.linhacal=linhacal;
clear i a b

nometest=[];
linhatest=[];
for i=1:inter
    a=not_selec{1,i};
    nometest=[nometest;a];
    b=not_selec{2,i};
    linhatest=[linhatest;b];  
end
info.nometest=nometest;
info.linhatest=linhatest+1;

Xcal=X2(linhacal,:);
ycal=y2(linhacal);
Xtest=X2(linhatest,:);
ytest=y2(linhatest);

ycal=[info.y_novo(menor);ycal;info.y_novo(maior)];
Xcal=[info.X_novo(menor,:);Xcal;info.X_novo(maior,:)];
info.nomecal=[info.FileName_novo(menor,:);info.nomecal;info.FileName_novo(maior,:)];
info.linhacal=[menor;info.linhacal+1;maior];

if graf==1
    plot(ycal,ycal,'*b'); hold on; plot (ytest,ytest,'*r')
    legend ('ycal','ytest'), grid on
    xlabel('Valores de referência')
    ylabel('Valores de referência')
    % xlabel ({prop ,'(referência)'})
    % ylabel ({prop , '(previsto)'})
   % vline ([info.intervalos])
    figure (2)
    [a,b]=hist(y,inter);
    h = bar(b,a);
    set(h,'facecolor','b')
    hold on
    [a2,b2]=hist(ycal,b);
    h2 = bar(b2,a2);
    set(h,'facecolor','r')
    legend ('y','ycal')
    % xlabel (prop)
    xlabel('Propriedade (y)')
    ylabel ('Número de amostras')
end
     
else
    fprintf('Este algoritmo não está pronto para trabalhar com réplicas.')
end
toc;

%% Subfunção para detectar réplicas
function [replicas]=rep(y)
if y(1,1)~=y(2,1)
    replicas.tipo=0;
    return
elseif y(1,1)==y(2,1)&& y(1,1)~=y(3,1)
    replicas.tipo=2;
elseif y(1,1)==y(2,1)&&y(1,1)==y(3,1)&&y(1,1)~=y(4,1)
    replicas.tipo=3;
elseif y(1,1)==y(2,1)&&y(1,1)==y(3,1)&&y(1,1)==y(4,1)&&y(1,1)~=y(5,1)
    replicas.tipo=4;
else
    error ('Este algoritmo não comporta mais de 4 réplicas')
end    

if rem(length(y)/replicas.tipo,1) ~= 0
    error ('Algumas amostras não contém o mesmo número de réplicas')
end
    
y2.sample=[];y2.sample=[y2.sample;1];
y2.numero=[];
 for ki=2:size(y,1)
     if sum(y(ki,:)==y(ki-1,:))==size(y,2)
        y2.numero=ki;
     else
        y2.sample=[y2.sample;ki];
    end
end
    aa2=y2.sample(2:end)-1;aa2=[aa2;y2.numero];
    replicas.sample_position=[y2.sample,aa2];  % vetor com as posições iniciais e finais para as réplicas.
    clear aa2 ki y2.sample