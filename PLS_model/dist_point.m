function vetor_dist = dist_point(X,xref,metrica)
% vetor_dist = dist_point(X,xref,metrica)
% Rotina para calcular a distância entre pontos no espaço.
% Esta rotina calcula a distância das amostras da matriz X à amostras xref.
%
% input:
%  X(n,m)   : matriz de dados, contendo n amostras com m variáveis. 
%   xref    : ponto para cálculo da distância. default mean(X).
% 'metrica' : métrica utilizada para calcular a distância entre pontos.
%    ++ 'euclidiana' : utilização da métrica euclidiana;
%    ++ 'standeuclid' : utilização da métrica euclidiana padronizada;
%    ++ 'mahalanobis' : utilização da métrica mahalanobis; default.
%
% output;
% vetor_dist: vetor contendo a distância de cada amostra da matrix X à
%             amostra xref.
% 
% Exemplo:
% vetor_dist = dist_point(X,xref,'mahalanobis');
%
%  Paulo R. Filgueiras   - 09/09/2014
%

if nargin==1
    xref=mean(X);
    metrica='mahalanobis';
elseif nargin==2
    metrica='mahalanobis';
end

[nx,mx] = size(X);

if strcmp(metrica,'euclidiana')
    vetor_dist=[];
    for ki=1:nx
        distancia=sum((X(ki,:)-xref)*(X(ki,:)-xref)'); % distância euclidiana quadrática
        vetor_dist=[vetor_dist;sqrt(distancia)]; % distância euclidiana
    end
elseif strcmp(metrica,'standeuclid')
    Xmean = X-repmat(mean(X),nx,1); % Matriz X centrada na média.
    covariance_X = Xmean'*Xmean/nx; % matrix covariance populacional.
    variance_X = diag(diag(covariance_X));   % matriz diagonal de variâncias
    vetor_dist=[];
    for ki=1:nx
        distancia=sum((X(ki,:)-xref)*inv(variance_X)*(X(ki,:)-xref)'); % distância euclidiana quadrática
        vetor_dist=[vetor_dist;sqrt(distancia)]; % distância euclidiana
    end
elseif strcmp(metrica,'mahalanobis')
    [scores,loads,varexp] = pca(X,nx-1);
    Ncomp = component_number(varexp);
    clear Xmean
    scores=scores(:,1:Ncomp);           % Matriz de scores dos dados X.
    scores_ref = xref*loads(:,1:Ncomp); % Scores do vetor de referência. 
    Xmean = scores-repmat(mean(scores),nx,1); % Matriz X centrada na média.
    covariance_scores = Xmean'*Xmean/nx; % matrix covariance populacional.
    vetor_dist=[];
    for ki=1:nx
        distancia=sum((scores(ki,:)-scores_ref)*inv(covariance_scores)*(scores(ki,:)-scores_ref)'); % distância euclidiana quadrática
        vetor_dist=[vetor_dist;sqrt(distancia)]; % distância euclidiana
    end
end

% hb = colorbar;   ylabel(hb,'Mahalanobis Distance')
% http://people.revoledu.com/kardi/tutorial/index.html


% Sub-rotinas
function [scores,loads,varexp] = pca(X,npc)
%  Análise por Componentes Principais por SVD
[n,m] = size(X);
if m < n
  cov = (X'*X);
  [u,s,v] = svd(cov);
  PCnumber = (1:m)';
else
  cov = (X*X');
  [u,s,v] = svd(cov);
  v = X'*v;
  for i = 1:n
    v(:,i) = v(:,i)/norm(v(:,i));
  end
  PCnumber = (1:n)';
end
individualExpVar = diag(s)*100/(sum(diag(s)));
varexp  = [PCnumber individualExpVar cumsum(individualExpVar)];
loads  = v(:,1:npc);
scores = X*loads;

function Ncomp = component_number(varexp)
% sub-rotina para determinar o número de componentes (PC) na redução da
% dimensão espacial (PCA) realizada na determinada da distância de
% mahalanobis.
% Condição: #Ncomp = # PC com variância explicada maior que 0.5%
PCs=find(varexp(:,2)>=0.5);
if length(PCs)==1
    Ncomp=2;
else
    Ncomp=max(PCs);
end

