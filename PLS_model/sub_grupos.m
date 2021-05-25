function SubGrupos = sub_grupos(X,N,MTR)
% Rotina para criar subgrupos de amostras. M�trica utilizada: dist�ncia
% euclidiana padronizada.
% input:
% X : matriz de espectros
% N : n�mero de subgrupos
% MTR: m�trica utilizada para determinar as dist�ncias: 
%   'seuclidean'; 'euclidean'
%
% output:
% SubGrupos : amostras dos subgrupos
% .dist : dist�ncia de cada amostra ao centro do conjunto de dados;
% .D2 : dist�ncias ordenadas crescentes;
% .D3 : amostras ordenadas crescentes;
% .Nint: n�mero inicial e final da ordem de cada amostra que formar� o
%       subgrupo;
% .Samplesi: amostras que ir�o compor os subgrupos.
%
% exemplo:
% SubGrupos = sub_grupos(X,10,'euclidean')
% SubModels=sub_pls_model(X,y,SubGrupos,4,'center');
% Subpred=sub_pls_pred(SubModels,Xt,yt);
%
% Paulo R. Filgueiras - 29/01/2016
%
if nargin==2;
    MTR = 'euclidean';
end
%----- determinar as as dist�ncias 
D2 = pdist([mean(X);X],MTR);
D3 = squareform(D2);
D.dist=D3(2:end,1);

[D.D2,D.D3]=sort(D.dist,'descend');

semente=D.D3(1);

D2 = pdist([X(semente,:);X],MTR);
D3 = squareform(D2);
D.dist=D3(2:end,1);

D.Nint=round(linspace(1,length(D.dist),N+1));
if ~(D.Nint(end)==length(D.dist)), D.Nint(end+1)=length(D.dist); end

for ki=1:length(D.Nint)-1
    a1=D.Nint(ki):D.Nint(ki+1);
    D.Samplesi{ki,1}=D.D3(a1);
end

SubGrupos=D;
SubGrupos.MTR=MTR;


