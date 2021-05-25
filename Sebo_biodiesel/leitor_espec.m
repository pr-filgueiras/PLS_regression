function [X,y,num,amostras,FileName]=leitor_espec(extensao,replica,planilha,aba,coluna)
% rotina para ler espectros
% input:
% extensao : extesão dos arquivos de leitura;
%  replica : caracteres finais da nomenclatura dos espectros utilizados
%            para identificar as réplicas;
% planilha : nome da planilha Excel onde encontra-se a propriedade de
%            interesse;
%    aba   : nome da aba da planilha Execel a ser lida;
%  coluna  : número da coluna a ser lida.
%
% output:
%     X   : matriz com os espectros;
%     y   : matriz com as propriedades de interesse;
%    num  : eixo-x;
% amostras: nome das amostras; 
% FileName: nome dos arquivos abertos com data de criação.
% 
% exemplo
% [X,y,num,amostras,FileName]=leitor_espec('ASC','1.ASC','biodiesel.xlsx','bio',2);
% [X,y,num,amostras,FileName]=leitor_espec('txt','.txt','amostras pca vitória.xlsx','amostra',2);
%
% Paulo R. Filgueiras  - 06/08/2014
%

% Leitura dos espectros na extensão .ascii
% extensao='ASC';  %alterar
files=dir(strcat('*.',extensao));
for ki=1:numel(files)
    m=numel(files(ki,1).name);
    FileName{ki,1}=char(files(ki,1).name);
    FileName{ki,2}=char(files(ki,1).date);
end

comp_amostra=length(FileName{1,1})-length(replica); % comprimento do nome da amostra

% Leitura dos espectros da planilha excel.
[sample_excel,pfq_excel]=xlsread(planilha,aba);
sample_excel=sample_excel(:,coluna-1);
% identificação das amostras em .ascii
pp2=[]; % Amostras com réplicas
for ki=2:size(pfq_excel,1)
    for kj=1:size(FileName,1)
        if sum(pfq_excel{ki,1}==FileName{kj,1}(1:comp_amostra))==comp_amostra
            pp1=pfq_excel{ki,1};
            pp2=[pp2;pp1];  % nome das amostras
        end
    end
end

% Abrindo os especros
ident=unique(pp2,'rows');
espectro_medio=[];
for ki=1:size(ident,1)
    aa2=[];
    for kj=1:size(FileName,1)
        if sum(ident(ki,:)==FileName{kj,1}(1:comp_amostra))==comp_amostra
            aa1=kj; aa2=[aa2;aa1];
        end
    end
    espectro_amostral=[];
    for kk=1:length(aa2)
        if strcmp(extensao,'ASC'), espectro_i=dlmread(FileName{aa2(kk),1},'\t',56,0);
        elseif strcmp(extensao,'txt'), espectro_i=load(FileName{aa2(kk),1});
        end
        espectro_amostral=[espectro_amostral,espectro_i(:,2)];
    end
    espectro_medio=[espectro_medio,mean(espectro_amostral,2)];
end
espectro_medio=espectro_medio';
num=espectro_i(:,1)';

% Abrindo as propriedades físico-químicas
spectros=[];
propriedades=[];
propriedade_fq=sample_excel;
for ki=1:length(ident)
    for kj=2:size(pfq_excel,1)
        if sum(ident(ki,:)==pfq_excel{kj,1})==comp_amostra && sum(isnan(propriedade_fq(kj-1,:)))==0
            spectros2=ki;  spectros=[spectros;spectros2];
            propriedades2=propriedade_fq(kj-1,:); propriedades=[propriedades;propriedades2];
        end
    end
end
spectros3=espectro_medio(spectros,:); % espectros das amostras selecionadas
ident2=ident(spectros,:); % identificação das amostras

X=spectros3;
y=propriedades;
amostras=ident;








