function model = cars_varsel(X,y,options)
% Seleção de variáveis pelo método CARS.
% Competitive Adaptative Reweighted Sampling Method.
%
% Exemplo
% model = cars_varsel(X,y,options)
%

tic;

if nargin ==2
options.Xpretreat={'msc'}; % pretreatment method.
options.vene = 5;   %k-fold.
options.vl = 10;    % Maximal Principle to extract
options.mcsr = 100; % Number of Monte Carlo Sampling runs.
options.runs = 50 ; % o método CARS será execultados n=50 vezes.
options.graficos = 1; % graphic output   
end

X2=X;
X=pretrat(X,[],options.Xpretreat);

matlabpool
parfor ki=1:options.runs
    F=carspls(X,y,options.vl,options.vene,'center',options.mcsr,1,1,1);
    F.W=[];
    RMSECVmin(ki,1)=F.RMSECV_min;
    CARS2{ki,1}=F;
end
matlabpool close

varsel=[]; % variáveis selecionadas com réplicas
for ki=1:numel(CARS2)
    varsel=[varsel,CARS2{ki,1}.vsel];
end
varsel=sort(varsel); 

varsel2=unique(varsel); % somente as variáveis selecionadas
for ki=1:size(X,2)
    if length(find(varsel2==ki))==1
        varsel3(ki) = length(find(varsel==ki));
    else
        varsel3(ki)=0; % variáveis selecionadas com frequências de seleção
    end
end


% output model
model.data = date;
model.time = toc;
model.CARS = CARS2; 
model.RMSECV = RMSECVmin; % valores mínimos de RMSECV.
model.varsel = varsel2; % mesmo comprimento da matrix X.
model.varsel2 = varsel3; % apenas as variáveis selecionadas.


if options.graficos==1
    figure(1)  
    subplot(2,1,1), stem(1:size(X,2),X2','Marker', 'none')
    ylabel('Absorbance','FontSize',14)
    subplot(2,1,2), bar(model.varsel2)
    xlabel('Variable','FontSize',14)
    ylabel('Frequency','FontSize',14)
    set(gcf,'Color','white')
end

