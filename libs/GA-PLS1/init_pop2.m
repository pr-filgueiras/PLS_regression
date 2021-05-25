function variable = init_pop2(nind,Nvar,var_init)
% variable = init_pop2(nind,Nvar,var_init)
% init_pop2: função que regula o número de variáveis selecionadas pelo GA
% nind: número de indivíduos.
% Nvar: número de genes (está relacionado ao número de variáveis espectrais).
% var_init: porcentagem de variáveis (Xcal) a serem utilizadas na população inicial
%
% Exemplo: variable = init_pop2(15,100,20);
%
% Paulo Roberto Filgueiras 30/12/2020
%
variable = zeros(nind,Nvar);
limite = var_init/100;
for ki=1:nind
    for kj=1:Nvar
        if rand<limite; variable(ki,kj)=1; end
    end
end
