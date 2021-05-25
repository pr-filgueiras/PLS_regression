function variable = init_pop2(nind,Nvar,var_init)
% variable = init_pop2(nind,Nvar,var_init)
% init_pop2: fun��o que regula o n�mero de vari�veis selecionadas pelo GA
% nind: n�mero de indiv�duos.
% Nvar: n�mero de genes (est� relacionado ao n�mero de vari�veis espectrais).
% var_init: porcentagem de vari�veis (Xcal) a serem utilizadas na popula��o inicial
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
