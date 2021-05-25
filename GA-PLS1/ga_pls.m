function modelo = ga_pls(Xcal,ycal,options)
% modelo = ga_pls(Xcal,ycal,options)
% Função para construção de um modelo GA-PLS




% Paulo Roberto Filgueiras  - 30/12/2020
H = waitbar(0,'Running...','name','GA-PLS variable optimizing... ');
for ki=1:options.Nga

variable = init_pop2(options.nind,options.Nvar,options.var_init); % Variáveis a serem utilizadas na população inicial. 
chrom = variable;

objv=obj_gapls(variable,Xcal,ycal,options.nvl,options.method,options.kfold);
trace = zeros(2,options.maxgen);  
gen = 0;
while gen<options.maxgen      
    fitnv=ranking(objv);                                          % Compute the fitness of each individual
    selch=select('sus',chrom,fitnv,options.ggap);                 % Choose the subset in terms of fitness and gap.
    selch=recombin('xovsp',selch,0.9);                            % Recombination
    selch=mut(selch);                                             % Mutation
    objvsel=obj_gapls(selch,Xcal,ycal,options.nvl,options.method,options.kfold);
    [chrom objv]=reins(chrom,selch,1,1,objv,objvsel);             % Reinsert to keep the population size unchanged
    gen=gen+1;                                                    % Generation plus 1.
    trace(1,gen)=min(objv);                                       % Mininum the optimized fitness
    trace(2,gen)=sum(objv)/length(objv);                          % Compute the average fitness.
end
waitbar(ki/options.Nga);

best_fitness=obj_gapls(chrom,Xcal,ycal,options.nvl,options.method,options.kfold);
[best_fitness,ordem]=sort(best_fitness);
best_var = chrom(ordem,:);
var_sel=find(best_var(1,:));


modelo.RMSECV(ki,1) = best_fitness(1);
modelo.best_var(ki,:) = best_var(1,:);
% modelo.var_sel = var_sel;
end
close(H);

