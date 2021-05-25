function class_param = calc_class_param(class_calc,class)

% calc_class_param calculates classification parameters:
% error rate, non-error rate, specificity, precision and sensitivity
%
% class_param = calc_class_param(class_calc,class);
%
% input:
% class           class vector [samples x 1]
% class_calc      calculated class vector [1 x n]
% 
% output:
% class_param     structure containing confusion matrix, error rate, non-error rate, 
%                 accuracy, specificity, precision and sensitivity
%
% The main routine is class_gui
%
% Note that a detailed HTML help is provided with the toolbox.
% See the HTML HELP files (help.htm) for futher details and examples
%
% Classification toolbox for MATLAB
% version 4.0 - October 2015
% Davide Ballabioabout
% Milano Chemometrics and QSAR Research Group
% http://michem.disat.unimib.it/chm/

num_class = max([max(class) max(class_calc)]);
nobj = size(class,1);

conf_mat = zeros(num_class,num_class+1);
for g = 1:num_class
    in_class = find(class==g);
    for k = 1:num_class
        conf_mat(g,k) = length(find(class_calc(in_class) == k));
    end
    conf_mat(g,num_class + 1) = length(find(class_calc(in_class) == 0));
end

% sensitivity, specificity, precision, class error, accuracy
accuracy = 0;
for g = 1:num_class
    if sum(conf_mat(:,g)) > 0
        precision(g)   = conf_mat(g,g)/sum(conf_mat(:,g));
        sensitivity(g) = conf_mat(g,g)/sum(conf_mat(g,1:num_class));
    else
        precision(g)   = 0;
        sensitivity(g) = 0;
    end
    in = ones(num_class,1); in(g) = 0;
    red_mat = conf_mat(find(in),1:num_class);
    specificity(g) = 0;
    for k = 1:size(red_mat,2)
        if k ~= g; specificity(g) = specificity(g) + sum(red_mat(:,k)); end;
    end
    if sum(sum(red_mat)) > 0
        specificity(g) = specificity(g)/sum(sum(red_mat));
    else
        specificity(g) = 0;
    end
    false_negative_rate(g) = 1 - sensitivity(g);
    false_positive_rate(g) = 1 - specificity(g);
    accuracy = accuracy + conf_mat(g,g);
end
accuracy = accuracy/sum(sum(conf_mat(:,1:num_class)));

% error rates
not_ass = sum(conf_mat(:,end))/nobj;
ner = mean(sensitivity);
er = 1 - ner;

class_param.conf_mat = conf_mat;
class_param.ner = ner;
class_param.er  = er;
class_param.accuracy  = accuracy;
class_param.not_ass = not_ass;
class_param.precision = precision;
class_param.sensitivity = sensitivity;
class_param.specificity = specificity;