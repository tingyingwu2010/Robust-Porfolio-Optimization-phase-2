clc;
clear all;


t = readtable('stock_price_list_bse_30.csv');

tab_name=string(t.Properties.VariableNames)';
tab_name=tab_name(2:end);



stock_prices = t{:,2:end};
log_returns = diff(log(stock_prices));
mu_vec = mean(log_returns)';

covariance_mat = cov(stock_prices);


% Markowitz optimization

    ra = 3;
    init_vec = rand(size(log_returns,2),1);
    
    obj_f_mark = @(x) (ra*x'*covariance_mat*x - mu_vec'*x);
    options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',8000);
    options.MaxFunctionEvaluations = 40000;

    [x_mark,fval,exitflag,output] = fmincon(obj_f_mark,init_vec,[],[],ones(1,size(log_returns,2)),1,[],[],[],options);
    
    
 
 % VaR minimization
 
    k = @(e) sqrt((1-e)/e);
    e=0.05;
    obj_f_var = @(x) (k(e)*sqrt(x'*covariance_mat*x) - mu_vec'*x);
    options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',8000);
    options.MaxFunctionEvaluations = 40000;

    [x_var,fval,exitflag,output] = fmincon(obj_f_var,init_vec,[],[],ones(1,size(log_returns,2)),1,[],[],[],options);
    
 format short;
 headings={'BSE_30_Company','Wts_Mark','Wts_VaR'};
 T_output=table(tab_name,x_mark,x_var,'VariableNames',headings);

 writetable(T_output,'output.csv','Delimiter',',');  

 
