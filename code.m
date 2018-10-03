
clc;
clear all;

table = readtable('final_list.csv');

stock_prices = table{:,2:end};

stock_prices = diff(log(stock_prices));

mu = mean(stock_prices);
mu = mu';
covariance = cov(stock_prices);

lambda = 1/100;
maxim = @(x) (lambda*x'*covariance*x - mu'*x);


init = rand(31,1);
init = init./sum(init);
options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',4000);
options.MaxFunctionEvaluations = 10000;
[x,fval,exitflag,output] = fmincon(maxim,init,[],[],ones(1,31),1,[],[],[],options);

risk_aversion = 1000:100:20000;
mean_vals_mark = [];
sd_vals_mark = [];
for i=1:size(risk_aversion,2)
    
    lambda = risk_aversion(1,i)
    init = rand(31,1);
    init = init./sum(init);
    maxim = @(x) (lambda*x'*covariance*x - mu'*x);
    options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',4000);
    options.MaxFunctionEvaluations = 10000;
    [x,fval,exitflag,output] = fmincon(maxim,init,[],[],ones(1,31),1,[],[],[],options);
    
    
    mean_vals_mark = [mean_vals_mark, mu'*x];
    sd_vals_mark = [sd_vals_mark, sqrt(x'*covariance*x)];
end

% figure(1)
% plot(sd_vals, mean_vals);
% figure(2)
% plot(100:1:15000,mean_vals);
% figure(3)
% plot(100:1:15000, sd_vals);

% Box uncertainty

mean_vals_box = [];
sd_vals_box = [];
alpha = 0.05;
for i=1:size(risk_aversion,2)
    
    lambda = risk_aversion(1,i)
    init = rand(31,1);
    init = init./sum(init);
    delta = abs(norminv(alpha/2)).*sqrt(diag(covariance))./size(stock_prices,1);
    maxim = @(x) (lambda*x'*covariance*x - mu'*x + delta'*abs(x));
    options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',4000);
    options.MaxFunctionEvaluations = 10000;
    [x,fval,exitflag,output] = fmincon(maxim,init,[],[],ones(1,31),1,[],[],[],options);
    
    
    mean_vals_box = [mean_vals_box, mu'*x];
    sd_vals_box = [sd_vals_box, sqrt(x'*covariance*x)];
    
end


% figure(1)
% plot(sd_vals, mean_vals);
% figure(2)
% plot(1000:100:15000,mean_vals);
% figure(3)
% plot(1000:100:15000, sd_vals);

% Ellipsoidal uncertainty

mean_vals_ellipsoid = [];
sd_vals_ellipsoid = [];

alpha = 0.05;
for i=1:size(risk_aversion,2)
    
    lambda = risk_aversion(1,i)
    init = rand(31,1);
    init = init./sum(init);
    delta = chi2inv(alpha,size(stock_prices,2));
    maxim = @(x) (lambda*x'*covariance*x - mu'*x + delta*sqrt(x'*covariance*x));
    options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',4000);
    options.MaxFunctionEvaluations = 10000;
    [x,fval,exitflag,output] = fmincon(maxim,init,[],[],ones(1,31),1,[],[],[],options);
    
    
    mean_vals_ellipsoid = [mean_vals_ellipsoid, mu'*x];
    sd_vals_ellipsoid = [sd_vals_ellipsoid, sqrt(x'*covariance*x)];
    
end


figure(1), hold on
plot(sd_vals_mark, mean_vals_mark);
plot(sd_vals_box, mean_vals_box);
plot(sd_vals_ellipsoid, mean_vals_ellipsoid);
legend('Markowitz','With box uncertainty','With ellipsoid uncertainty');
ylabel('return');
xlabel('Standard Deviation')
hold off

risk_free = log(1.06)/365;
mark = (mean_vals_mark - risk_free)./sd_vals_mark;
box = (mean_vals_box - risk_free)./sd_vals_box;
ellipsoid = (mean_vals_ellipsoid - risk_free)./sd_vals_ellipsoid;




% figure(2)
% plot(1000:100:15000,mean_vals);
% figure(3)
% plot(1000:100:15000, sd_vals);

% lower_mu = zeros(size(stock_prices,2),1);
% upper_mu = zeros(size(stock_prices,2),1);
% lower_cov = zeros(size(stock_prices,2),size(stock_prices,2));
% upper_cov = zeros(size(stock_prices,2),size(stock_prices,2));
% 
% iterations = 8000;
% mu_matrix = zeros(size(stock_prices,2),iterations);
% cov_matrix = zeros(size(stock_prices,2),size(stock_prices,2),iterations);
% rows = 1:size(stock_prices,1);
% 
% for i=1:iterations
%     
%    temp_rows = datasample(rows,size(stock_prices,1));
%    temp_data = stock_prices(temp_rows,:);
%    mu_matrix(:,i) = mean(temp_data)';
%    cov_matrix(:,:,i) = cov(temp_data);

    temp_data = mvnrnd(mu,covariance,size(stock_prices,1));
    mu_matrix(:,i) = mean(temp_data)';
    cov_matrix(:,:,i) = cov(temp_data);
    
% end
% 
% for i = 1:size(stock_prices,2)
% 
%     temp = sort(mu_matrix(i,:));
%     lower_mu(i,1) = temp(iterations*0.025);
%     upper_mu(i,1) = temp(iterations*0.975);
%     
% end
% 
% 
% for i = 1:size(stock_prices,2)
% 
%     for j=1:size(stock_prices,2)
%         
%        temp = sort(cov_matrix(i,j,:));
%        lower_cov(i,j) = temp(iterations*0.025);
%        upper_cov(i,j) = temp(iterations*0.975);
%       
%     end
%     
% end

% Joint UNcertainty
% alpha = 0.05;
% for i=10000:100:15000
%     
%     lambda = i
%     init = rand(31,1);
%     maxim = @(x) (lambda*x'*upper_cov*x - lower_mu'*x);
%     options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',4000);
%     options.MaxFunctionEvaluations = 10000;
%     [x,fval,exitflag,output] = fmincon(maxim,init,[],[],ones(1,31),1,[],[],[],options);
%     
%     
%     mean_vals = [mean_vals, mu'*x];
%     sd_vals = [sd_vals, sqrt(x'*covariance*x)];
%     
% end
% 
% 
% figure(1)
% plot(sd_vals, mean_vals);
% figure(2)
% plot(10000:100:15000,mean_vals);
% figure(3)
% plot(10000:100:15000, sd_vals);
% 

