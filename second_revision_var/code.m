
clc;
clear all;

% Read the csv file (Change it to "final_list100.csv" for BSE100)
table = readtable('./data_related/final_list30.csv');
% table = readtable('./data_related/final_list100.csv');

stock_prices = table{:,2:end};
stock_prices = diff(log(stock_prices));
mu = mean(stock_prices);
mu = mu';
covariance = cov(stock_prices);


% Uncomment if needed to use the simulated data

% rng default  % For reproducibility
% % m=1000;   % If #simulations is 1000
% m=size(stock_prices,1); % If #simulations is same as market data
% temp_data = mvnrnd(mu,covariance,m);
% stock_prices=temp_data;
% mu = mean(stock_prices);
% mu = mu';
% covariance = cov(stock_prices);


risk_aversion=2:0.25:4;

% Implementing Markowitz model

mean_vals_mark = [];
sd_vals_mark = [];

for i=1:size(risk_aversion,2)
    
    lambda = risk_aversion(1,i)
    init = rand(size(stock_prices,2),1);
    init = init./sum(init);
    maxim = @(x) (lambda*x'*covariance*x - mu'*x);
    options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',4000);
    options.MaxFunctionEvaluations = 20000;
    A=eye(size(stock_prices,2));
    b=zeros(size(stock_prices,2),1);
    A=A.*(-1);
    [x,fval,exitflag,output] = fmincon(maxim,init,A,b,ones(1,size(stock_prices,2)),1,[],[],[],options);
    mean_vals_mark = [mean_vals_mark, mu'*x];
    sd_vals_mark = [sd_vals_mark, sqrt(x'*covariance*x)];
end



% Implementing Box model.

mean_vals_box = [];
sd_vals_box = [];
alpha = 0.05;

for i=1:size(risk_aversion,2)
    
    lambda = risk_aversion(1,i)
    init = rand(size(stock_prices,2),1);
    init = init./sum(init);
    delta = abs(norminv(alpha/2)).*sqrt(diag(covariance))./size(stock_prices,1);
    maxim = @(x) (lambda*x'*covariance*x - mu'*x + delta'*abs(x));
    options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',4000);
    options.MaxFunctionEvaluations = 20000;
    A=eye(size(stock_prices,2));
    b=zeros(size(stock_prices,2),1);
    A=A.*(-1);
    [x,fval,exitflag,output] = fmincon(maxim,init,A,b,ones(1,size(stock_prices,2)),1,[],[],[],options);
    mean_vals_box = [mean_vals_box, mu'*x];
    sd_vals_box = [sd_vals_box, sqrt(x'*covariance*x)];
    
end


% Implementing Ellipsoidal model.

mean_vals_ellipsoid = [];
sd_vals_ellipsoid = [];

for i=1:size(risk_aversion,2)
    
    lambda = risk_aversion(1,i)
    init = rand(size(stock_prices,2),1);
    init = init./sum(init);
    delta = sqrt(chi2inv(1-alpha,size(stock_prices,2)));
    cov_err=covariance./size(stock_prices,1);
    maxim = @(x) (lambda*x'*covariance*x - mu'*x + delta*sqrt(x'*cov_err*x));
    options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',4000);
    options.MaxFunctionEvaluations = 20000;
    A=eye(size(stock_prices,2));
    b=zeros(size(stock_prices,2),1);
    A=A.*(-1);
    [x,fval,exitflag,output] = fmincon(maxim,init,A,b,ones(1,size(stock_prices,2)),1,[],[],[],options);
    
    mean_vals_ellipsoid = [mean_vals_ellipsoid, mu'*x];
    sd_vals_ellipsoid = [sd_vals_ellipsoid, sqrt(x'*covariance*x)];
    
end

% Implementing Separable model

% Calculating Bounds by using Bootstrap algorithm.

lower_mu = zeros(size(stock_prices,2),1);
upper_mu = zeros(size(stock_prices,2),1);
lower_cov = zeros(size(stock_prices,2),size(stock_prices,2));
upper_cov = zeros(size(stock_prices,2),size(stock_prices,2));

iterations = 8000;
mu_matrix = zeros(size(stock_prices,2),iterations);
cov_matrix = zeros(size(stock_prices,2),size(stock_prices,2),iterations);
rows = 1:size(stock_prices,1);

for i=1:iterations
    
%     For Non parametric bootstrap
   temp_rows = datasample(rows,size(stock_prices,1));
   temp_data = stock_prices(temp_rows,:);
   mu_matrix(:,i) = mean(temp_data)';
   cov_matrix(:,:,i) = cov(temp_data);

%    For parametric bootstrap
%     temp_data = mvnrnd(mu,covariance,size(stock_prices,1));
%     mu_matrix(:,i) = mean(temp_data)';
%     cov_matrix(:,:,i) = cov(temp_data);
    
end

for i = 1:size(stock_prices,2)

    temp = sort(mu_matrix(i,:));
    lower_mu(i,1) = temp(iterations*0.025);
    upper_mu(i,1) = temp(iterations*0.975);
    
end


for i = 1:size(stock_prices,2)

    for j=1:size(stock_prices,2)
        
       temp = sort(cov_matrix(i,j,:));
       lower_cov(i,j) = temp(iterations*0.025);
       upper_cov(i,j) = temp(iterations*0.975);
      
    end
    
end

mean_vals_var = [];
sd_vals_var = [];

for i=1:size(risk_aversion,2)
    
    lambda = risk_aversion(1,i)
    init = rand(size(stock_prices,2),1);
    maxim = @(x) (lambda*x'*upper_cov*x - lower_mu'*x);
    options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',4000);
    options.MaxFunctionEvaluations = 20000;
    A=eye(size(stock_prices,2));
    b=zeros(size(stock_prices,2),1);
    A=A.*(-1);
    [x,fval,exitflag,output] = fmincon(maxim,init,A,b,ones(1,size(stock_prices,2)),1,[],[],[],options);
    
    mean_vals_var = [mean_vals_var, mu'*x];
    sd_vals_var = [sd_vals_var, sqrt(x'*covariance*x)];
    
end

% Plotting the efficient frontier

mark_size = 5;
F=figure(1); hold on;
box on
grid on
plot(sd_vals_mark, mean_vals_mark,'-o','markers',mark_size);
plot(sd_vals_box, mean_vals_box,'-s','markers',mark_size);
plot(sd_vals_ellipsoid, mean_vals_ellipsoid,'k-.','markers',mark_size);
plot(sd_vals_var, mean_vals_var,'-*','markers',mark_size);
lgd = legend('Vanilla Markowitz','With Box uncertainty','With Ellipsoid uncertainty', 'With Separable uncertainty');
lgd.Location = 'southeast';
ylabel('Return');
xlabel('Standard Deviation')
hold off

% change the name of the files and folder accordingly.
saveas(F,'./JPEG/bse30_market/ef_ideal_range.jpeg');
saveas(F,'./EPSWTs/bse30_simulated/ef_ideal_range.eps','epsc');

risk_free = log(1.06)/365;
mark = (mean_vals_mark - risk_free)./sd_vals_mark;
box_set = (mean_vals_box - risk_free)./sd_vals_box;
ellipsoid = (mean_vals_ellipsoid - risk_free)./sd_vals_ellipsoid;
joint_var = (mean_vals_var - risk_free)./sd_vals_var;

% Comment the snippets to construct the tables.

% format short;
% Tab=zeros(5,14);
% id_ra_range=2:0.5:4;
% [tf,loc]=ismember(risk_aversion,id_ra_range);
% idx=[1:length(risk_aversion)];
% idx=idx(tf);
% idx=idx(loc(tf));
% Tab(1:end,1)=risk_aversion(idx)';
% Tab(1:end,2)=ones(length(id_ra_range),1).*risk_free;
% Tab(1:end,3)=mean_vals_mark(idx)';
% Tab(1:end,4)=sd_vals_mark(idx)';
% Tab(1:end,5)=mean_vals_box(idx)';
% Tab(1:end,6)=sd_vals_box(idx)';
% Tab(1:end,7)=mean_vals_ellipsoid(idx)';
% Tab(1:end,8)=sd_vals_ellipsoid(idx)';
% Tab(1:end,9)=mean_vals_var(idx)';
% Tab(1:end,10)=sd_vals_var(idx)';
% Tab(1:end,11)=mark(idx)';
% Tab(1:end,12)=box_set(idx)';
% Tab(1:end,13)=ellipsoid(idx)';
% Tab(1:end,14)=joint_var(idx)';
% 
% 
% 
% headings={'Risk_aversion','Riskfree','Mark_u','Mark_sig','Box_u','Box_sig','Ellip_u','Ellip_sig','Sep_u','Sep_sig','Mark_SR','Box_SR','Ellip_SR','Sep_SR'};
% 
% Avg=zeros(1,6);
% Avg(1,1)=mean(Tab(1:end,11));
% Avg(1,2)=mean(Tab(1:end,12));
% Avg(1,3)=mean(Tab(1:end,13));
% Avg(1,4)=mean(Tab(1:end,14));
% Avg(1,5)=Avg(1,3)-Avg(1,1);
% Avg(1,6)=Avg(1,4)-Avg(1,1);
% 
% 
% % change the names of the files and folders accordingly.
% tab_loc='./tables/bse30_market/tab_ideal_range.csv';
% headings=strjoin(headings, ',');
% fid_tab=fopen(tab_loc,'w'); 
% fprintf(fid_tab,'%s\n',headings);
% fclose(fid_tab);
% dlmwrite(tab_loc,Tab,'-append','delimiter', ',', 'precision', 3);
% 
% 
% 
% avg_headings={'Mark','Box','Ellip','Sep','Diff_Ellip_Mark','Diff_Sep_Mark'};
% % change the names of the files and folders accordingly.
% avg_loc='./tables/bse30_market/avg_ideal_range.csv';
% avg_headings = strjoin(avg_headings, ',');
% fid_avg = fopen(avg_loc,'w'); 
% fprintf(fid_avg,'%s\n',avg_headings);
% fclose(fid_avg);
% dlmwrite(avg_loc,Avg,'-append','delimiter', ',', 'precision', 3);


% Plotting the Sharpe Ratios

F=figure(2); hold on;
box on
grid on
plot(risk_aversion, mark,'-o');
plot(risk_aversion, box_set,'-s');
plot(risk_aversion, ellipsoid,'k-.');
plot(risk_aversion, joint_var,'-*');
lgd = legend('Vanilla Markowitz','With Box uncertainty','With Ellipsoid uncertainty', 'With Separable uncertainty');
lgd.Location = 'southeast';
ylabel('Sharpe Ratio');
xlabel('Risk Aversion');

% change the names of the files and folders accordingly.
saveas(F,'./JPEG/bse30_market/sr_ideal_range.jpeg');
saveas(F,'./EPSWTs/bse30_market/sr_ideal_range.eps','epsc');
hold off



