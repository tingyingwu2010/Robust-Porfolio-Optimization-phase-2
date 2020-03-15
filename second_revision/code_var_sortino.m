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
% cnt = 0;
% for i = 1:size(stock_prices,2)
%     ks_vec=stock_prices(:,i);
%     m_ks=mean(ks_vec);
%     sig_ks=std(ks_vec);
%     ks_vec=ks_vec-m_ks;
%     ks_vec=ks_vec/sig_ks;
%     [h,p] = jbtest(ks_vec);
%     if (h==1)
%         cnt = cnt+1;
%     end
%     
% end
% disp(cnt)
% return;

% Uncomment if needed to use the simulated data

rng default  % For reproducibility
% m=1000;   % If #simulations is 1000
m=size(stock_prices,1); % If #simulations is same as market data
temp_data = mvnrnd(mu,covariance,m);
stock_prices=temp_data;
mu = mean(stock_prices);
mu = mu';
covariance = cov(stock_prices);

n_stocks = size(stock_prices,2);
Y = min(stock_prices - mu',0);
semi_covar = zeros(n_stocks, n_stocks);
covar = zeros(n_stocks, n_stocks);
for i = 1:n_stocks
    temp = Y(:,i);
    temp_c = stock_prices(:,i) - mu';
    temp_mat = temp.*Y;
    temp_mat_c = temp_c.*(stock_prices-mu');
    tot_sum = sum(temp_mat)./(size(stock_prices,1)-1);
    tot_sum_c = sum(temp_mat_c)./(size(stock_prices,1)-1);
    semi_covar(i,:) = tot_sum;
    covar(i,:) = tot_sum_c;
end

k = @(e) sqrt((1-e)/e);
%k = @(e) -1*norminv(e);
e_range  = 0.0001:5*10^-3:0.1;
% e_range = 0.0001:10^(-4):0.01

mean_vals_base = [];
sd_vals_base = [];

init = rand(size(stock_prices,2),1);
init = init./sum(init);
for i=1:size(e_range,2)
    
    e = e_range(1,i)
%     init = rand(size(stock_prices,2),1);
%     init = init./sum(init);
%     delta = sqrt(chi2inv(1-alpha,size(stock_prices,2)));
%     cov_err=covariance./size(stock_prices,1);
    maxim = @(x) (k(e)*sqrt(x'*covariance*x) - mu'*x);
    options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',4000);
    options.MaxFunctionEvaluations = 20000;
    A=eye(size(stock_prices,2));
    b=zeros(size(stock_prices,2),1);
    A=A.*(-1);
    [x,fval,exitflag,output] = fmincon(maxim,init,A,b,ones(1,size(stock_prices,2)),1,[],[],[],options);
    
    mean_vals_base = [mean_vals_base, mu'*x];
    sd_vals_base = [sd_vals_base, sqrt(x'*semi_covar*x)];
%     sqrt(x'*semi_covar*x)
%     sqrt(x'*covariance*x)
%     sqrt(x'*covar*x)
    
end
%% Worst Case VaR

mean_vals_wvar = [];
sd_vals_wvar = [];
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


for i=1:size(e_range,2)
    
    e = e_range(1,i)
%     init = rand(size(stock_prices,2),1);
    maxim = @(x) (k(e)*sqrt(x'*upper_cov*x) - lower_mu'*x);
    options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',4000);
    options.MaxFunctionEvaluations = 20000;
    A=eye(size(stock_prices,2));
    b=zeros(size(stock_prices,2),1);
    A=A.*(-1);
    [x,fval,exitflag,output] = fmincon(maxim,init,A,b,ones(1,size(stock_prices,2)),1,[],[],[],options);
    
    mean_vals_wvar = [mean_vals_wvar, mu'*x];
    sd_vals_wvar = [sd_vals_wvar, sqrt(x'*semi_covar*x)];
    
end


%% 



risk_free = log(1.06)/365;
base = (mean_vals_base - risk_free)./sd_vals_base;
wvar = (mean_vals_wvar - risk_free)./sd_vals_wvar;


F=figure(2); hold on;
box on
grid on
plot(e_range, base,'-o');
plot(e_range, wvar,'-s');
lgd = legend('VaR','WVaR');
lgd.Location = 'southeast';
ylabel('Sortino Ratio');
xlabel('\epsilon');

% change the names of the files and folders accordingly.
% saveas(F,'./JPEGs/bse30_simulated/sr_exact_cheb.jpeg');
saveas(F,'./EPSs/bse30_simulated/sr_exact_cheb.eps','epsc');
hold off

format short;

id_ra_range=0.0001:20*10^-3:0.1;
Tab=zeros(length(id_ra_range),8);
[tf,loc]=ismember(e_range,id_ra_range);
idx=[1:length(e_range)];
idx=idx(tf);
idx=idx(loc(tf));
Tab(1:end,1)=e_range(idx)';
Tab(1:end,2)=ones(length(id_ra_range),1).*risk_free;
Tab(1:end,3)=mean_vals_base(idx)';
Tab(1:end,4)=sd_vals_base(idx)';
Tab(1:end,5)=mean_vals_wvar(idx)';
Tab(1:end,6)=sd_vals_wvar(idx)';
Tab(1:end,7)=base(idx)';
Tab(1:end,8)=wvar(idx)';



headings={'ConfidenceLevel','RiskfreeRate','BaseVar_u','BaseVar_sig','WorstCaseVar_u','WorstCaseVar_sig','BaseVar_SR','WorstCaseVar_SR'};

Avg=zeros(1,3);
Avg(1,1)=mean(base);
Avg(1,2)=mean(wvar);

Avg(1,3)=Avg(1,2)-Avg(1,1);



% change the names of the files and folders accordingly.
tab_loc='./tables/bse30_simulated/tab_exact_cheb.csv';
headings=strjoin(headings, ',');
fid_tab=fopen(tab_loc,'w'); 
fprintf(fid_tab,'%s\n',headings);
fclose(fid_tab);
dlmwrite(tab_loc,Tab,'-append','delimiter', ',', 'precision', 3);



avg_headings={'BaseVar_SR','WorstCaseVar_SR','Diff_SR'};
% change the names of the files and folders accordingly.
avg_loc='./tables/bse30_simulated/avg_exact_cheb.csv';
avg_headings = strjoin(avg_headings, ',');
fid_avg = fopen(avg_loc,'w'); 
fprintf(fid_avg,'%s\n',avg_headings);
fclose(fid_avg);
dlmwrite(avg_loc,Avg,'-append','delimiter', ',', 'precision', 3);



% mark_size = 5;
% F=figure(1); hold on;
% box on
% grid on
% plot(sd_vals_base, mean_vals_base,'-s','markers',mark_size);
% plot(sd_vals_wvar, mean_vals_wvar,'-o','markers',mark_size);
% lgd = legend('Base VaR','Worst case VaR');
% lgd.Location = 'southeast';
% ylabel('Return');
% xlabel('Standard Deviation')
% 
% saveas(F,'./JPEGs/bse30_simulated/ef_exact_cheb.jpeg');
% saveas(F,'./EPSs/bse30_simulated/ef_exact_cheb.eps','epsc');
% 
% hold off


