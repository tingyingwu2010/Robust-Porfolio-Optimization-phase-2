clc;
clear all;

% Read the csv file (Change it to "final_list100.csv" for BSE100)
table = readtable('./data_related/final_list100.csv');
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

% rng default  % For reproducibility
% %m=1000;   % If #simulations is 1000
% m=size(stock_prices,1); % If #simulations is same as market data
% temp_data = mvnrnd(mu,covariance,m);
% stock_prices=temp_data;
% mu = mean(stock_prices);
% mu = mu';
% covariance = cov(stock_prices);

k = @(e) sqrt((1-e)/e);
%k = @(e) -1*norminv(e);
e_range  = 0.001:5*10^(-3):0.1;
% e_range = 0.0001:10^(-4):0.01

mean_vals_base = [];
sd_vals_base = [];

%init = rand(size(stock_prices,2),1);
%init = init./sum(init);
for i=1:size(e_range,2)
    
    e = e_range(1,i);
    N=size(stock_prices,2);
    S=size(stock_prices,1);
    
    
    %init=rand(N+1,1);
%     A_f=[stock_prices,ones(S,1)];
%     min_obj_f=@(y)(y(N+1)+(1/(S*e))*sum(max(A_f*y,zeros(S,1))) );
%     A=(-1)*[eye(N),zeros(N,1)];
%     b=zeros(N,1);
%     A_eq=[ones(1,N),0];
%     b_eq=1;
    
    
    
    
    
    
    options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',6000);
    options.MaxFunctionEvaluations = 40000;
    
    init=rand(N+S+2,1);
    min_obj_f = @(y) (y(N+S+2));
    A=zeros(N+2*S+1,N+S+2);
    b=zeros(N+2*S+1,1);
    
    A(1:N,1:N)=(-1)*eye(N);
    
    A(N+1,N+1)=1;
    A(N+1,N+S+2)=-1;
    A(N+1,(N+2):(N+S+1))=(1/(S*e))*ones(1,S);
    
    A((N+2):(N+S+1),N+1)=(-1)*ones(S,1);
    A((N+2):(N+S+1),(N+2):(N+S+1))=(-1)*eye(S);
    A((N+2):(N+S+1),1:N)=(-1)*stock_prices;
    
    A((N+S+2):(N+2*S+1),(N+2):(N+S+1))=(-1)*eye(S);
    
    A_eq=zeros(1,N+S+2);
    A_eq(1,1:N)=ones(1,N);
    b_eq=1;
    
    
    [y,fval,exitflag,output] = fmincon(min_obj_f,init,A,b,A_eq,b_eq,[],[],[],options);
    x=y(1:N);
    mean_vals_base = [mean_vals_base, mu'*x];
    sd_vals_base = [sd_vals_base, sqrt(x'*covariance*x)];
    
end




mark_size = 5;
F=figure(1); hold on;
box on
grid on
plot(sd_vals_base, mean_vals_base,'-s','markers',mark_size);
%plot(sd_vals_wvar, mean_vals_wvar,'-o','markers',mark_size);
lgd = legend('Base CVaR');
lgd.Location = 'southeast';
ylabel('Return');
xlabel('Standard Deviation')
saveas(F,'ef_100_cvar_cons.jpeg');
%saveas(F,'./JPEGs/bse30_simulated/ef_exact_cheb.jpeg');
%saveas(F,'./EPSs/bse30_simulated/ef_exact_cheb.eps','epsc');

hold off


risk_free = log(1.06)/365;
base = (mean_vals_base - risk_free)./sd_vals_base;
%wvar = (mean_vals_wvar - risk_free)./sd_vals_wvar;


F=figure(2); hold on;
box on
grid on
plot(e_range, base,'-o');
%plot(e_range, wvar,'-s');
lgd = legend('Base CVaR');
lgd.Location = 'southeast';
ylabel('Sharpe Ratio');
xlabel('\epsilon(Confidence level)');
saveas(F,'sr_100_cvar_cons.jpeg');

% change the names of the files and folders accordingly.
%saveas(F,'./JPEGs/bse30_simulated/sr_exact_cheb.jpeg');
%saveas(F,'./EPSs/bse30_simulated/sr_exact_cheb.eps','epsc');
hold off

