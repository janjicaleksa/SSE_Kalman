clear
close all
clc
%% a
load('dom2_zad1.csv');
opservations_pom = dom2_zad1;
x = opservations_pom(1,:); %prvi eksperiment

Nr = 100; 
N = 10;

syms theta
logp = - sum(x-theta) - sum(exp(-(x-theta))); %log verodostojnost
dlogp_analiticki = N - sum(exp(-(x-theta))); % prvi izvod analiticki
dlogp_numericki = diff(logp,theta); % prvi izvod numericki
d2logp_analiticki = - sum(exp(-(x-theta))); % drugi izvod analiticki
d2logp_numericki = diff(dlogp_numericki,theta); % drugi izvod numericki

figure(1)
    fplot(logp);
    title('Log-verodostojnost'); xlabel('\theta');
    grid on
figure(2)
    fplot(dlogp_analiticki);
    hold all;
    fplot(dlogp_numericki);
    title('Prvi izvod log-verodostojnosti'); xlabel('\theta');
    grid on;
    legend('analiticki','numericki','Location','Best');
figure(3)
    fplot(d2logp_analiticki);
    hold on;
    fplot(d2logp_numericki);
    title('Drugi izvod log-verodostojnosti'); xlabel('\theta');
    grid on;
    legend('analiticki','numericki','Location','Best');

%% c    
estimations = zeros(1,100);

for i = 1:Nr
    x = opservations_pom(i,:);
    logp = - sum(x-theta) - sum(exp(-(x-theta)));
    dlogp = N - sum(exp(-(x-theta)));
    d2logp = - sum(exp(-(x-theta)));
    
    estimations(i) = Newton_method(2,dlogp,d2logp,10^(-3),20);
end

mi_est = mean(estimations);
sigma_est = sqrt(var(estimations)); 
xosa=0:0.01:5;
y = normpdf(xosa,mi_est,sigma_est);

figure(4)
plot(xosa,y,'k-','LineWidth',2);
hold all;
histogram(estimations,'Normalization','pdf');
xline(mi_est,'r--','LineWidth',2)
xline(mi_est - sigma_est,'g--','LineWidth',2)
xline(mi_est + sigma_est,'b--','LineWidth',2)
grid on;
legend('normalna fgv','histogram','\mu','\mu - \sigma','\mu + \sigma');