clear all;0
close all;
clc;

y = load("gnss_data.csv");

T = 1;
A = [1 T T^2/2;0 1 T;0 0 1];
B = [T^2/2; T; 0];
U = [0;0;1];
H = [1 0 0];

mi_a = [5*ones(1,30) 0*ones(1,40) -5*ones(1,30)];
sigma_a = [5*ones(1,30) 1*ones(1,40) 5*ones(1,30)];
var_a = (sigma_a/3).^2; %smatramo 99.7% sigurnost

Q = U*U'*var_a(1); %kovarijaciona matrica modela

sigma_u_poc = 10;
C = sigma_u_poc^2*ones(1,length(y)); %vektor varijanse suma merenja
for i=1:length(y)
    if(i>=40) && (i<45)
        sigma_u_poc = sigma_u_poc + 10;
        C(i) = sigma_u_poc^2;
    end
    if(i>=45) && (i<=55)
        C(i) = 10^12;
    end
    if(i>55) && (i<=60)
        sigma_u_poc = sigma_u_poc - 10;
        C(i) = sigma_u_poc^2;
    end
end

%inicijalizacija i pocetna estimacija
s_est = [0; 0; 0];
M_est = eye(3); %prilicno smo sigurni da se krece iz nulte tacke,...
                %sa nultom pocetnom brzinom i nultim odstupanjem
                %ubrzanja

s_estimirano = zeros(3, length(y) + 1); %matrica vektora estimiranih stanja
M_estimirano = zeros(3, length(y) + 1); %matrica vektora kovarijansi greske estimacije
K_pojacanje = zeros(3,length(y)); %matrica vektora Kalmanovog pojacanja

s_estimirano(:,1) = s_est;
M_estimirano(:,1) = [M_est(1,1);M_est(2,2);M_est(3,3)];

for i = 1:length(y)
    %predikcija
    s_pred = A*s_est + B*mi_a(i);
    M_pred = A*M_est*A' + U*U'*var_a(i);  
    
    %estimacija
    if (i>=45) && (i<=55)
        K = 0;
        s_est = s_pred;
        M_est = M_pred;
    else
        K = M_pred*H'*inv(H*M_pred*H' + C(i));
        s_est = s_pred + K*(y(i) - H*s_pred);
        M_est = (eye(3) - K*H)*M_pred;
    end
    
    K_pojacanje(:,i) = K;
    s_estimirano(:, i+1) = s_est;
    M_estimirano(:, i+1) = [M_est(1,1);M_est(2,2);M_est(3,3)];
end

M_2sigma = 2*sqrt(M_estimirano);
t = 1:100;

figure(1)
plot (t, y);
hold all;
plot(t, s_estimirano(1,2:end));
hold all; 
plot(t, s_estimirano(1,2:end)+M_2sigma(1,2:end));
hold all;
plot(t, s_estimirano(1,2:end)-M_2sigma(1,2:end));
grid on;
xlabel('t[s]');
ylabel('x[m]');
title('Pozicija');
legend('stvarna','estimirana','estimirana + 2\sigma','estimirana - 2\sigma')

figure(2)
plot(t, s_estimirano(2,2:end))
hold all;
plot(t, s_estimirano(2,2:end)+M_2sigma(2,2:end));
hold all;
plot(t, s_estimirano(2,2:end)-M_2sigma(2,2:end));
grid on;
xlabel('t[s]');
ylabel('v[m/s]');
title('Brzina');
legend('estimirana','estimirana + 2\sigma','estimirana - 2\sigma')

figure(3)
plot(t, s_estimirano(3,2:end)+mi_a(1,:));
hold all;
plot(t, s_estimirano(3,2:end)+mi_a(1,:)+M_2sigma(3,2:end));
hold all;
plot(t, s_estimirano(3,2:end)+mi_a(1,:)-M_2sigma(3,2:end));
grid on;
xlabel('t[s]');
ylabel('a[m/s^2]');
title('Ubrzanje');
legend('estimirano','estimirano + 2\sigma','estimirano - 2\sigma')

figure(4)
plot(t, K_pojacanje(1,:));
hold all;
plot(t,K_pojacanje(2,:));
hold all
plot(t,K_pojacanje(3,:));
grid on;
xlabel('t[s]');
ylabel('K');
title('Kalmanovo pojacanje');
legend('pozicija','brzina','ubrzanje');