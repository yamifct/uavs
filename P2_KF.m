% crazyflie load usd card log csv file for sensor analysis
clear all

% available files:
csvfilename_motors_off = 'L2Data1.csv';
csvfilename_hover = 'L2Data2.csv';
csvfilename_motion_x = 'L2Data3.csv';
csvfilename_motion_y = 'L2Data4.csv';
csvfilename_motion_z = 'L2Data5.csv';
csvfilename_motion_inf = 'L2Data6.csv';

% read file
csvfilename =csvfilename_motion_z ;
array = dlmread(csvfilename,',',1,0);
%T = table2array(readtable(csvfilename)); % Matlab only

g=9.81;

% get data from table
time = array(:,1)'*1e-3;
pos = array(:,2:4)'; % [m]
vel = array(:,5:7)'; % [m/s]
lbd = array(:,8:10)'; % [deg]
gyro = array(:,11:13)'; % [deg/s]
acc = array(:,14:16)'*g; % [m/s^2]
baro_asl = array(:,17)'; % [m]
% lighthouse = array(:,18:20))'; % [m]

% convert date to print format
t = time - time(1);

% plot data
initPlots;
%crazyflie_show_sensor_data(t,acc,gyro,baro_asl,pos,vel,lbd);


%Filtro de Kalman Linear
%% Analizar L2Data5

g = 9.81; 

% Matrizes do Sistema
A = [0 0  1  1  0;
     0 0  0  0  0;
     0 0  0  0  0;
     0 0  0  0  0;
     0 0  0  0  0];

C = [-g  0  0  0  0;
      0  g  0  0  0;
      0  0  1  0  1;
      0  0  0  1  0];

% Carregamento dos Dados 


% Ruídos do processo e da medição
Q = 0.01 * eye(5); % Covariância do ruído do processo
cov_ax=cov(acc(1,:));
cov_ay=cov(acc(2,:));
cov_az=cov(acc(3,:));
cov_gx=cov(gyro(1,:));
cov_gy=cov(gyro(2,:));
cov_gz=cov(gyro(3,:));
R=diag([cov_ax, cov_ay, cov_gx, cov_gy]); % Covariância do ruído da medição 

% Medidas 
y=[acc(1,:)
   acc(2,:)
   gyro(1,:)
   gyro(2,:)];
% Estado inicial e covariância inicial
x_est = zeros(5, 1); % Estado inicial (zero)
P = eye(5); % Covariância inicial do estado
N=length(acc);

for k = 1:N
    % Predição
    x_pred = A * x_est;
    P_pred = A * P * A' + Q;
    
    % Atualização
    K = P_pred * C' * inv(C * P_pred * C' + R);
    x_est = x_pred + K * (y(:, k) - C * x_pred);
    P = (eye(5) - K * C) * P_pred;
    
    estimates(:, k) = x_est; % Armazenamento das estimativas
end

figure;
subplot(2, 1, 1);
plot(estimates(1, :), 'b', 'LineWidth', 2);
hold on;
plot(lbd(2,:)*pi/180, 'r');
title('Pitch-L2Data5');
legend('Estimado - KF', 'Real');

subplot(2, 1, 2);
plot(estimates(2,:), 'b', 'LineWidth', 2);
hold on;
plot(lbd(1, :)*pi/180, 'r');
title('Roll-L2Data5');
legend('Estimado - KF', 'Real');
%Filtro de Kalman estendido


%Tensor de inércia - encontrado na parte 1
J =(10^-9)*[ 10910.69, 0, 0 
                0, 11037.19, 0
                0, 0, 21098.42];
JJ=diag(J);
J1=JJ(1);
J2=JJ(2);
J3=JJ(3);

% Matrizes do Sistema
A = [0 0  1  1  0;
     0 0  0  0  0;
     0 0  0  0  0;
     0 0  0  0  0;
     0 0  0  0  0];

Bw = eye(5);
D=eye(6);

Q = 0.01 * eye(5); % Covariância do ruído do processo

R=diag([cov_ax, cov_ay, cov_az, cov_gx, cov_gy, cov_gz]); % Covariância do ruído da medição
% Estado inicial e covariância inicial
x_est = zeros(5, 1); % Estado inicial (zero)
x_pred=x_est;

P = 1*eye(5); % Covariância inicial do estado
N=length(acc);


% Implementação do Filtro de Kalman

y=[acc(1,:)
   acc(2,:)
   acc(3,:)
   gyro(1,:)
   gyro(2,:)
   gyro(3,:)];

for k = 1:N
    p=gyro(1,k);
    q=gyro(2,k);
    r=gyro(3,k);

   theta=asin(-acc(1,k)/sqrt(acc(1,k)^2+acc(2,k)^2+acc(3,k)^2));
   phi=atan2(acc(1,k),acc(3,k));

    A  = [0 -q*sin(phi)-r*cos(phi) cos(phi) 0 -sin(phi)
          (q*sin(phi)+r*cos(phi))*sec(theta)^2 (q*cos(phi)-r*sin(phi))*tan(theta) sin(phi)*tan(theta) 1 cos(phi)*tan(theta)
          0 0 0 r/J2*(J3-J1) p/J2*(J3-J1)
          0 0 r/J1*(J2-J3) 0 q/J1*(J2-J3)
          0 0 1/J3*(J1*p - r*J2) q/J3*J1 q/J3*J2]; 


   C=[-9.81*cos(theta) 0 0 0 0
    -9.81*sin(theta)*sin(phi) 9.81*cos(theta)*cosd(phi) 0 0 0
    -9.81*sin(theta)*cos(phi) -9.81*cos(theta)*sind(phi) 0 0 0
    0 0 1 0 0
    0 0 0 1 0
    0 0 0 0 1];

    % Predição
    x_pred = A * x_est;
    P_pred = A * P * A' + Bw*Q*Bw';
    % Atualização
    K = P_pred * C' *(C * P_pred * C' + D*R*D')^-1;
    P = (eye(5) - K * C) * P_pred;
    x_est = x_pred + K * (y(:,k) - C * x_pred);

    estimates1(:, k) = x_est; % Armazenamento das estimativas
end


figure;
subplot(2, 1, 1);
plot(estimates1(1, :), 'b', 'LineWidth', 2);
hold on;
plot(lbd(2,:)*pi/180, 'r');
title('Pitch-L2Data5');
legend('Estimado - EKF', 'Real');

subplot(2, 1, 2);
plot(estimates1(2,:), 'b', 'LineWidth', 2);
hold on;
plot(lbd(1, :)*pi/180, 'r');
title('Roll-L2Data5');
legend('Estimado - EKF', 'Real');




