clc
clear all
close all
%matlabrc

% ParÃ¡metros de simulaciÃ³n.
J = 100e-3     % Momento de inercia
T = 20;	% Tiempo de simulacion
dt = 0.01     % Intervalo de muestreo
N = T/dt   ;     % Indice maximo para estados discretos
ts = 0:dt:T-dt; % Vector de tiempos discretos


K =1; 
Td =sqrt(2/5);

GA =  9.80665;       % Gravity constant, m/s^2

D2R = (pi/180);     % degrees to radians
R2D = (180/pi);     % radians to degrees

%% Matrices de covarianza Q y R
Q = dt^2*[0.000104 0;0 0.49]; 
R = 0.0043;               

%Matrices:
F = [1 -dt;0 1];
G = [dt;0];
L = dt;
H = [1 0]; 


A3DMGX1.arw      = 3.5   .* ones(1,3);     % Angle random walks [X Y Z] (deg/root-hour)
A3DMGX1.arrw     = zeros(1,3);           % Angle rate random walks [X Y Z] (deg/root-hour/s)
A3DMGX1.vrw      = 0.4 .* ones(1,3);     % Velocity random walks [X Y Z] (m/s/root-hour) 
A3DMGX1.vrrw     = zeros(1,3);           % Velocity rate random walks [X Y Z] (deg/root-hour/s)
A3DMGX1.gb_sta   = 0.7   .* ones(1,3);   % Gyro static biases [X Y Z] (deg/s) 
A3DMGX1.ab_sta   = 10  .* ones(1,3);     % Acc static biases [X Y Z] (mg) 
A3DMGX1.gb_dyn   = 0.02 .* ones(1,3);   % Gyro dynamic biases [X Y Z] (deg/s)
A3DMGX1.ab_dyn   = 0.2 .* ones(1,3);     % Acc dynamic biases [X Y Z] (mg)
A3DMGX1.gb_corr  = 100 .* ones(1,3);     % Gyro correlation times [X Y Z] (seconds)
A3DMGX1.ab_corr  = 100 .* ones(1,3);     % Acc correlation times [X Y Z] (seconds)
A3DMGX1.freq     = 1/dt;             % IMU operation frequency [X Y Z] (Hz)
A3DMGX1.m_psd    = 0 .* ones(1,3);   % Magnetometer noise density [X Y Z] (mgauss/root-Hz) TBD

% ref time is used to simulate IMU sensors
A3DMGX1.t = ts;                       % IMU time vector
%dt = mean(diff(A3DMGX1.t));              % IMU sampling interval

imu = imu_si_errors(A3DMGX1, dt);       % IMU manufacturer error units to SI units.

wb = gyro_gen (N, imu,dt);   % Generation of gyro in the body frame
imu.wb = wb;
fb = acc_gen(N, imu, dt);    %Generation of acc in the body frame
fb(:) = asin(fb(:)/GA);
imu.fb = fb;
AngStd=std(fb)
% Vector de estado inicial.
q0 = [0;0]; % primer indice = posición - segundo indice = vel. angular

% Vector de estados e inicializaciÃ³n.
q_True = zeros(2, N);
q_True(:, 1) = q0;
q = zeros(2, N);
q(:, 1) = q0; %todas las filas de la columna 1 las iguala a q0
% Vector de acciones de control = los torques pero en forma discreta
u = zeros(1, N);
u_True = zeros(1, N);
% Acciones de control
%u(1, 1:1:N)=1;
% SerÃ­a mejor: 
%u = ones(1, N1);


% Matrices del sistema de estados discretizado.
A = [1, dt; 0, 1];
B = [0; dt/J];

% Angulo y velocidad Ideales
v = ones(1, N);
u(1, 1) = (v(1,1)-q(1,1)+Td*(v(1,1)-q(1,1))/dt)*K;
u_True(1, 1) = (v(1,1)-q_True(1,1)+Td*(v(1,1)-q_True(1,1))/dt)*K; 
for j= 1:(N-1)
    q_True(:, j+1) = A*q_True(:, j) + B*u_True(:,j);
    u_True(:, j+1) = (v(:,j+1)-q_True(1,j+1)+Td*(v(:,j+1)-q_True(1,j+1)-(v(:,j)-q_True(1,j)))/dt)*K;
end


%Angulo y velocidad Medidas
for j= 1:(N-1)
    q(:, j+1) = A*q(:, j) + B*u(:,j);
    q(1, j+1) = q(1, j+1) + fb(j);  
    q(2, j+1) = q(2, j+1) + wb(j);  
    u(:, j+1) = (v(:,j+1)-q(1,j+1)+Td*(v(:,j+1)-q(1,j+1)-(v(:,j)-q(1,j)))/dt)*K;
end

zk=q(1,:); 
uk = q(2,:);                                                  
%Primera estimación
xh(1,1) = zk(1,1);
xh(2,1) = dt*randn(1,1)*imu.gb_psd(1);
P = eye(2);                                                                  
                                                                                                                                                            
for i = 2 : size(zk,2)                                                           
    %Prediccion
    xh(:,i) = F* xh(:,i-1) + G*uk(:,i-1);
    P_a = F*P*F' + Q; 
    Kk = P_a* H' * inv( H*P_a*H' + R ); 
    %Correccion
    xh(:,i) =  xh(:,i) + Kk .* (zk(:,i) - H*xh(:,i)); 
    P = P_a - Kk * H * P_a; 
end

    % Colores
    blue = [0, 0.4470, 0.7410];
    orange = [0.8500, 0.3250, 0.0980];
    gray= ones(1,3) * 0.75;
    
    % Line width
    lw = 1.5;

%Plot del angulo Medido
subplot(2, 1, 1);
plot(ts, q(1, :), 'color', blue,'LineWidth',lw)
legend("Posición angular con ruido")
grid on; xlabel('t [s]'); ylabel('Ángulo (\theta) [rad]');
title('Posicion angular')

%Plot de la vel  medida
subplot(2, 1, 2);
plot(ts, q(2, :), 'color', orange,'LineWidth',lw)
legend("Velocidad angular con ruido")
grid on; xlabel('t [s]'); ylabel('vel. angular (\omega) [rad/s]');
title('Velocidad angular')

%Plot de la estimación por Kalman del angulo
f2 = figure();
hold on; 
grid on;
plot(ts, q(1, :), 'color',blue,'LineWidth',lw) 
plot(ts,xh(1,:),'color',orange,'LineWidth',lw)
legend("Posición angular medida", "Posición angular estimada")
grid on; xlabel('t [s]'); ylabel('q1 (\theta)');
title('Posicion angular')
