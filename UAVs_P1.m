%% Project UAVs - Crazyflie drone control

clear 
clc
initPlots;
 
% Earth parameters  
g_T = 9.81; 
beta = 0.1;

% Drone details
m = 0.027;
J =(10^-9)*[ 14854.94, 0, 0 
                0, 14850.49, 0
                0, 0, 28625.5]; %Values calculated with solidworks
zI = [0;0;1];

    %Thrusts positions
           l = 5/100;
           lz = 0;
           alpha = 45;
           p1 = [l*cosd(alpha);-l*cosd(alpha);lz];
           p2 = [-l*cosd(alpha);-l*cosd(alpha);lz];
           p3 = [-l*cosd(alpha);l*cosd(alpha);lz];
           p4 = [l*cosd(alpha);l*cosd(alpha);lz];

%1. Crazyflie modeling

% simulate nonlinear system
Dt = 0.01;
t = 0:Dt:10;
Nsim = length(t);
nx = 12;
ny = 4; 
x0 = zeros(nx,1);
x = zeros(nx,Nsim);
y = zeros(ny,Nsim);
a=zeros(3,Nsim); %auxiliary variable
fa=a;np=a;fp=a;f1=a;f2=a;f3=a;f4=a;n1=a;n2=a;n3=a;n4=a;

for k=1:Nsim
    % prepare variables:
    p  = x(1:3,k);
    v   = x(4:6,k);
    lbd = x(7:9,k);
    om  = x(10:12,k);
    R = Euler2R(lbd);

    % grandezas - Ver referenciais;
    % Força de arrasto
    d1 = beta;
    d2 = beta;
    d3 = beta;
    Dd = diag([d1, d2,d3]);
    fa(:,k) = - R*Dd*R'*v; 
    
    %Força propulsiva e momento
    fp(:,k) =  m*g_T*R'*zI;
    f1(:,k) = fp(:,k)/4;
    f2(:,k) = f1(:,k);
    f3(:,k) = f1(:,k); 
    f4(:,k) = f1(:,k);
    T1 = f1(3); T2 = f2(3); T3 = f3(3); T4 = f4(3);     
                
           % moments (rotor)
           cQ = 0.1;
           cT = 5;
           n1(:,k) = -cQ*(T1/cT)*[0;0;1]; 
           n2(:,k) = cQ*(T2/cT)*[0;0;1];
           n3(:,k) = -cQ*(T3/cT)*[0;0;1];
           n4(:,k) = cQ*(T4/cT)*[0;0;1];

     % moments
     nf1 = n1+cross(p1,f1(:,k));
     nf2 = n2+cross(p2,f2(:,k));
     nf3 = n3+cross(p3,f3(:,k));
     nf4 = n4+cross(p4,f4(:,k));
     np(:,k) = nf1(:,k) + nf2(:,k) + nf3(:,k) + nf4(:,k);
    
    % compute state derivative:
    p_dot = R*v;
    v_dot = -skew(om)*v + g_T*R'*zI + fa + 1/m*(fp(:,k));
    lbd_dot = Euler2Q(lbd)*om;
    om_dot = -J\skew(om)*J*om + J\(np(:,k));
    
    x_dot = [p_dot;v_dot(:,k);lbd_dot;om_dot];

    % integrate state
    x(:,k+1) = x(:,k) + Dt*x_dot; 

    % state variable computation
    v(1:3,k)   = x(4:6,k);
    lbd(1:3,k) = x(7:9,k);
    om(1:3,k)  = x(10:12,k);
    % compute current output:
    C = [   eye(3)    , zeros(3)   , zeros(3) , zeros(3)
        zeros(1,3), zeros(1,3) , zI'      , zeros(1,3)  ];
  
    y(:,k) = C*x(:,k);
end   

% plot results
figure(1);
plot(t,y(1,:),'-.',t,y(2,:),'--',t,y(3,:),'-',t,y(4,:),':');
title('Nonlinear model simulation')
xlabel('$$Time [s]$$');
ylabel('$$Output (t)$$');
grid on;
legend('$$ x(t)$$',...
        '$$ y(t)$$',...
        '$$ z(t)$$',...
        '$$ yaw(t)$$',...
       'Location','northwest');


%   Linearization

% Hover conditions
A1 = [   zeros(3), eye(3)     , zeros(3)  , zeros(3)
        zeros(3), -Dd         , skew(g_T*zI), zeros(3)
        zeros(3), zeros(3)   , zeros(3)  , eye(3)
        zeros(3), zeros(3)   , zeros(3)  , zeros(3)  ];

B = [       zeros(3,1),      zeros(3,3);
            1/m*zI,     zeros(3,3);
            zeros(3,1), zeros(3,3);
            zeros(3,1), inv(J)];

D = zeros(4);
% Results
sys1 = ss(A1,B,C,D);
Te  = g_T*m;
u_L1 = [Te;0.01;0.01;0.01]*(t>=0);
y_L1 = lsim(sys1,u_L1',t,x0)';

figure(2);
plot(t,y_L1(1,:),'-.',t,y_L1(2,:),'--',t,y_L1(3,:),'-',t,y_L1(4,:),':');
title('Hover linearized model simulation')
xlabel('$$Time [s]$$');
ylabel('$$Output (t)$$');
grid on;
legend('$$  xL(t)$$',...
        '$$ yL(t)$$',...
        '$$ zL(t)$$',...
        '$$ yaw(t)$$',...
       'Location','northwest');

%Decomposition into modes - Compute left and right eigenvectors 
[V_1,D_1,W_1] = eig(A1);

%Since there are 12 state-space variables, there are 12 modes.
C_t1 = zeros(ny,nx);
O_b1 = zeros(ny,nx);
for i=1:nx
    C_t1(:,i)= W_1(:,i)'*B;
    O_b1(:,i)= C*V_1(:,i); 
end

%Check controlability of each mode
disp('1st linearized model: modes controlability check')
disp(' ');
for i=1:nx   
% Check if all elements of the vector are zeros
    zc1 = all(C_t1(:,i) == 0);
       % Display result
        if zc1
            disp("Mode "+ num2str(i) + " from OP1 linearized model is uncontrolable");
        else
            disp("Mode "+ num2str(i) + " from OP1 linearized model is controlable");
        end
end

%Check observability of each mode
disp(' ');
disp('1st linearized model: modes observability check')
disp(' ');
for i=1:nx   
% Check if all elements of the vector are zeros
     zo1 = all(O_b1(:,i) == 0);
    % Display result
        if zo1
           disp("Mode "+ num2str(i)+ " from OP1 linearized model is unobservable");

        else
            disp("Mode "+ num2str(i)+ " from OP1 linearized model is observable");
        end
 end
disp(' ');
% Stability, controlability and observability 
disp('1st linearized model: model stability, observability and controlability check')
disp(' ');
%Stability
        % Create a vector with all eigenvalues
        D1 = diag(D_1);
        % Find unique values in the vector
        un_eg1 = unique(D1);
        % Compare the length of the D1 vector with the length of the unique values
        if length(D1) ~= length(un_eg1)
            disp('Model 1 is unstable');
        else
            disp('Model 1 is stable');
        end

% Observability
        mode_obs1 = obsv(A1,C);
        if rank(mode_obs1) == size(A1,1)
          disp('Model 1 is observable');
        else
           disp('Model 1 is unobservable')  
        end

 %Controlability
           
        mode_ctrl1 = ctrb(A1,B);
        if rank(mode_ctrl1)==size(A1,1)
          disp('Model 1 is controllable');
        else
           disp('Model 1 is uncontrollable') 
        end

% Horizontal flight conditions
    % Auxiliary Matrices and variables
    theta_eq = 5*pi/180;
    Vex = 5;
    Ve = [Vex;0;0];
    Ry_e=[cos(theta_eq), 0, sin(theta_eq)
          0,1, 0
          -sin(theta_eq), 0, cos(theta_eq)];
    
    Pl = [0, -sin(theta_eq)*Vex, 0
          0, 0, cos(theta_eq)*Vex
          0, -cos(theta_eq)*Vex,0];
    
    Nv = [0, -cos(theta_eq)*g_T, 0 
         g_T*cos(theta_eq), 0, 0
          0, g_T*sin(theta_eq),0];
    
    Ql = [1,0,0
          0,1,0
          0,0,1/cos(theta_eq)];
    
    % A matrix
    A2 = [      zeros(3), Ry_e,                      Pl,     zeros(3)
                zeros(3), -beta*eye(3),   Nv,  skew(Ve)
                zeros(3), zeros(3),                     zeros(3),      Ql
                zeros(3), zeros(3),                     zeros(3),      zeros(3) ];
% Results
sys2 = ss(A2,B,C,D);
u_L2 = [Te;0;0.1;0]*(t>=0);
y_L2 = lsim(sys2,u_L2',t,x0)';

% plot results
figure(3);
plot(t,y_L2(1,:),'-.',t,y_L2(2,:),'--',t,y_L2(3,:),'-',t,y_L2(4,:),':');
title('Horizontal light Linearized model simulation')
xlabel('$$Time [s]$$');
ylabel('$$Output (t)$$');
grid on;
legend('$$  xL(t)$$',...
        '$$ yL(t)$$',...
        '$$ zL(t)$$',...
        '$$ yaw(t)$$',...
       'Location','northwest');

%Decomposition into modes - Compute left and right eigenvectors 
[V_2,D_2,W_2] = eig(A2);

%Since there are 12 state-space variables, there are 12 modes.
C_t2 = zeros(ny,nx);
O_b2 = zeros(ny,nx);
for i=1:nx
    C_t2(:,i)= W_2(:,i)'*B;
    O_b2(:,i)= C*V_2(:,i); 
end

%Check controlability of each mode
disp('2nd linearized model: modes controlability check')
disp(' ');
for i=1:nx   
% Check if all elements of the vector are zeros
    zc2 = all(C_t2(:,i) == 0);
       % Display result
        if zc2
            disp("Mode "+ num2str(i) + " from OP2 linearized model is uncontrolable");
        else
            disp("Mode "+ num2str(i) + " from OP2 linearized model is controlable");
        end
end

%Check observability of each mode
disp(' ');
disp('2nd linearized model: modes observability check')
disp(' ');
for i=1:nx   
% Check if all elements of the vector are zeros
     zo2 = all(O_b2(:,i) == 0);
    % Display result
        if zo2
           disp("Mode "+ num2str(i)+ " from OP2 linearized model is unobservable");%, ' from? the OP2 linearized mode is uncontrolable'

        else
            disp("Mode "+ num2str(i)+ " from OP2 linearized model is observable");
        end
 end
disp(' ');
% Stability, controlability and observability 
disp('2nd linearized model: model stability, observability and controlability check')
disp(' ');
%Stability
        % Create a vector with all eigenvalues
        D2 = diag(D_2);
        % Find unique values in the vector
        un_eg2 = unique(D2);
        % Compare the length of the D2 vector with the length of the unique values
        if length(D2) ~= length(un_eg2)
            disp('Model 2 is unstable');
        else
            disp('Model 2 is stable');
        end

% Observability
        mode_obs2 = obsv(A2,C);
        if rank(mode_obs2) == size(A2,1)
          disp('Model 2 is observable');
        else
           disp('Model 2 is unobservable')  
        end

 %Controlability
           
        mode_ctrl2 = ctrb(A2,B);
        if rank(mode_ctrl2)==size(A2,1)
          disp('Model 2 is controllable');
        else
           disp('Model 2 is uncontrollable') 
        end

   % Closed loop stability with root-locus
   G_T_Pz = tf(1/m,[1 d3 0]);
   figure (4)
   rlocus(G_T_Pz) 
   title('From total thrust to z position')
   G_Phi_Py =  tf(g_T,[1 d2 0]);
   figure (5)
   rlocus(G_Phi_Py) 
    title('From roll angle to y position')



