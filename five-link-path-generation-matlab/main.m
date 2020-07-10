clear
clc
global N m1 m2 m3 m4 m5 l1 l2 l4 l3 l5 timer init_joint_angles T step_time
N = 10;m1 = 1;m2 = 1;m3 = 1;m4 = 1;m5 = 1;l1 = 0.4;l2 = 0.4;l5 = 0.4;l3=0.6;l4=0.4; timer = 0;
tic
T = 1;

dmax = 20;
umax = 10;
d = 1;
h = T/(N);
t = 0:h:T;
step_time=length(t);
t0=t;
% t = t(1:(length(t)-1));
init_joint_angles = [...
    -0.3; % stance leg tibia angle
    0.7; % stance leg femur angle
    0.0; % torso angle
    -0.5; % swing leg femur angle
    -0.6]; % swing leg tibia angle

 %init_joint_angles = [-0.4889;-0.1571;0.2095;-0.5063;-0.0174];
fini_joint_angles = [init_joint_angles([5;4;3;2;1])];...init_joint_angles([10;9;8;7;6])]; 
% fini_joint_angles = [-deg2rad(28);-deg2rad(9);deg2rad(12);-deg2rad(29);-deg2rad(1) ...
%                      ;0.1;-1.09;0;-0.16;1.03];
% init_joint_angles = [-0.3;0.7;0;-0.5;-0.6; 0;0;0;0;0];
% fini_joint_angles = [0.6;-0.5;0;0.7;-0.3; 0;0;0;0;0];
x0=initial_guess(t,init_joint_angles,fini_joint_angles);
x = initial_guess(t,init_joint_angles,fini_joint_angles);
xf= zeros(15,11);
% x = [(0:N-1)*(d)*h/T;(0:N-1)*pi*h/T, zeros(1,3*N)];
   
    lb = [...-dmax*ones(1,N) ; 
     (-pi/2)*ones(5,N+1) ;  -30*ones(5,N+1); zeros(1,N+1);-umax*ones(4,N+1)];
    ub = [... dmax*ones(1,N) ;
      (pi/2)*ones(5,N+1) ; 30*ones(5,N+1) ; zeros(1,N+1);umax*ones(4,N+1)];  
  options = optimoptions('fmincon','MaxFunctionEvaluations',5000,'Display','iter-detailed','PlotFcn',@optimplotconstrviolation);
    x = fmincon(@callme,x,[],[],[],[],lb,ub,@colons,options);%function handle is paper wala f(x)
    xf=xf+x;
  for i = 1:1:9 %steps
    options = optimoptions('fmincon','MaxFunctionEvaluations',5000,'Display','iter-detailed','PlotFcn',@optimplotconstrviolation);
    x = fmincon(@callme,x,[],[],[],[],lb,ub,@colons,options);%function handle is paper wala f(x)
    xn=[x(1:5,end);x(6:10,end)];
    xp=[x(1:5,end-1);x(6:10,end-1)];
    xt = impact(xn,xp);
    x(6:10,end) = xt;
    
%          xt = x;   
%          g=x(1,:);
%          x(1,:)=x(5,:);
%          x(5,:)=g;
%          g=x(2,:);
%          x(2,:)=x(4,:);
%          x(4,:)=g;
%          g=x(6,:);
%          x(6,:)=x(10,:);
%          x(10,:)=g;
%          g=x(7,:);
%          x(7,:)=x(9,:);
%          x(9,:)=g;
%          g=x(11,:);
%          x(11,:)=x(15,:);
%          x(15,:)=g;
%          g=x(12,:);
%          x(12,:)=x(14,:);
%          x(14,:)=g;
          xf=[xf,x];
    t=[t,(i+t0)];
  end
toc

% % subplot(3,1,1);
% % plot(t,x(1:5,:));
% % subplot(3,1,2);
% % plot(t,x(6:10,:));
% % subplot(3,1,3);
% % plot(t,x(11:14,:));
% % subplot(3,1,1);
% % plot(t,quadratic_spline(x(1:N),x(2*N+1:3*N),t));
% % subplot(3,1,2);
% % plot(t,quadratic_spline(x(N+1:2*N),x(3*N+1:4*N),t));
% % subplot(3,1,3);
% % plot(t,linear_spline(x(4*N+1:5*N),t));
% % %stem(x(4*N+1:5*N));