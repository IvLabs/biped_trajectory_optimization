function x = initial_guess(t,i0,f0)
global N timer%m1 m2 m5 l1 l2 l5;
k = 0;

% r0c1=[l1*sin(theta1);l1*cos(theta1)]/2;
% r1c2=[l2*sin(theta1+theta2);l2*cos(theta1+theta2)]/2;
% r2c3=[l2*sin(theta1+theta2+theta3);l2*cos(theta1+theta2+theta3)]/2;
% r3c4=[l1*sin(theta1+theta2+theta3+theta4);l1*cos(theta1+theta2+theta3+theta4)]/2;
% r2c5=[l5*sin(theta1+theta2+theta5);l5*cos(theta1+theta2+theta5)]/2;
%     
% r1c1=-r0c1;
% r2c2=-r1c2;
% r3c3=-r2c3;
% 
% r0c1=[l1*sin(ftheta1);l1*cos(ftheta1)]/2;
% r1c2=[l2*sin(ftheta1+ftheta2);l2*cos(ftheta1+ftheta2)]/2;
% r2c3=[l2*sin(ftheta1+ftheta2+ftheta3);l2*cos(ftheta1+ftheta2+ftheta3)]/2;
% r3c4=[l1*sin(ftheta1+ftheta2+ftheta3+ftheta4);l1*cos(ftheta1+ftheta2+ftheta3+ftheta4)]/2;
% r2c5=[l5*sin(ftheta1+ftheta2+ftheta5);l5*cos(ftheta1+ftheta2+ftheta5)]/2;
% 

% xcom = ((r0c1*9) + (r1c3*7) + (r2c3*3) + (r3c4) + (r2c5))/5;
% con_torque = k*ones(5,1)*t;
xi = i0(1:5);
xf = f0(1:5);
% x0 = (xi) + (((xf-xi)/2)*(t.^2)) + (((xf-xi)/2)*t);
x0 = (xi)+((xf-xi)*t);
% xi_dot = i0(6:10);
% xf_dot = f0(6:10);

% x0_dot = xi_dot + (xf_dot - xi_dot)*t;
x0_dot = (xf-xi)*ones(1,N+1);
% con_torque = (xf_dot-xi_dot)*ones(1,N+1);
con_torque = zeros(5,N+1);
x = [x0;x0_dot;con_torque]; 
% timer = timer + toc
end