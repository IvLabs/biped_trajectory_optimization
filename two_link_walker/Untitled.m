clc
clear
N = 50;
T = 1;
tic
h = T/N;
y = 0:h:T;
i = 0:N;
teta = 0:(2*pi/N):2*pi;
theta1_ini = 110;
theta2_ini = 220;
theta1_final = 70;
theta2_final = 140;
%x0 = [ ((theta1_ini-i*((theta1_ini-theta1_final)/N))*(pi/180)).*sin(teta); ((theta2_ini-i*((theta2_ini-theta2_final)/N))*(pi/180)).*sin(teta); zeros(1,N+1); zeros(1,N+1); zeros(1,N+1)];
x0 = [ theta1_ini*sin(teta) ; theta2_ini*sin(teta);zeros(1,N+1); zeros(1,N+1); zeros(1,N+1)] ;
umax = 250;
options = optimoptions('fmincon','MaxFunctionEvaluations',10000);
%options = optimoptions('fmincon','MaxFunctionEvaluations',10000,'Display','iter-detailed','PlotFcn',@optimplotconstrviolation);
lb = [ zeros(1,N+1); 0*ones(1,N+1);  -inf*ones(1,N+1); -inf*ones(1,N+1); -umax*ones(1,N+1)];
ub = [ pi*ones(1,N+1); 2*pi*ones(1,N+1); inf*ones(1,N+1); inf*ones(1,N+1); umax*ones(1,N+1)];
x = fmincon(@callme,x0,[],[],[],[],lb,ub,@colons,options);

for k = 1:N+1
    if abs(sin(x(1,k)+x(2,k)))>0.99998
        t = x(1,k+1:end);
        x(1,k+1:end) = x(2,k+1:end);
        x(2,k+1:end) = x(2,k+1:end)-t;
        t = x(3,k+1:end);
        x(3,k+1:end) = x(4,k+1:end);
        x(4,k+1:end) = x(4,k+1:end) - t;
        %x(5,k+1:end) = -x(5,k+1:end);
        display(k);
    end
end
        

%{
figure
plot(y,x(1,:));
xlabel('time');
ylabel('Theta 1');


figure
plot(y,x(2,:));
xlabel('time');
ylabel('Theta2');
%}
% animation
l1 = 1;
l2 = 1;
fixedPoint = [0, 0];
axis(gca,'equal');
axis([-2 2 -1 2]);

for t = 1:N+1
    H = [l1*cos(x(1,t)), l1*sin(x(1,t))];
    swingPoint = [l1*cos(x(1,t))+l2*cos(x(1,t)-x(2,t)),l1*sin(x(1,t))+l2*sin(x(1,t)-x(2,t))];
    stanceLeg = line([fixedPoint(1) H(1)], [fixedPoint(2) H(2)]);
    swingLeg = line([H(1) swingPoint(1)], [H(2) swingPoint(2)]);
    horizontalLine = line([-2 2], [0 0]);
    pause(0.1);
    delete(stanceLeg);
    delete(swingLeg);
end

toc
