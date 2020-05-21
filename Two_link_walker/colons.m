function [c,ceq] = colons(x)
    T = 1;
    N = 50;
    h = T/N;
    l1 = 1;
    l2 = 1;
    c = [];
    theta1_ini = 110;
    theta2_ini = 220;
    theta1_final = 70;
    theta2_final = 140;
    ceq = [x(1,1)-theta1_ini*(pi/180); x(2,1)-theta2_ini*(pi/180); x(3,1); x(4,1); x(1,N+1)-theta1_final*(pi/180); x(2,N+1)-theta2_final*(pi/180); x(3,N+1); x(4,N+1)];
    for k = 1:N
        c = [c; -(l1*sin(x(1,k))+l2*sin(x(1,k)-x(2,k)))]; %condition for swing point to be above ground
        state = [(x(1:2,k)); x(3:4,k)];
        nextState = [(x(1:2,k+1)); x(3:4,k+1)];
        tau = x(5,k);
        tauNext = x(5,k+1);
        xdot = dynamics(state,tau);
        xdotNext = dynamics(nextState,tauNext);
        ceq = [ceq; nextState - state - ((h/2)*(xdot + xdotNext))];
    end
    