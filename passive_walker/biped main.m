clear
clc

    global p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 psi;
    
    mH = 10;
    m1 = 5;
    m2 = 5;
    a = 0.5;
    b = 0.5;
    l = a + b;
    g= 9.81;
    
    p1 = mH*l^2 + m1*b^2 + m2*l^2;
    p2 = m2*l*a;
    p3 = m2*a^2;
    p4 = (m1*b + m2*l + mH*l)*g;
    p5 = m2*a*g;
    
    p6 = m1*l^2 + mH*l^2 + m2*b^2;
    p7 = m1*a*l;
    p8 = m1*a^2;
    p9 = m1*a*b;
    p10 = m2*b*l + mH*l^2 + m1*b*l;
    p11 = m2*b*a;
    
    psi = 4.1255*pi/180;
    %psi = 3.1*pi/180;
    
q0 = [(12.53*pi/180); (-18.53*pi/180)];
Dq0 = [-1;1.5];

x0 = [q0; Dq0];

MAP_SELECTOR = 1;

SIMULATION_TIME = 15;
DRAW_INTERVAL =  0.05;
T0=0;
tspan=[T0 SIMULATION_TIME];

state_space = [];
time = [];
impacts = [];
last_impact = 0;
gait_period = [];

current_time = T0;

while (current_time < SIMULATION_TIME)
    
    options = odeset('Events', @impact_event);
    [tout, xout, event_time, event_state, event_id] = ode45(...
        @equations_of_motion, tspan, x0, options);
    xout(:,1)=wrapToPi(xout(:,1));
    xout(:,2)=wrapToPi(xout(:,2));
    disp(xout)
    if ~isempty(event_id) && event_time(end) == tout(end)
        impact_time = tout(end);
        impact_index = length(time) + length(tout);
        impacts = [impacts; impact_time, impact_index];
        gait_period = [gait_period; impact_time - last_impact];
        last_impact = impact_time;
        x0 = impact_map(xout(end,:));
        
        %fprintf('**************\n');
        %fprintf('Impact at time = %0.2f\n', impact_time);
        %fprintf('Index = %d\n', impact_index);
        %fprintf('q1 = %0.2f\n', xout(end,1)*180/pi);
        %fprintf('q2 = %0.2f\n', xout(end,2)*180/pi);
    end
    
    time = [time; tout];
    state_space = [state_space; xout];
    
    tspan = [time(end) max(SIMULATION_TIME, time(end))];
    current_time = time(end);
end
%x0=wrapToPi(x0);

plot(state_space(:,1),state_space(:,3),state_space(:,2),state_space(:,4))
%animate_walker(time, impacts, DRAW_INTERVAL, state_space, figure);
