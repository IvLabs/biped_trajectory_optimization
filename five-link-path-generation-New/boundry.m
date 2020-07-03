function [c,ceq] = boundry(x)

%initial and final state ground
global N m1 m2 m3 m4 m5 l1 l2 l4 l3 l5 timer init_joint_angles

c = [];

r01=[-l1*sin(x(1,1));l1*cos(x(1,1))];
r12=r01+[-l2*sin(x(2,1));l2*cos(x(2,1))];
r23=r12+[l3*sin(x(3,1));l3*cos(x(3,1))];
r24=r12+[l4*sin(x(4,1));-l4*cos(x(4,1))];
r45=r24+[l5*sin(x(5,1));-l5*cos(x(5,1))];

p01=[-l1*sin(x(1,end));l1*cos(x(1,end))];
p12=p01+[-l2*sin(x(2,end));l2*cos(x(2,end))];
p23=p12+[l3*sin(x(3,end));l3*cos(x(3,end))];
p24=p12+[l4*sin(x(4,end));-l4*cos(x(4,end))];
p45=p24+[l5*sin(x(5,end));-l5*cos(x(5,end))];
% xi = [-deg2rad(1);-deg2rad(29);deg2rad(0);-deg2rad(9);-deg2rad(28)];
% xf = [-deg2rad(28);-deg2rad(9);deg2rad(12);-deg2rad(29);-deg2rad(1)];
xi = init_joint_angles;
% xi_dot=[-0.3954;0.2865;0.1712;-0.2940;0.3922];
xf = [xi([5;4;3;2;1])];
% xf_dot = [xi_dot([5;4;3;2;1])];
ceq = [...r45(2);p45(2); ...
     ...
    x(1:10,1) - [xi;((xf-xi)/N)];
    x(1:10,end) - [xf;impact(x(1:10,end),x(1:10,end-1))]] ;...
%     (p45(1)) - abs(r45(1))];

for i = 2:1:N-2
    e01=[-l1*sin(x(1,i));l1*cos(x(1,i))];
    e12=e01+[-l2*sin(x(2,i));l2*cos(x(2,i))];
    e23=e12+[l3*sin(x(3,i));l3*cos(x(3,i))];
    e24=e12+[l4*sin(x(4,i));-l4*cos(x(4,i))];
    e45=e24+[l5*sin(x(5,i));-l5*cos(x(5,i))];

    c = [c;-e45(2);...; ...
%         -(e45(1)) + abs(r45(1))]; ... %symm walking];
        ];
end