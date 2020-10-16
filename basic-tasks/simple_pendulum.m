% damping,gravity,length and length of animation
b = -0.5; g = -9.8; l = 2; tl = 100;

%expression for angle function in time
syms r(t) 
%angular velocity
r_dot = diff(r,t,1); 
%diffrential equation governing the simple pendulum
equ = diff(r,t,2) == ((b * r_dot)+((g/l) * r)); 
%initial conditions
con1 = r(0) == pi ; con2 = r_dot(0) == 0;
conds = [con1 con2];

%solving the diffrential equation to get the expression of solution
r_sol = dsolve(equ, conds);
r_solsim = simplify(r_sol);

%time data of animation
t = 0:0.1:tl;
%getting numerical values of theta for each instant of animation
theta = double(subs(r_solsim,t));

%plot for theta vs t (uncomment for viewing)
%plot(t,theta)

axis([-4,4,-4,4])

for i = 1:1:length(theta)
    x = l*sin(theta(i));
    y = -l*cos(theta(i));
    rod = line([0 x],[0 y],'LineWidth',2,'Marker','o','MarkerFaceColor','r');
    ball = viscircles([x y],0.09);
    pause(0.1)
    delete(rod)
    delete(ball)
end