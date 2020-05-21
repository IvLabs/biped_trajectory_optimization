function xdot = dynamics(x, tau)
    xdot = zeros(4,1);
    l1 = 1;
    l2 = 1;
    m1 = 1;
    m2 = 1;
    g = 9.8;
    lc1 = l1/2;
    lc2 = l2/2;
  
    if abs(sin(x(1)+x(2)))>0.99998
        t = x(1);
        x(1) = x(2);
        x(2) = pi-x(2)+t;
        t = x(3);
        x(3) = x(4);
        x(4) = -x(4) + t;
        display(x);
        %x(5,k+1:end) = -x(5,k+1:end);
    end
   
    I1 = m1*(l1^2)/12;
    I2 = m1*(l2^2)/12;
    H11 = m1*(lc1^2) + I1 + m2*(l1^2 + lc2^2 + 2*l1*lc2*cos(x(1)));
    H22 = m2*(lc2^2) + I2;
    H12 = m2*(lc2^2 + l1*lc2*cos(x(2))) + I2;
    h = m2*l1*lc2*sin(x(2));
    G1 = m1*lc1*g*cos(x(1)) + m2*g*(lc2*cos(x(1)+x(2)) + l1*cos(x(1)));
    G2 = m2*g*lc2*cos(x(1) + x(2));    
    xdot(1) = x(3);
    xdot(2) = x(4);
    xdot(4) = (tau + ((h*H11)/H12)*(x(3)^2) + G2*H11 + 2*h*x(3)*x(4) - G1)/((-H11*H22)/H12 + H12);
    xdot(3) = -(H22*xdot(4) + h*(x(3)^2) + G2)/H12;