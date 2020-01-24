function [Dx] = equation_of_motion(t,x)
    global p1 p2 p3 p4 p5;
    
    q1 = x(1);
    q2 = x(2);
    Dq1 = x(3);
    Dq2 = x(4);
    Dq = [Dq1; Dq2];
    disp(size(Dq))
    
    m11 = p1;
    m12 = -p2*cos(q1 - q2);
    m21=m12;
    m22 = p3;
    c11 = 0;
    c12 = -Dq2*sin(q1 - q2) *p2;
    c21 = -Dq1*sin(q1 - q2) *p2;
    c22 = 0;
    g1 = -sin(q1) *p4;
    g2 = sin(q2) *p5;
    
    M = [m11, m12; m21, m22];
    C = [c11, c12; c21, c22];
    G = [g1; g2];
    disp(size(C))
    
    DDq = M\(-C*Dq-G);
    Dx = [Dq; DDq];
end