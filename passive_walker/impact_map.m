function [x_plus] = impact_map(x_minus)
    global p6 p7 p8 p9 p10 p11;
    
    q1 = x_minus(1);
    q2 = x_minus(2);
    Dq1 = x_minus(3);
    Dq2 = x_minus(4);
    
    Qp = [-cos(q1-q2) * p7 + p6 p8- cos(q1-q2)*p7;...
        -cos(q1 - q2)*p7 p8;];
    Qn = [cos(q1-q2) * p10 - p9 -p11; -p9, 0];
    
    Dq_plus = Qp\Qn*[Dq1; Dq2];
    
    q_plus = [q2; q1];
    x_plus = [q_plus; Dq_plus];
end
    