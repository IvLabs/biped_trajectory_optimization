function [value, isterminal, direction] = impact_event(t,x)
    global psi;
    
    q1 = x(1);
    q2 = x(2);
    
    Dq2 = x(4);
    
    value = cos(q1+psi) - cos(q2+psi);
    %if abs(value1)<=0.009
        %value=0;
    %else
        %value=value1;
    %end
    %value1=round(value, 2);
    direction = -1;
    
    if Dq2 <= 0
        isterminal = 1;
    else
        isterminal = 0;
    end
end