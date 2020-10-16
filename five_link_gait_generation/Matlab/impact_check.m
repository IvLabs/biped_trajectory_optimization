function [xf,stop] = impact_check(x)
    global N
    for k = 1:1:N
        state = [(x(1:5,k));x(6:10,k)];
        next_state = [(x(1:5, k+1));x(6:10,k+1)];
        [xf,flag] = impact(next_state,state);
        if flag == true
            
            stop = flag;
            break;
        else
            stop = false;
        end
    end
    
end