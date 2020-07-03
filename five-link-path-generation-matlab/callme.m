function result = callme(x)
    global N timer T;
    
    result = 0;
    h = T/N;
    
    for k = 1:1:N
        result = result + ((h/2)*(x(11,k+1)^2 + x(11,k)^2))+((h/2)*(x(12,k+1)^2 + x(12,k)^2))+((h/2)*(x(13,k+1)^2 + x(13,k)^2))+((h/2)*(x(14,k+1)^2 + x(14,k)^2))+((h/2)*(x(15,k+1)^2 + x(15,k)^2));
        %result = result + ((h/2)*(x(12,k+1)^2 + x(12,k)^2));
        %result = result + ((h/2)*(x(13,k+1)^2 + x(13,k)^2));
        %result = result + ((h/2)*(x(14,k+1)^2 + x(14,k)^2));
        %result = result + ((h/2)*(x(15,k+1)^2 + x(15,k)^2));

    end
    
%     timer = timer + (toc/1000);
end