function result = callme(x)
    N = 50;
    T = 1;
    h = T/N;
    result = 0;
    for k = 1:N
        result = result + ((h/2)*(x(5,k+1)^2 + x(5,k)^2));   
    end