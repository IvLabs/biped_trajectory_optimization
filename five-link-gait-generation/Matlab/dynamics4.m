function xdot = dynamics4(x,tau)

    global m1 m2 m3 m4 m5 l1 l2 l3 l4 l5 timer
  
    m = [m1 m2 m3 m4 m5];
    l = [l1 l2 l3 l4 l5];
    n = zeros(5);
    p = zeros(5,5);
    q = zeros(5,5);
    g = zeros(5);
    D = zeros(5,5);
    H = zeros(5,5);
    G = zeros(5,1);
    f = 1:1:5;
    
    for i = 1:length(f)    
        if (i == 3) || (i == 5)
            n(i) = 0;
        else
            for b = i+1:1:5
                n(i) = m(b)*l(i);
            end
        end
        
        for j = 1:1:5
            
            if (i == j)
                p(i,j) = (m(i)*(l(i)^2)/3) + (n(i)*l(i));
            elseif (i == 3) && (j > i)
                p(i,j) = 0;
            else
                p(i,j) = m(j)*(l(j)/2)*l(i) + m(j)*l(i);
            end
            
            if (i < 4 && j < 4) || (i > 3 && j > 3)
                q(i,j) = x(i) - x(j);
            else
                q(i,j) = x(i) + x(j);
            end
            
            D(i,j) = p(i,j)*cos(q(i,j));
            H(i,j) = p(i,j)*sin(q(i,j));
        end
        
        if i < 4
            g(i) = -(m(i)*(l(i)/2) + n(i))*(9.81);
        else
            g(i) = (m(i)*(l(i)/2) + n(i))*(9.81); 
        end
        
        G(i) = g(i)*sin(x(i));
    end
    
    T = [0;tau(2);tau(3);tau(4);tau(5)];
    
    dot = (((T) - (G) - (H*(x(6:10).^2)))\D);
    xdot = [x(6:10);dot.'];
%     timer = timer + toc

end