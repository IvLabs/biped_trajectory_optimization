    function [xn,flag] = impact(x,xp)

global N m1 m2 m3 m4 m5 l1 l2 l4 l3 l5 timer T
%    xp = xp.';
r01=[l1*sin(x(1));l1*cos(x(1))];
r12=r01+[l2*sin(x(2));l2*cos(x(2))];
r23=r12+[l3*sin(x(3));l3*cos(x(3))];
r24=r12+[l4*sin(x(4));-l4*cos(x(4))];
r45=r24+[l5*sin(x(5));-l5*cos(x(5))];

p01=[l1*sin(xp(1));l1*cos(xp(1))];
p12=p01+[l2*sin(xp(2));l2*cos(xp(2))];
p23=p12+[l3*sin(xp(3));l3*cos(xp(3))];
p24=p12+[l4*sin(xp(4));-l4*cos(xp(4))];
p45=p24+[l5*sin(xp(5));-l5*cos(xp(5))];

m = [m1 m2 m3 m4 m5];
l = [l1 l2 l3 l4 l5];

        %compute jacobian
        J = zeros(2,5);
        for i = 1:1:2
            for j = 1:1:5
                    if x(j) ~= xp(j)
                        J(i,j) = (r45(i)-p45(i))/(x(j) - xp(j));
                    else
                        J(i,j) = (r45(i)-p45(i))/(x(j) - xp(j) + 0.05);
                    end
                
            end
        end
        Jt = J.';
%         Impulses = zeros(7,1);
        % ImpaVel = zeros(5,1);

        %compute W
        W = zeros(5,5);
        k = zeros(5,5);
        q = zeros(5,5);
        for i = 1:1:5
            for j = 1:1:5
                if i == j
                    if i == 3
                       k(i,j) = (m(i)*(l(i)^2)/3); 
                    else
                       k(i,j) = -(m(i)*(l(i)^2)/6);
                    end
                elseif j>i
                    k(i,j) = 0;
                else
                    if (i == 2 || i == 5) && (j == i-1)
                        k(i,j) = (-1/2)*(m(i)+m(j))*l(i)*l(j);
                    elseif (i == 3) && (j < i)
                        k(i,j) = m(i)*l(i)*l(j)/2;
                    elseif (i == 4 || i == 5) && (j == 3)
                        k(i,j) = -m(j)*l(i)*l(j)/2;
                    else
                        const = 0;     
                        for p = (j+1):1:(i-1)
                            const = const + (m(p)*l(j)*l(i));
                        end
                        k(i,j) = (-m(j)*l(i)*l(j)/2)- (m(i)*l(i)*l(j)/2) - const;
                    end

                    if (i < 4 && j < 4) || (i > 3 && j > 3)
                        q(i,j) = x(i) - x(j);
                    else
                        q(i,j) = x(i) + x(j);
                    end
                end
                W(i,j) = k(i,j)*(cos(q(i,j)));
            end    
        end

        Impulses = [(W*(x(6:10)));0;0]\[W,Jt;J,zeros(2,2)];
        t = [0,0,0,0,-1; ...
             0,0,0,-1,0; ...
             0,0,1,0,0; ...
             0,-1,0,0,0; ...
             -1,0,0,0,0;];
         
         
         
%         xn = [T,zeros(5,5);zeros(5,5),T]*[x(1:5);Impulses(1:5).'];

        xn = t*Impulses(1:5).';
              

end