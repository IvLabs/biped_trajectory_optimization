function [c,ceq] = colons(x)
    global N timer;
    [c,ceq] = boundry(x);
    h = 1/N;
    
        
%     terrian is flat
    for k = 1:1:N
        %[xt,flag]=impact_check(x); 
        %if (flag==true)
        %    x=xt;
        %end
        state = [(x(1:5, k));x(6:10,k)];
        next_state = [(x(1:5, k+1));x(6:10,k+1)];
        u = x(11:15, k);
        u_next = x(11:15, k+1);
        x_dot = dynamics5(state,u);
        x_dot_next = dynamics5(next_state,u_next);
        ceq = [ceq;next_state - state - ((h/2)*(x_dot+x_dot_next))];
                  
     end
%         ceq = [ceq , next_state - state - ((h/2)*(x_dot+x_dot_next))];
    end
%     ceq = [ceq ,[wrapToPi(x(1:5,1));x(6:10,1)]- ...
%         [-deg2rad(1.18);-deg2rad(31.29);deg2rad(0);-deg2rad(10.19);-deg2rad(28.14); ...
%                      1.03;-0.16;0;-1.09;0.1]...
%         ,[wrapToPi(x(1:5,N));x(6:10,N)]-[-deg2rad(28.14);-deg2rad(10.19);-deg2rad(30);-deg2rad(31.29);-deg2rad(1.18); ...
%                      0.1;-1.09;0;-0.16;1.03] ...
%                      ];
%     timer = timer + (toc/1000);  
% fprintf('%i seconds \n',timer);
