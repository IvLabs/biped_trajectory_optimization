function [movie_storage,winsize] = animate_walker(time,...
                                        impacts,draw_interval,q,fig)
   global psi;
   
    q1 = q(:,1); 
    q2 = q(:,2);

  
 static_origins = walker_origins(q1,q2); 
 offset_origins = travel_offset(impacts, static_origins);

 
 o1x = offset_origins(:,1); 
 o1y = offset_origins(:,2); 
 ohx = offset_origins(:,3); 
 ohy = offset_origins(:,4); 
 o2x = offset_origins(:,5); 
 o2y = offset_origins(:,6);
 oex = offset_origins(:,7); 
 oey = offset_origins(:,8); 
 ocx = offset_origins(:,9); 
 ocy = offset_origins(:,10);
 
 if draw_interval == 0 
     time_indexes = 1:length(time); 
 else
     time_indexes(1) = 1; 
     index = find(time > 0,1); 
     i = 1;

     
      while ~isempty(index) 
           time_indexes(i) = index; 
            index = find(time >=i *draw_interval,1); 
            
            i = i+1; 
      end
  end
 max_x = max(max([ocx, oex, o1x, ohx, o2x])); 
 max_y = max(max([ocy, oey, o1y, ohy, o2y])); 
 min_x = min(min([ocx, oex, o1x, ohx, o2x])); 
 min_y = min(min([ocy, oey, o1y, ohy, o2y]));

 axis_scale = [min_x, max_x min_y max_y];

 winsize = get(fig, 'Position'); 
 winsize(1:2) = [0 0]; 
 numframes = length(time_indexes); 
 movie_storage = moviein(numframes,fig,winsize); 
 set(fig,'NextPlot','replacechildren')

 set(0,'defaulttextinterpreter','Tex')
 
  for j=1:numframes 
       index = time_indexes(j);

        clf;

        hold on; 
        axis(axis_scale); 
        axis manual;

         
        plot(o1x(index), o1y(index), 'or'); 
        plot(ohx(index), ohy(index), 'or'); 
        plot(o2x(index), o2y(index), 'og');

         
        line([ocx(index) ohx(index)],[ocy(index) ohy(index)],... 
            'LineWidth',3, 'Color', 'r'); 
        line([ohx(index) oex(index)],[ohy(index) oey(index)],... 
            'LineWidth',3, 'Color', 'g');

         
        line([axis_scale(1); axis_scale(2)], ... 
            [-axis_scale(1)*tan(psi) -axis_scale(2)*tan(psi)],... 
            'LineWidth',3, 'Color', 'k');
        
         
        time_instant = ['\ittime \rm= ' num2str(time(index), '%.2f') 's']; 
        text(axis_scale(2),axis_scale(4), time_instant,... 
            'VerticalAlignment','top',... 
            'HorizontalAlignment','right',... 
            'FontSize',16)

         
        model_data = {sprintf('\\itq 1 \\rm= %.2f',q(index,1)*180/pi);... 
        sprintf('\\itq 2 \\rm= %.2f', q(index,2)*180/pi)}; 
        text(axis_scale(1),axis_scale(3), model_data,... 
            'VerticalAlignment','bottom',... 
            'HorizontalAlignment','left',... 
            'FontSize',16)

        
        hold off;

        movie_storage(:,j) = getframe(fig,winsize); 
  end
end