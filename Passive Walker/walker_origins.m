function origins = walker_origins(q1,q2) 
  mH = 10; 
  m1 = 5;  
  m2 = 5;  
  a = 0.5; 
  b = 0.5; 
  l = a + b;
  g = 9.81; 


o1x = -sin(q1) * b;
o1y = cos(q1) * b; 
ohx = -sin(q1) * a - sin(q1) * b; 
ohy = cos(q1) * a + cos(q1) * b; 
o2x = -l * sin(q1) + a * sin(q2); 
o2y = l * cos(q1) - a * cos(q2); 
oex = sin(q2) * b - l * sin(q1) + a * sin(q2); 
oey = -cos(q2) * b+l* cos(q1) - a * cos(q2); 
o1 = [o1x, o1y]; 
oh = [ohx, ohy]; 
o2 = [o2x, o2y]; 
oe = [oex, oey];
origins = [o1, oh, o2, oe]; 
end