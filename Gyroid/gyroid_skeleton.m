figure, 
syms x y z;
a = 2;

h=10.0*cos(x/a)*sin(y/a)+10.0*cos(y/a)*sin(z/a)+10.0*cos(z/a)*sin(x/a)-.5*cos(2*x/a)*cos(2*y/a)-.5*cos(2*y/a)*cos(2*z/a)-.5*cos(2*z/a)*cos(2*x/a); 
hold on 
ezimplot3(h-10.23,[-2*pi 2*pi])

x =  - x/a;
y =  - y/a;
z =  - z/a;
h2= 10.0*cos(x)*sin(y)+10.0*cos(y)*sin(z)+10.0*cos(z)*sin(x)-.5*cos(2*x)*cos(2*y)-.5*cos(2*y)*cos(2*z)-.5*cos(2*z)*cos(2*x);
ezimplot3(h2-10.23,'blue')