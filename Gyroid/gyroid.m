clc; clear
N=96;
a=2;

u = linspace(0,2*pi, N);
g1 = zeros(N,N,N);
g2 = g1;

for ix = 1:N
    for iy = 1:N
        for iz = 1:N
            x = u(ix);
            y = u(iy);
            z = u(iz);
            % note: below still works without the 2x terms but is just less
            % smoothly approaches the maximum 15
            g1(ix,iy,iz) = 10.0*cos(x/a)*sin(y/a)+10.0*cos(y/a)*sin(z/a)+10.0*cos(z/a)*sin(x/a)-.5*cos(2*x/a)*cos(2*y/a)-.5*cos(2*y/a)*cos(2*z/a)-.5*cos(2*z/a)*cos(2*x/a);
            g2(ix,iy,iz) = -10.0*cos(x/a)*sin(y/a)-10.0*cos(y/a)*sin(z/a)-10.0*cos(z/a)*sin(x/a)-.5*cos(2*x/a)*cos(2*y/a)-.5*cos(2*y/a)*cos(2*z/a)-.5*cos(2*z/a)*cos(2*x/a);
        end
    end
end

g3 = max(g1,g2);
for x=1:15
    f(x)=mean(mean(mean(g3>x)))
end
plot(1:15,f,'s-')
hold on
