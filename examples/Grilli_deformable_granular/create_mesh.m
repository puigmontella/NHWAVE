M = 250;
N = 10;
dx = 0.025;

x = ((1:M)-0.5)*dx;


h0 = 0.33;
theta = 35*pi/180;
x_s = 0.6;
Ls = 0.9;
x_e = x_s + Ls;

depth = zeros(N,M);

for i = 1:M
    if x(i) < x_s
        depth(:,i) =-0.0901+ x(i)*tan(theta);
    else
        depth(:,i) = h0;
    end
end

dlmwrite('depth.txt', depth, 'delimiter', ' ', 'precision', 6);
