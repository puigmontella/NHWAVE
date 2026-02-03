clear all
close all

fdir='results/';


load('d61');

tm = d61(:,1);
gauge2 = d61(:,2)/1000;
windowSize = 5;
gauge2 = filter(ones(1,windowSize)/windowSize,1,gauge2);
gauge3 = d61(:,3)/1000;
windowSize = 5;
gauge3 = filter(ones(1,windowSize)/windowSize,1,gauge3);
gauge4 = d61(:,4)/1000;
windowSize = 5;
gauge4 = filter(ones(1,windowSize)/windowSize,1,gauge4);


e1=load([fdir 'probe_0001']);
e2=load([fdir 'probe_0002']);
e3=load([fdir 'probe_0003']);

t0 = 1.70/1.12;
b = 0.395;

subplot(311)
plot((e1(:,1)+0.05)/t0,e1(:,2)/b,'b-')
hold on
plot(tm/t0,gauge2/b,'b--')
axis([0 2.5 -0.06 0.07])
ylabel('\eta/b')

subplot(312)
plot((e2(:,1)+0.1)/t0,e2(:,2)/b,'b-')
hold on
% load('gauge2');
plot(tm/t0,gauge3/b,'b--')
axis([0 2.5 -0.20 0.10])
ylabel('\eta/b')

subplot(313)
plot((e3(:,1)+0.04)/t0,e3(:,2)/b,'b-')
hold on
plot(tm/t0,gauge4/b,'b--')
% grid
axis([0 2.5 -0.06 0.07])
ylabel('\eta/b')
xlabel('t/t0')






Ha = load([fdir 'depth_00002']);   % size: N x M
Ha2 = load([fdir 'depth_00022']);   % size: N x M
Ha3 = load([fdir 'depth_00024']);   % size: N x M
[N, M] = size(Ha);

dx = 0.02;   % use your DX
dy = 0.02;   % use your DY

x = ((1:M) - 0.5) * dx;   % cell-center x
y = ((1:N) - 0.5) * dy;   % cell-center y
[X, Y] = meshgrid(x, y);




figure
subplot(311)
imagesc(x, y, Ha); axis xy equal tight
colorbar
title('Ha (two-layer slide thickness) at step 00001')
xlabel('x'); ylabel('y');


subplot(312)
imagesc(x, y, Ha2); axis xy equal tight
colorbar
title('Ha (two-layer slide thickness) at step 00001')
xlabel('x'); ylabel('y');

subplot(313)
imagesc(x, y, Ha3); axis xy equal tight
colorbar
title('Ha (two-layer slide thickness) at step 00001')
xlabel('x'); ylabel('y');



row = 45;

figure
plot(x, Ha(row,:), 'b-',  x, Ha2(row,:), 'r--',  x, Ha3(row,:), 'k-.');
grid on
axis([0 10 0 1.5])
ylabel('Ha (m)')
xlabel('x (m)')
legend('Ha\_00002','Ha\_00022','Ha\_00102', 'Location','best');
title(sprintf('Ha cross-section at j=%d', row));
