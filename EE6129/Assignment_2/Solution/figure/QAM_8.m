% Create 32-QAM modulator
clear all;
M = 4;
msg = randi([0,M-1],1,80); % 八进制，80个符号
% figure(1);stem(msg);
msg1 = pskmod(msg,M,pi/4); % psk调制
scatterplot(msg1,1,0,'r*'); % 画星座图
title("QPSK");
grid on;
hold on;
rectangle('Position',[-1, -1, 2, 2],'Curvature',[1, 1]);axis equal; % 画圆