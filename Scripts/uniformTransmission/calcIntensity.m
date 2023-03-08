L = 250;
D = 50.8;
x = 0:0.01:D;

recpAng = rad2deg(atan(x/L)+atan((D-x)/L));

figure();
plot(x, recpAng)