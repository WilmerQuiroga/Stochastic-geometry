r=320; %Rayon du disk en metres 
xx0=0; yy0=0; %centre of disk
%Simulate Poisson Processus:
lambda=0.01; %intensity m^-2
areaTotal=pi*r^2;
numbPoints=poissrnd(areaTotal*lambda);

theta=2*pi*(rand(numbPoints,1)); %angular coordinates
rho=r*sqrt(rand(numbPoints,1)); %radial coordinates
%Convert from polar to Cartesian coordinates
[xx,yy]=pol2cart(theta,rho); %x/y coordinates of Poisson points
C=[xx,yy];
th = 0:pi/50:2*pi;
xunit = r * cos(th) + xx0;
yunit = r * sin(th) + yy0;

figure(1)
%Plotting
scatter(xx,yy);
hold on
plot(xunit, yunit);
scatter(0,0,'filled','red');
hold off
xlabel('x');ylabel('y');
axis square;