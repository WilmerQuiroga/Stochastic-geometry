r=320; %Rayon du disk en metres 
xx0=0; yy0=0; %centre of disk
p=0.01;
K = 10^6;
w=180*10^3;
c=162*10^3;
%Simulate Poisson Processus:
Count=0;
%Ps=zeros(1,20);
Sinicial=194;
lambda=0.01*p;

%{
for s=195:1:215
 for j=1:10000
  lambda=0.01*p; %intensity m^-1
  areaTotal=pi*r^2;
  numbPoints=poissrnd(areaTotal*(lambda));
  theta=2*pi*(rand(numbPoints,1)); %angular coordinates
  rho=r*sqrt(rand(numbPoints,1)); %radial coordinates
  %Convert from polar to Cartesian coordinates
  [xx,yy]=pol2cart(theta,rho); %x/y coordinates of Poisson points
  C=[xx,yy];
  xx=xx+xx0;
  yy=yy+yy0;
  F=0;
for i=1:1:numbPoints 
    O=[0 0];
    Y=[C(i,1) C(i,2)];
    dis=pdist2(Y,O);
    m=exprnd(1);
    F=F+min(7,ceil(c/(w*log2(1+(K*m)/(dis^2.8)))));
end 
    if F>s
       Count=Count+1;
    end
 end
Ps(s-Sinicial)=Count/10000;
Count=0;
end
%}

EFm=EF*gamma((2/2.8)+1);

%{
F=0;
for i=1:1:numbPoints 
    O=[0 0];
    Y=[C(i,1) C(i,2)];
    dis=pdist2(Y,O);
    F=F+ceil(c/(w*log2(1+K/(dis^2.8))));
end 
beta=7;
alpha=beta*r*sqrt(lambda*pi);
rk=zeros(1,7);

for i=1:1:6
 rk(i)=(K/(2^(c/(w*i))-1))^(1/2.8); 
end  
rk(7)=r;
EF=0;
for i=1:1:7   
 if i==1
 EF=EF+lambda*pi*i*(rk(i)^2); 
 else
 k1=i-1;
 EF=EF+lambda*pi*i*(rk(i)^2-rk(k1)^2); 
 end
end  

PS1=zeros(1,20);
for i=0:1:20
    if(i==0)
    PS1(1)=PS(160);
    else
    PS1(i)=PS(160+i);
    end
end
%}    
beta=7;
alpha=beta*r*sqrt(lambda*pi);

y=215-EFm:1:235-EFm;
t=exp(-((y/beta)+(alpha^2/beta^2)).*log(1+((beta*y)/alpha^2))+y/beta);
figure (1)
plot(y,t);
hold on
plot(y,Ps);
hold off
xlabel('S');ylabel('Probabilite');
legend({'Equation','P(F>E[F]+y)'},'Location','northwest');

%{
th = 0:pi/50:2*pi;
xunit = r * cos(th) + xx0;
yunit = r * sin(th) + yy0;


figure(2)
%Plotting
scatter(xx,yy);
hold on
plot(xunit, yunit);
scatter(0,0,'filled','red');
hold off
xlabel('x');ylabel('y');
axis square;
%}