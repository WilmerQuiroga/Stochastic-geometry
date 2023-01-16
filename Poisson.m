
xCarre=10;%Dimensions Carre en Kilometres
yCarre=10;%Dimensions Carre en Kilometres 
areaTotal=xCarre*yCarre;


lambda=0.4; %intensity km^-1
numbPoints=poissrnd(areaTotal*lambda);%Nombre des Points Poisson. 
xx=xCarre*(rand(numbPoints,1));
yy=yCarre*(rand(numbPoints,1));
C = [xx yy];

dist=zeros(1,numbPoints);
distm=zeros(1,numbPoints);
for i=1:1:numbPoints
   X=[C(i,1) C(i,2)];
   for j=1:1:numbPoints
      if j~=i
        Y=[C(j,1) C(j,2)];   
        dist(j)=pdist2(X,Y);
      end
   end
   distm(i)=min(dist);
end
Moyenne1=mean(distm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calcul de SIR
P=unifrnd(0,1,numbPoints,1); % Puissance suit une loi uniforme sur [0, 1]
%%%%%%%%%%%%%%%%%%%%%%% Loi exponentielle Fading.
AlphaFad = exprnd(1,101,101,numbPoints);
%%%%%%%%%%%%%%%%%%%%%%%%%%% Shadowing Loi normal-log  mu=0 (mean of logarithmic values) and sigma=2 (standard deviation of logarithmic values)
AlphaShad  = lognrnd(0,2,101,101,numbPoints);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Calcul de SIR

S=zeros(101,101,numbPoints);
SIR=zeros(101,101);
dist2=zeros(1,101);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calcul de la puissance reçue S(i,j) 
for i=1:1:101       
    for j=1:1:101
           y=[i j];
         for s=1:1:numbPoints  %Nombre des Station de base.
            J=[C(s,1) C(s,2)]*10; %Position de Station de base j
            dist2(j)=pdist2(y,J);
           if dist2(j)> 0.05
             S(i,j,s)=(P(s)/(dist2(j)^2.5))*AlphaFad(i,j,s)*AlphaShad(i,j,s);
           end 
         end 
    end    
end     
dist3=zeros(1,numbPoints);
PuissR=zeros(1,numbPoints);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calcul de SIR(i)
for i=1:1:101       
       for j=1:1:101
             Summ=0;
             y=[i j];
           for s=1:1:numbPoints
             %J=[C(s,1) C(s,2)]*10;
             %dist3(s)=pdist2(y,J);
             PuissR(s)=S(i,j,s);
             Summ=Summ+S(i,j,s);
           end
            %[a,b]=min(dist3);
            [a,b]=max(PuissR);
            SIR(i,j)=10*log10(S(i,j,b)./(Summ-S(i,j,b)));
       end    
end   

%[vx,vy] = voronoi(yy,xx);
xv = [0 0 1.25 1.209 0];
yv = [0 2.89 2.15  0 0];

xv1 = [2.65 3.618 3.55 2.56 2.65];
yv1 = [9.37 8.92 7.62 8.28 9.37];

xv2 = [7.20 8.00 7.61 7.068 7.206];
yv2 = [6.07 6.13 5.26 5.74 6.07];

xv3 = [5.87 6.74 8.03 6.44 5.87];
yv3 = [8.52 6.79 7.81 9.04 8.52];


xvm=zeros(6,5);
for i=1:1:4
 for j=1:1:5   
   if i==1  
     xvm(i,j)= xv(j); 
   elseif i==2
     xvm(i,j)= xv1(j);
   elseif i==3
    xvm(i,j)= xv2(j);   
   else
     xvm(i,j)= xv3(j);
   end
 end
end
yvm=zeros(6,5);
for i=1:1:4
 for j=1:1:5   
   if i==1  
     yvm(i,j)= yv(j); 
   elseif i==2
     yvm(i,j)= yv1(j);
   elseif i==4
      yvm(i,j)= yv2(j); 
   else
     yvm(i,j)= yv3(j);
   end
 end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xxv = [3.39 3.40 4.15 4.88  4.58 3.39];
yyv = [2.75  4.08 4.17 2.77 2.55 2.75];

xxv1 = [4.12 4.88 5.58 5.96 4.54 4.12];
yyv1 = [1.75 2.77 2.66 1.17 0  1.75];

xxv2 = [4.88 4.15 4.38 5.80 5.66 4.88];
yyv2 = [2.77 4.17 4.61 4.75 2.72 2.77];

xxv3 = [2.97 4.31 4.38 4.15 3.30 2.97];
yyv3 = [5.59 4.95 4.61 4.17 4.09 5.59];


xxvm=zeros(6,6);
for i=1:1:4
 for j=1:1:6   
   if i==1  
     xxvm(i,j)= xxv(j); 
   elseif i==2
     xxvm(i,j)= xxv1(j);
   elseif i==3
    xxvm(i,j)= xxv2(j);   
   else
    xxvm(i,j)= xxv3(j);
   end
 end
end
yyvm=zeros(6,6);
for i=1:1:4
 for j=1:1:5   
   if i==1  
     yyvm(i,j)= yyv(j); 
   elseif i==2
     yyvm(i,j)= yyv1(j);
   elseif i==4
      yyvm(i,j)= yyv2(j); 
   else
     yyvm(i,j)= yyv3(j);
   end
 end
end


NumUCellule=zeros(1,100); % Nombre d'Utilisateurs par Cellule
Count=0;
for i=1:1:200
lambda=5; %intensity km^-1
numbPointsU=poissrnd(areaTotal*lambda);%Nombre des Points Poisson Utilisateurs. 
xxu=xCarre*(rand(numbPointsU,1));
yyu=yCarre*(rand(numbPointsU,1));
U = [xxu yyu];

for j=1:1:4
   Count=Count+1;
  [in,on] = inpolygon(xxu,yyu,xvm(j,:),yvm(j,:));
  NumUCellule(Count)=numel(xxu(in));
end

for s=1:1:4
   Count=Count+1;
  [in,on] = inpolygon(xxu,yyu,xxvm(s,:),yyvm(s,:));
  NumUCellule(Count)=numel(xxu(in)); 
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda=5; %intensity km^-1
numbPointsU=poissrnd(areaTotal*lambda);%Nombre des Points Poisson Utilisateurs. 
xxu=xCarre*(rand(numbPointsU,1));
yyu=yCarre*(rand(numbPointsU,1));
U = [xxu yyu];

%Calcul de SIR
P=unifrnd(0,1,numbPoints,1); % Puissance suit une loi uniforme sur [0, 1]
%%%%%%%%%%%%%%%%%%%%%%% Loi exponentielle Fading.
AlphaFad = exprnd(1,numbPointsU,numbPointsU,numbPoints);
%%%%%%%%%%%%%%%%%%%%%%%%%%% Shadowing Loi normal-log  mu=0 (mean of logarithmic values) and sigma=2 (standard deviation of logarithmic values)
AlphaShad  = lognrnd(0,2,numbPointsU,numbPointsU,numbPoints);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Calcul de SIR
S=zeros(numbPointsU,numbPointsU,numbPoints);
SIR=zeros(numbPointsU,numbPointsU);
dist2=zeros(1,numbPointsU);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calcul de la puissance reçue S(i,j) 
     
for j=1:1:numbPointsU  %Nombre des Utilisateurs.
          y=[U(j,1) U(j,2)]; % Position de Utilisateur
         for s=1:1:numbPoints  %Nombre des Station de base.
            J=[C(s,1) C(s,2)]; %Position de Station de base j
            dist2(j)=pdist2(y,J);
           if dist2(j)> 0.05
             S(U(j,1),U(j,2),s)=(P(s)/(dist2(j)^2.5))*AlphaFad(U(j,1),U(j,2),s)*AlphaShad(U(j,1),U(j,2),s);
           end 
         end 
end       
dist3=zeros(1,numbPoints);
PuissR=zeros(1,numbPoints);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calcul de SIR(i)
for j=1:1:numbPointsU  %Nombre des Utilisateurs.     
             Summ=0;
             y=[U(j,1) U(j,2)]; % Position de Utilisateur
           for s=1:1:numbPoints
             J=[C(s,1) C(s,2)]*10;
             dist3(s)=pdist2(y,J);
             %PuissR(s)=S(i,j,s);
             Summ=Summ+S(i,j,s);
           end
            %[a,b]=min(dist3);
            [a,b]=max(PuissR);
            SIR(i,j)=10*log10(S(U(j,1),U(j,2),b)./(Summ-S(U(j,1),U(j,2),b)));   
end   

%{
figure (2)
x=0:0.1:10;
y=0:0.1:10;
[X,Y]=meshgrid(x,y);
Z=SIR;
surf(X,Y,Z); %Map SIR
hold on
scatter(yy,xx,'filled','red') % Graphique de la position des stations de base 
xlim([0 10])
ylim([0 10])
voronoi(yy,xx,'white') %Graphique des cellules formées par un pavage de Voronoï
xlim([0 10])
ylim([0 10])
hold off
colormap default;
colorbar;
view(2); 
%}

figure (3)
%Plotting
scatter(yy,xx,'filled','red') % Graphique de la position des stations de base 
xlim([0 10])
ylim([0 10])
hold on
scatter(yyu,xxu,'c','blue') % Graphique de la position des stations de base 
xlim([0 10])
ylim([0 10])
scatter(vx,vy)
xlim([0 10])
ylim([0 10])
plot(xvm(2,:),yvm(2,:)) % polygon
voronoi(yy,xx) %Graphique des cellules formées par un pavage de Voronoï
xlim([0 10])
ylim([0 10])
hold off
xlabel('x');ylabel('y');

figure (4)
histogram(NumUCellule);
title('Histogramme du nombre d’utilisateurs par cellule')
