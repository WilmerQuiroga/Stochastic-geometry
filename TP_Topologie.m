%%%% Vecteurs Lambda,B0,B1,Xc pour Tracer B0,B1 et Xc en fonction de Lambda
VecLambda=linspace(0,200,201);
VecB0=zeros(1,length(VecLambda));
VecB1=zeros(1,length(VecLambda));
VecXc=zeros(1,length(VecLambda));
%%%%%%%%%%%%%%%%%%%%%Moyenne
N=5; %Nombre Points pour faire la moyenne 
VecB0Moy=zeros(1,N);
VecB1Moy=zeros(1,N);
VecXcMoy=zeros(1,N);
%%%%%%%%%%%%%%% moyenne sur 100 simulations par valeur de Lambda
for i=1:1:length(VecLambda)
    for j=1:1:length(VecB0Moy)
        [VecB0Moy(j),VecB1Moy(j),VecXcMoy(j)]=SimplexCalcul(VecLambda(i));
    end
    VecB0(i)=mean(VecB0Moy);
    VecB1(i)=mean(VecB1Moy);
    VecXc(i)=mean(VecXcMoy);
end

figure(2)
plot(VecLambda,VecB0);
hold on
plot(VecLambda,VecB1);
plot(VecLambda,VecXc);
hold off
title('B0,B1 et Xc en fonction de Lambda')
legend({'B0 vs Lamdda','B1 vs Lamdda','Xc vs Lamdda'},'Location','southwest')


function [B0,B1,Xc]=SimplexCalcul(Lambda)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Poisson Processus
xCarre=1;%Dimensions Carre en Kilometres
yCarre=1;%Dimensions Carre en Kilometres 
areaTotal=xCarre*yCarre;
lambda=Lambda; %intensity km^-1
numbPoints=poissrnd(areaTotal*lambda);%Nombre des Points Poisson. 
xx=xCarre*(rand(numbPoints,1));
yy=yCarre*(rand(numbPoints,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1-Simplex
Simplex1x=zeros(1,2);
Simplex1y=zeros(1,2);
count=1;
for i=1:1:numbPoints
   X=[xx(i) yy(i)];
   for j=i+1:1:numbPoints
       if (j>i)
           Y=[xx(j) yy(j)];
           distance=pdist2(X,Y);
           if(distance<0.1)
               Simplex1x(count,1)=xx(i);
               Simplex1y(count,1)=yy(i);
               Simplex1x(count,2)=xx(j);
               Simplex1y(count,2)=yy(j);
               count=count+1;
           end
       end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2-Simplex
Simplex2x=zeros(1,3);
Simplex2y=zeros(1,3);
count2=1;
Simplex1xLen=count-1;
for i=1:1:Simplex1xLen
  X=[Simplex1x(i,1) Simplex1y(i,1)];
  Y=[Simplex1x(i,2) Simplex1y(i,2)];
  for j=1:1:numbPoints
         Z=[xx(j) yy(j)];
         distance2=pdist2(X,Z);
         distance3=pdist2(Y,Z);
         
         if (X(1)~=Z(1)&&X(2)~=Z(2))&&(Y(1)~=Z(1)&&Y(2)~=Z(2))&&(isempty(find(Simplex2x==Z(1))))
           if distance2<0.1 && distance3<0.1
             Simplex2x(count2,1)=Simplex1x(i,1);
             Simplex2y(count2,1)=Simplex1y(i,1);
             Simplex2x(count2,2)=Simplex1x(i,2);
             Simplex2y(count2,2)=Simplex1y(i,2);
             Simplex2x(count2,3)=xx(j);
             Simplex2y(count2,3)=yy(j);
             Simplex2x(count2,4)=Simplex1x(i,1);
             Simplex2y(count2,4)=Simplex1y(i,1);
             count2=count2+1;
           end
         end
  end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calcul Matrux d1
Simplex2xLen=count2-1;
d1=zeros(numbPoints,Simplex1xLen);
for i=1:1:Simplex1xLen
   Px1=[Simplex1x(i,1) Simplex1y(i,1)];
   Px2=[Simplex1x(i,2) Simplex1y(i,2)];
   for j=1:1:numbPoints
       X=[xx(j,1) yy(j,1)];
       if Px1(1)==X(1) &&  Px1(2)==X(2)
          d1(j,i)=-1;
       end 
       if Px2(1)==X(1) &&  Px2(2)==X(2)
          d1(j,i)=1;
       end 
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calcul Matrux d2
d2=zeros(Simplex1xLen,Simplex2xLen);
for i=1:1:Simplex2xLen
   Px1=[Simplex2x(i,1) Simplex2y(i,1)];
   Px2=[Simplex2x(i,2) Simplex2y(i,2)];
   Px3=[Simplex2x(i,3) Simplex2y(i,3)];
   
   for j=1:1:Simplex1xLen
       X1=[Simplex1x(j,1) Simplex1y(j,1)];
       X2=[Simplex1x(j,2) Simplex1y(j,2)];
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Px1->Px2
       if Px1(1)==X1(1) && Px1(2)==X1(2) && Px2(1)==X2(1) && Px2(2)==X2(2)
          d2(j,i)=1;
       end 
       if Px1(1)==X2(1) && Px1(2)==X2(2) && Px2(1)==X1(1) && Px2(2)==X1(2)
          d2(j,i)=-1;
       end 
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Px2->Px3
       if Px2(1)==X1(1) &&  Px2(2)==X1(2) && Px3(1)==X2(1) &&  Px3(2)==X2(2)
          d2(j,i)=1;
       end 
       if Px2(1)==X2(1) &&  Px2(2)==X2(2) && Px3(1)==X1(1) &&  Px3(2)==X1(2)
          d2(j,i)=-1;
       end 
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Px3->Px1
       if Px3(1)==X1(1) &&  Px3(2)==X1(2) && Px1(1)==X2(1) &&  Px1(2)==X2(2)
          d2(j,i)=1;
       end 
       if Px3(1)==X2(1) &&  Px3(2)==X2(2) && Px1(1)==X1(1) &&  Px1(2)==X1(2)
          d2(j,i)=-1;
       end 
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% B0,B1,Xc
B0=numbPoints-rank(d1);
[f, g]=size(null(d1));
B1=g-rank(d2);
Xc=numbPoints-Simplex1xLen+Simplex2xLen;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot
%{
figure(1)
scatter(xx,yy)
hold on
%for i=1:1:length(Simplex1x)
%   plot(Simplex1x(i,:),Simplex1y(i,:)); 
%end
for i=1:1:Simplex2xLen
  plot(Simplex2x(i,:),Simplex2y(i,:)); 
end
title('2-Simplexes')
hold off
%}
end