PBs=sfr_lte; %Data Operateur
xx=PBs(:,1);
yy=PBs(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calcul Lambda
nombrePoints=length(xx);
count1=0;
for i=1:1:nombrePoints
    if(PBs(i,1)>-5000 && PBs(i,1)<5000)
        if (PBs(i,2)>-2000 && PBs(i,2)<2000)
            count1=count1+1;
        end
    end 
end
Lambdam=count1/(4000*10000);%m^-2
Lambdakm=count1/(10*4) %km^-2

figure(1)
scatter(xx,yy);
rectangle('Position',[-5000 -2000 5000*2 2000*2]);
O=[0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PBs2=zeros(count1,2);
NombrePointsExte=nombrePoints-count1;
PBs3=zeros(NombrePointsExte,2);
t=1;
for i=1:1:nombrePoints
    if(PBs(i,1)>-5000 && PBs(i,1)<5000)
        if (PBs(i,2)>-2000 && PBs(i,2)<2000)
            PBs2(t,1)=PBs(i,1);
            PBs2(t,2)=PBs(i,2);
            t=t+1;
        end
    end 
end
Distances1=zeros(1,nombrePoints);
Distances2=zeros(1,length(PBs2));
for i=1:1:length(PBs2)
   X=[PBs2(i,1) PBs2(i,2)];
   for j=1:1:nombrePoints
       Y=[PBs(j,1) PBs(j,2)];
       if(X(1)~=Y(1) && X(2)~=Y(2) )
        Distances1(j)=pdist2(X,Y);
       end
   end
   Distances2(i)=min(Distances1);
end
count=0;
R1=0:0.1:300;
F=0;
for i=1:1:length(R1)
   for j=1:1:length(PBs2)
     if(Distances2(j)<R1(i))
       count=count+1;
     end 
   end
   F(i)=count/(length(PBs2));
   count=0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=1;
for i=1:1:nombrePoints
    if(PBs(i,1)<=-5000||PBs(i,1)>=5000||PBs(i,2)<=-2000 || PBs(i,2)>=2000)
            PBs3(t,1)=PBs(i,1);
            PBs3(t,2)=PBs(i,2);
            t=t+1;
    end
end

Distances3=zeros(1,nombrePoints);
Distances4=zeros(1,NombrePointsExte);

for i=1:1:NombrePointsExte
   X=[PBs3(i,1) PBs3(i,2)];
   for j=1:1:nombrePoints
       Y=[PBs(j,1) PBs(j,2)];
       if(X(1)~=Y(1) && X(2)~=Y(2))
        Distances3(j)=pdist2(X,Y);
       end
   end
   Distances4(i)=min(Distances3);
end
count=0;
G=0;
for i=1:1:length(R1)
   for j=1:NombrePointsExte
     if(Distances4(j)<R1(i))
       count=count+1;
     end 
   end
   G(i)=count/(length(PBs3));
   count=0;
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Functions F et G

%R4=0:0.1:300;
%Beta=0.6;
%Count5=1;
%for i=1:1:length(R4)
%   J1(i)=1/(1-Beta +Beta*exp(-((Lambda*pi*R4(i)^2)/Beta)));
%end
%J2=1./(1-Beta+Beta*exp(-((Lambdakm*pi*(0.001*R4).^2)/Beta)));


J=(1-G)./(1-F);
R4=0:0.1:300;
Betam=0.01:0.01:1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Méthode MMSE
err=zeros(1,length(Betam));
for i=1:1:length(Betam)
  J1=1./(1-Betam(i) +Betam(i)*exp(-((Lambdakm*pi*(0.001*R4).^2)/Betam(i))));
  err(i)=immse(J,J1);
end
[a,b]=min(err);
BetaG=Betam(b)
JF=1./(1-BetaG+BetaG*exp(-((Lambdakm*pi*(0.001*R4).^2)/BetaG)));

figure(3)
plot(R4,JF);
hold on
plot(R1,J);
hold off
xlabel('Operateur: SFR Lte');ylabel('J(R)')
title('Functions Je(R) et Jt(R)')
legend({'J Théorique','J expérimental'},'Location','northwest')
%}

%{
Distances=zeros(1,nombrePoints);
Distances2=zeros(1,nombrePoints);

for i=1:1:nombrePoints
   X=[PBs(i,1) PBs(i,2)];
   for j=1:1:nombrePoints
       Y=[PBs(j,1) PBs(j,2)];
       if(i~=j)
        Distances(j)=pdist2(X,Y);
       end
   end
   Distances2(i)=min(Distances);
end

count=0;
G=zeros(1,700);
Count5=1;
for R=0:0.1:700
   for i=1:1:nombrePoints
     if(Distances2(i)<R)
       count=count+1;
     end 
   end
   G(Count5)=count/nombrePoints;
   count=0;
   Count5=Count5+1;
end

figure(1)
plot(R1,G)
xlabel('R');ylabel('G(R)')
title('Function G(R)')
%}

%Calcul de la funtion F(r)
%t=1;
%distance=zeros(count1,count1);
%Count3=0;
%Radius=1:1:700;
%F=zeros(1,700);
%{
for R=1:1:700
  for i=1:count1
     X=[PBs2(i,1) PBs2(i,2)];
    for j=1:count1
        if i~=j
          Y=[PBs2(j,1) PBs2(j,2)];
          distance=pdist2(X,Y);
          if distance<R
              Count3=Count3+1;
              break
          end
        end 
    end
  end
  F(R)=Count3/count1;
  Count3=0;
end
%}

%{
Count3=0;
for R=1:1:700
  for i=1:nombrePoints
     X=[PBs(i,1) PBs(i,2)];
    for j=1:nombrePoints
        if i~=j
          Y=[PBs(j,1) PBs(j,2)];
          distance=pdist2(X,Y);
          if distance<R
              Count3=Count3+1;
              break
          end
        end 
    end
  end
  F(R)=Count3/nombrePoints;
  Count3=0;
end

%{
figure(2)
plot(Radius,F)
xlabel('R');ylabel('F(R)')
title('Function F(R)')

R=1:1:700;
Beta=0.6;
Lambda=9.3;
J=1./(1-Beta +Beta*exp(-(Lambda*pi*R.^2/Beta)));
plot(R,J);

%}

%{

% Calcul de la funtion G
for R=1:1:700
  for i=1:1:nombrePoints
     X=[PBs(i,1) PBs(i,2)];
     for j=1:1:nombrePoints
       if(PBs(i,1)>-6000 && PBs(i,1)<6000 && j~=i)
            if (PBs(i,2)>-2500 && PBs(i,2)<2500)
               Y=[PBs(j,1) PBs(j,2)];   
               distance=pdist2(X,Y);
               if distance<R
                   count3=count3+1;
               end
            end 
        end 
     end
  end
  f(R)=count3/count1;
  count3=0;
end

plot(r,F);

%}

%}
%{

dist=zeros(1,nombrePoints);
distm=zeros(1,nombrePoints);
for i=1:1:nombrePoints
   X=[PBs(i,1) PBs(i,2)];
   for j=1:1:nombrePoints
      if j~=i
        Y=[PBs(j,1) PBs(j,2)];   
        dist(j)=pdist2(X,Y);
      end
   end
   distm(i)=min(dist);
end
Moyenne1=mean(distm);
Count=0;
Voisines=zeros(1,nombrePoints);
distt=0;
for i=1:nombrePoints
    X=[PBs(i,1) PBs(i,2)];
    for j=1:1:nombrePoints
      if (j~=i)
        Y=[PBs(j,1) PBs(j,2)];   
        distt=pdist2(X,Y);
      end
      if distt<Moyenne1
          Count=Count+1;
      end
      
    end
   Voisines(i)=Count;
   Count=0;
end
count2=0;
for i=1:nombrePoints
    if(Voisines(i)==3)
     count2=count2+1;  
    end
end 
Genibre=zeros(count2,2);
Count3=1;
for i=1:nombrePoints
    if(Voisines(i)==3 || Voisines(i)==4)
      Genibre(Count3,1)=PBs(i,1);
      Genibre(Count3,2)=PBs(i,2);
      Count3=Count3+1;
    end
end 
xx1=Genibre(:,1);
yy1=Genibre(:,2);

figure(2)
scatter(xx1,yy1);

%}



