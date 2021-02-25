
clear all
close all
clc
%%%%%%%%%%%%%%%%% SAR Diffusion Model for Drug Distribution%%%%%%%%%%%%%%%
%%%% Yusuf Al-Husaini, Alec Frosh, Adam Mason, Ashrak Ullah%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  ALL YOU HAVE TO DO IS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ENTER THE RESPECTIVE VALUE FOR alpha FROM TABLE 1 %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% AND THE TOTAL LENGTH OF TIME UNDER "T" %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%    THE REST IS DONE     %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%
alpha1=0.1;
%%%%%%%%%%%%
T=200;   %Total length of time, months%
%%%%%%%%%%%%

%Grid range%
a=0;
b=20;
c=0;
d=20;


n =80;  %number of discrete point in i and j%


%h=2;    

h=(b-a)/n;  %Length of discrete point divisions%



%Onto a 200x200km grid at a scale of 1:10 the towns are placed, each with..%
%..a differing a*exp^(-bx^2) and a unique initial density of addicts%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OXFORD%
xo=4.4;      %x-coordinate%
yo=9.4;     %y-coordinate %
%IC distribution of susceptible for "(ao)*Exp[-(bo)*x^2]"%
ao=3.537;   
bo=6.975;
aao=0.0001;      %Addicted density at t=0%
%This indicates 373 addicted people in Oxford initially
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CAMBRIDGE%
xc=13.8;   %x-coordinate%
yc=14.5;   %y-coordinate%
%IC distribution of susceptible for "(ac)*Exp[-(bc)*x^2]"%
ac=3.135;  
bc=7.147;
aac=0.0002;     %Addicted density at t=0%
%This indicates 373 addicted people in Cambridge initially
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%LONDON%
xl=12.2;   %x-coordinate %
yl=7;    %y-coordinate%
%IC distribution of susceptible for "(al)*Exp[-(bl)*x^2]"%
al=6.285;
bl=1.16*10^-17.12;
aal=0.001;    %Addicted density at t=0%
%This indicates 373 addicted people in London initially
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Diffusion EquationCoefficients, spatial expansion of the density, per time step%
delta1=0.01;
delta2=0.01;
delta3=0.01;



hh=(b-a)/n;
x=a:hh:b;
y=c:hh:d;
k= (((b-a)/n)^2 )    /(4*max(delta1(:)));

%alpha1=1;
%alpha2(i,j)= %1/4;
            
%SAR model parameters%
for i=1:n+1            %"For loop" applied to every cell in the i,j matrix%
    for j=1:n+1
        alpha(i,j) = alpha1*((373)/(8*10^4));%Probabilty of meeting an addict%
                %alpha1*((1/2)*exp(-0.2*((x(i)-xo).^2 + (y(j)-yo).^2  ))...
                %+(1/2)*exp(-0.2*((x(i)-xc).^2 + (y(j)-yc).^2  ))... 
                %+(1/10)*exp(-0.2*((x(i)-xl).^2 + (y(j)-yl).^2  )) )  ;     %Proportion of S joining A per time step%    
       beta(i,j) = 0.1*(1-alpha1); 
       gamma(i,j) = 0.1*(0.4+(0.2*alpha1));  %Proportion of R joining S per time step%
        
    end
end



%INITIAL CONDITIONS: Defined for each town$
for i=1:n+1
    for j=1:n+1
        S0(i,j)=ao*exp(-bo*((x(i)-xo).^64 + (y(j)-yo).^64  ))...  %Oxford
                +ac*exp(-bc*((x(i)-xc).^64 + (y(j)-yc).^64  ))... %Cambridge
                +al*exp(-bl*((x(i)-xl).^64 + (y(j)-yl).^64  ));   %London
        A0(i,j)=aao*exp(-((x(i)-xo).^64 + (y(j)-yo).^64  ))...
            +aac*exp(-((x(i)-xc).^64 + (y(j)-yc).^64  ))...
            +aal*exp(-((x(i)-xl).^64 + (y(j)-yl).^64  ));
        R0(i,j)=0;
    end
end


S=S0;
Snew=S0;

A=A0;
Anew=A0;

R=R0;
Rnew=R0;

%"for loops" for discretized Diffusion Equation adapted for the SAR model%
for m=1:T/k   %for time steps 1 to total time divided by k times steps%
    for i=2:n      %for all i from 2 until n %
        for j=2:n  %for all j from 2 until n %
            Snew(i,j)=S(i,j)+k*(-alpha(i,j)*S(i,j)*A(i,j) ...
                +delta1*((S(i-1,j)-2*S(i,j)+S(i+1,j))+S(i,j-1)-2*S(i,j)+S(i,j+1))/h^2);
            Anew(i,j)=A(i,j)+k*(alpha(i,j)*S(i,j)*A(i,j) ...
                +delta2*((A(i-1,j)-2*A(i,j)+A(i+1,j))+A(i,j-1)-2*A(i,j)+A(i,j+1))/h^2 ...
                +gamma(i,j)*R(i,j)*A(i,j));
            Rnew(i,j)=R(i,j)+k*(beta(i,j)*A(i,j)+delta3*(((R(i-1,j)-2*R(i,j)...
                +R(i+1,j))/h^2)+((R(i,j-1)-2*R(i,j)+R(i,j+1))/(h^2)))-gamma(i,j)*R(i,j)*A(i,j));
        end
    end
    
    %BOUNDARY CONDITIONS%
    for j=2:n
        
        Snew(1,j)=Snew(2,j);      %left boundary
        Snew(n+1,j)=Snew(n,j);    %right boundary
        
        Anew(1,j)=Anew(2,j);      %left boundary
        Anew(n+1,j)=Anew(n,j);    %right boundary
        
        Rnew(1,j)= Rnew(2,j);      %left boundary
        Rnew(n+1,j)=Rnew(n,j);    %right boundary
    end
    for i=2:n
        Snew(i,1)=Snew(i,2);      %bottom boundary
        Snew(i,n+1)= Snew(i,n);    %top boundary 
        
        Anew(i,1)=Anew(i,2);      %bottom boundary
        Anew(i,n+1)=Anew(i,n);    %top boundary 
        
        Rnew(i,1)= Rnew(i,2);      %bottom boundary
        Rnew(i,n+1)=Rnew(i,n);    %top boundary 
    end
    A=[ao,ac,al];
    Z=max(A);            %Limit of z-axis in the plots%
 
    
    
    %f = figure('pos',[0 200 1920 800]);    %For creating the animation 1/3%
    Col=3;
    subplot(1,3,1)   %out of 3 plots this plot will be on left%
    mesh(x,y,Snew');  %Produces a mesh grid with Snew as function of x and y%
    
    title(['Susceptable:   ','time: ',num2str(m*k),' Months ( ',num2str(round(trapz(a:hh:b,trapz(a:hh:b,Snew'))*(10^2)*1000)),' people)'])
    shading interp;
    zlim([0 Z])     %Range of the z-axis is [0,Z]%
    caxis([0 Col])
    xlabel('Dist,x (10km)');
    ylabel('Dist,y (10km)');
    zlabel('Density, People in thousands per km^2');
    drawnow
    
    subplot(1,3,2)     %out of 3 plots this plot will be in the middle%
    mesh(x,y,Anew');    %Produces a mesh grid with Anew as function of x and y%
    title(['Addicted:   ','time: ',num2str(m*k),' Months ( ',num2str(round(trapz(a:hh:b,trapz(a:hh:b,Anew'))*(10^2)*1000)),' people)'])
    shading interp;
    caxis([0 Col])
    zlim([0 Z/1000])
    xlabel('Dist,x (10km)');
    ylabel('Dist,y (10km)');
    zlabel('Density, People in thousands per km^2');
    drawnow
   
    subplot(1,3,3)     %out of 3 plots this plot will be on the right%
    mesh(x,y,Rnew');   %Produces a mesh grid with Rnew as function of x and y%
    colorbar;
    title(['Recovered:   ','time: ',num2str(m*k),' months ( ',num2str(round(trapz(a:hh:b,trapz(a:hh:b,Rnew'))*(10^2)*1000)),' people)'])
    shading interp;
    caxis([0 Col])
    zlim([0 Z/1000])
    xlabel('Dist,x (10km)');
    ylabel('Dist,y (10km)');
    zlabel('Density, People in thousands per km^2');
    drawnow
  
    %F(m) = getframe(f);     %For creating the animation 2/3%
    
    pause=1;
    S=Snew;
    A=Anew;
    R=Rnew;
    
end      %Includes all within "m" timesteps% 

%movie2avi(F,'test.avi')
  %For creating the animation 3/3%
