function [xnew,u,tau]=rk4(x,h,J,u,wgs84,t)
    u=Bdotcontroller(t,x,wgs84);
    [h1 tau]=orbitaldynamicsB(x,J,u);
    [h2 tt]=orbitaldynamicsB(x+h/2*h1,J,u);
    [h3 tt]=orbitaldynamicsB(x+h/2*h2,J,u);
    [h4 tt]=orbitaldynamicsB(x+h*h3,J,u);
    xnew=x+h/6*(h1+2*h2+2*h3+h4);
    xnew(1:4)=xnew(1:4)/norm(xnew(1:4));
end
%%
clear; close all;
a=417000+6378137; %m
e=0.002243;
i=51.6448; %deg
o=36.2327; %deg
w=27.8710; %deg
v=332.2560; %deg
m=5; %kg
Lx=0.3; %m
Ly=0.1001; %m
Lz=0.1; %m
Ix=m/12*(Ly^2+Lz^2);
Iy=m/12*(Lx^2+Lz^2);
Iz=m/12*(Lx^2+Ly^2);
sortedI=sort([Ix Iy Iz]);
J=diag(sortedI);
[r0,v0]=orb2rv(a,e,deg2rad(i),deg2rad(o),deg2rad(w),deg2rad(v)); %(km, km/s)

tf=1000; %final time (s)
n=10000; %numsteps
h=tf/n;
t=linspace(0,tf,n);

nturns=1; %number of turns
I=10; %current (A)
A=0.1/1000; %enclosed area (km^2)  

q0=[1;0;0;0];
omega0=0.1*randn(3,1);
rho0=zeros(3,1);
u=[0;0;0]; %delete if including wheels
x0=[q0;omega0;rho0;r0;v0];

options=odeset('RelTol',1e-6,'AbsTol',1e-4);
% wgs84=wgs84Ellipsoid('meter');
wgs84=0;

xmat=zeros(16,n);
u=zeros(3,n);
tau=zeros(6,n);
xmat(:,1)=x0;
for k=1:n-1
    [xmat(:,k+1) u(:,k+1) tau(:,k+1)]=rk4(xmat(:,k),h,J,rho0,wgs84,t(k));
end
% [tmat,xmat]=ode45(@(t,x) orbitaldynamicsB(x,J,u),t,x0,options);

%% plotting 
figure
xlim([-7000,7000]);ylim([-7000,7000]);zlim([-7000,7000]);
hold on
plot3(xmat(11,:),xmat(12,:),xmat(13,:));
plot3(0,0,0,'o','Color','b','MarkerSize',50,'MarkerFaceColor','#D9FFFF');
title("Orbit around Earth (size of earth not to scale)");xlabel('x');ylabel('y');zlabel('z');
%quaternion plot
figure;hold off;
plot(xmat(1:4,:)');legend('q1','q2','q3','q4');xlabel("Timestep");ylabel("q value");title("Attitude quaternion over time");subtitle(['Time = ', num2str(tf), ' s']);
%angular momentum
% plot(t,xmat(:,11:13));
%omegaplot
figure;hold on;plot(xmat(5,:));plot(xmat(6,:));plot(xmat(7,:));legend('\omega_x','\omega_y','\omega_z');






% %controllers %@@ how do you do these? q and omega?
% function tao=hPDcontroller(x,hd,hdotd) %for wheel
%     kp=0.6;
%     kd=0.13;
%     h=x(1:4);
%     hdot=x(5:7); %@@ h part of state? where get hdot?
%     herror=h-hd; 
%     hdoterror=hdot-hdotd;
%     tao=-kp*herror-kd*hdoterror;
% end
% function tao=bPDcontroller(x,Bd,Bdotd) %for magnetorquers
%     kp=0.6;
%     kd=0.13;
%     B=x(1:4);
%     Bdot=x(5:7);
%     Berror=h-hd; 
%     Bdoterror=hdot-hdot;
%     tao=-kp*Berror-kd*Bdoterror;
% end

%@@ how to implement with rk4 since 2 controllers?
% xPD(:,1)=x0;
% for i=1:(n-1)
%     umag=steering(axes,bPDcontroller(xPD(:,i),));
%     for j=1:(n-1)
%         uwheel=steering(axes,hPDcontroller(xPD(:,i),));
%     end
%     xPD(:,i+1)=rk4(xPD(:,i),h,J,u,w,axes);
% end

%get Bcross detumbling to work and then plot by tomorrow

%gyrostat with one wheel's angular momentum, 3 magnetorquers, will have 2 control loops
% 1 for  magnetorquer to point (fixed-magnitude) ang momentum vector at a target
%inner loop adjusts 1 wheel to despin bus
%align ang momentum axis with wheel axis
% 2 pd loops