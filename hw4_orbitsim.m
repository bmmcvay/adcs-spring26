%% ENVIRONMENTAL TORQUES
clear; close all;
a=417+6378.137; %km
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
J=diag(sortedI); %kgm^2

[r0,v0]=orb2rv(a,e,deg2rad(i),deg2rad(o),deg2rad(w),deg2rad(v));
tf=50000; %final time (s)
n=500000; %numsteps
t=linspace(0,tf,n);

%ICs
omega0=0.1*randn(3,1); %(rad/s)
rho0=[0;0;0];
q0=[1;0;0;0];
x0=[r0;v0;q0;omega0;rho0];
u=[0;0;0];

options=odeset('RelTol',1e-3,'AbsTol',1e-2);
[tmat,xmat]=ode45(@(t,x) orbitaldynamics(t,x,J,u),t,x0,options);

%orbit plot
figure
xlim([-7000,7000]);ylim([-7000,7000]);zlim([-7000,7000]);
hold on
plot3(xmat(:,1),xmat(:,2),xmat(:,3));
plot3(0,0,0,'o','Color','b','MarkerSize',50,'MarkerFaceColor','#D9FFFF');
title("Orbit around Earth (size of earth not to scale)");xlabel('x');ylabel('y');zlabel('z');
%quaternion plot
figure;hold off;
plot(tmat,xmat(:,7:10));legend('q1','q2','q3','q4');xlabel("Timestep");ylabel("q value");title("Attitude quaternion over time");subtitle(['Time = ', num2str(tf), ' s']);
angle=atan2d(norm(xmat(12:13,1)),xmat(11,1));
