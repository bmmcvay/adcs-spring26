%% FUNCTION JAIL
function [xdot tau]=dynamics(x,J,u)
    q=x(1:4);
    omega=x(5:7);
    rho=x(8:10);
    r=x(11:13);
    v=x(14:16);
    m=5; %(kg)
    mu=3.986*10^14; %(m^3/s^2)
    range=norm(r);
    q=q/norm(q);

    %drag force
    Cd=2.2;
    density0=2.72e-12; %(kg/m^3)
    SH=582000; %scale height (m) (from SMAD)
    density=@(h) density0*exp((-h-6778137)/SH);
    normalvectors=[1 -1 0 0 0 0;0 0 1 -1 0 0;0 0 0 0 1 -1]; %vector normal to each face
    centroid=[0.15 -0.15 0 0 0 0;0 0 0.05 -0.05 0 0;0 0 0 0 0.05 -0.05]; %distances from centroid to face (m)
    area=[0.01;0.01;0.03;0.03;0.03;0.03]; %area of each face (m^2) +-x,y,z
    Fdrag=zeros(3,6);
    FdragN=zeros(3,6);
    tau_drag=zeros(3,6);
    vB=Q(q)'*v; %inertial to body (m/s)
    for i=1:6
        S=dot(vB/norm(vB),normalvectors(:,i))*area(i); %body
        if sign(S)==-1
            S=0;
        end
        Fdrag(:,i)=0.5*density(range)*vB.^2*Cd*S; %body (kgm/s^2)
        tau_drag(:,i)=cross(centroid(:,i),Fdrag(:,i)); 
        FdragN(:,i)=Q(q)*Fdrag(:,i); %rotate from body to inertial
    end

    %gravity gradient
    tau_g=3*mu/(r'*r)^(5/2)*cross(r,J*r); %kgm/s^2

    tau=[sum(tau_drag,2);tau_g];

    qdot=0.5*G(q)*omega;
    omegadot=-J\(hat(omega)*(J*omega+rho)-u-tau_g-sum(tau_drag,2));
    rhodot=u;
    rdot=v;
    vdot=-mu*r/range^3-sum(FdragN,2)/m;
    xdot=[qdot;omegadot;rhodot;rdot;vdot];
end
function [xnew tau]=rk4(x,h,J,u)
    [h1 tau]=dynamics(x,J,u);
    [h2 tt]=dynamics(x+h/2*h1,J,u);
    [h3 tt]=dynamics(x+h/2*h2,J,u);
    [h4 tt]=dynamics(x+h*h3,J,u);
    xnew=x+h/6*(h1+2*h2+2*h3+h4);
    xnew(1:4)=xnew(1:4)/norm(xnew(1:4));
end
function xp=state_prediction(x,u,h)
    q=x(1:4);
    beta=x(5:7);
    omega=u-beta;
    dq=expq(0.5*h*omega);
    xp=[L(q)*dq;beta];
end
function A=state_prediction_deriv(x,u,h)
    q=x(1:4);
    beta=x(5:7);
    omega=u-beta;
    dq=expq(0.5*h*omega);
    qn=L(q)*dq;
    A=[G(qn)'*R(dq)*G(q) -0.5*h*G(qn)'*G(q);zeros(3,3) eye(3,3)];
end
function z=innovation(x,y,m)
    q=x(1:4);
    zk=zeros(3,m);
    yk=reshape(y,4,m);
    for k=1:m
        zk(:,k)=G(q)'*yk(:,k);
    end
    z=zk(:);
end
function C=innovation_deriv(x,y,m)
    q=x(1:4);
    C=zeros(3*m,6);
    yk=reshape(y,4,m);
    for k=1:m
        C((k-1)*3+1:(k-1)*3+3,1:3)=H()'*R(yk(:,k))*T()*G(q);
    end
end
function u=steering(axes,tao)
    u=axes'*((axes*axes'+1e-3*eye(3))\tao);
end
%% SIMULATION (MEKF + PD) 
%put in filter's state into pd, (q0 and omega minus bias)
close all;clear
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
J=diag(sortedI);
mu=3.986*10^14; %m^3/s^2

%parameters
h=0.1; %10 Hz
n=55000;
tf=n*h;
tspan=linspace(0,tf,n);

%sample random initial attitude & angular velocity
q0=0.1*randn(4,1);
q0=q0/norm(q0);
omega0=0.1*randn(3,1);
rho0=zeros(3,1);
[r0,v0]=orb2rv(a,e,deg2rad(i),deg2rad(o),deg2rad(w),deg2rad(v)); %(km, km/s)
r0=r0*1000; %convert to (m)
v0=v0*1000; %convert to (m)
x0=[q0;omega0;rho0;r0;v0];
     
% Noise
m=1; %number of star trackers
W=0.001*eye(3*m,3*m); %measurement noise
V=0.000001*eye(6,6); %process noise
     
%Simulate n time steps
xtraj=zeros(16,n);
xtraj(:,1)=x0;

%RW controller stuff
axes=[1 0 0;0 1 0;0 0 1]; %normalized RW axes
qd=[1;0;0;0]; %desired attitude
omegad=[0;0;0]; %desired angular velocity

% FILTER (MEKF + PD)
% Filter Initialization
xfilt=zeros(7,n);
xfilt(1:4,1)=q0+0.01*randn(4,1); %initial "guess"
xfilt(1:4,1)=xfilt(1:4,1)/norm(xfilt(1:4,1));
ytraj=zeros(4*m,n);
gyro=zeros(3,n);
bias=zeros(3,n);
u=zeros(3,n);
bias(:,1)=0.3*randn(3,1); %initial bias
M=0.001*eye(3,3)+eye(3,3); %mounting error
M=expm(hat(0.01*randn(3,1)))*M*(expm(hat(0.01*randn(3,1))))';
taumat=zeros(6,n);

options=odeset('RelTol',1e-6,'AbsTol',1e-4);
wgs84=wgs84Ellipsoid('meter');

P=zeros(6,6,n);
P(:,:,1)=0.5*eye(6,6); %initial covariance guess
     
for k=1:(n-1)
    %Controller
    u(:,k)=steering(axes,Bdotcontroller(tspan(k),xtraj(:,k),wgs84));
    %u(:,k)=steering(axes,PDcontroller([xtraj(1:4,k);xtraj(5:7,k)],qd,omegad));

    %True state
    [xtrajk tauk]=rk4(xtraj(:,k),h,J,u(:,k));
    xtraj(:,k+1)=xtrajk;
    taumat(:,k+1)=tauk;
     
    %Sensors
    %gyro
    gyro(:,k+1)=M*xtraj(5:7,k)+bias(:,k)+sqrt(V(4:6,4:6))*randn(3,1);
    bias(:,k+1)=bias(:,k)+0.0001*randn(3,1); %random walk
    %star tracker
    bore=[0.7071 -0.7071;0.7071 0.7071;0 0]; %unit vector directions of boresight vectors
    Wab=0.00019393; %about-boresight error (rad)
    Wc=0.00002909; %cross-boresight error (rad)
    for i=1:m
        Wstar=Wab*bore(:,i)*bore(:,i)'+Wc*(eye(3,3)-bore(:,i)*bore(:,i)');
        wk=Wstar*randn(3,1); %@@ was sqrt(Wstar) but there were negative values in W
        yk(:,i)=L(xtraj(1:4,k+1))*expq(wk); %noise
    end
    ytraj(:,k+1)=yk(:);

    %Prediction
    xpred=state_prediction(xfilt(:,k),gyro(:,k+1),h); %either gyro timestep works
    A=state_prediction_deriv(xfilt(:,k),gyro(:,k+1),h);
    Ppred=A*P(:,:,k)*A'+V;

    %Innovation
    z=innovation(xpred,ytraj(:,k+1),m);
    C=innovation_deriv(xpred,ytraj(:,k+1),m);
    S=C*Ppred*C'+W;

    %Kalman Gain
    K=Ppred*C'/S;

    %Update
    dx=-K*z;
    phi=dx(1:3);
    dq=[sqrt(1-phi'*phi);phi];
    xfilt(:,k+1)=[L(xpred(1:4))*dq;xpred(5:7)+dx(4:6)];
    
    P(:,:,k+1)=(eye(6,6)-K*C)*Ppred*(eye(6,6)-K*C)'+K*W*K';
end

error=L(qd)'*xfilt(1:4,:);
rmserror=rms(error);

%% plotting - xMEKF
%quaternion
figure
subplot(2,2,1)
plot(xfilt(1,:));
hold on;
plot(xtraj(1,:));
title("q_0");legend("filter","actual");
xlabel('Timestep');ylabel('q')
hold off;
subplot(2,2,2)
plot(xfilt(2,:));
hold on;
plot(xtraj(2,:));
title("q_1");legend("filter","actual");
xlabel('Timestep');ylabel('q')
hold off;
subplot(2,2,3)
plot(xfilt(3,:));
hold on;
plot(xtraj(3,:));
title("q_2");legend("filter","actual");
xlabel('Timestep');ylabel('q')
hold off;
subplot(2,2,4)
plot(xfilt(4,:));
hold on;
plot(xtraj(4,:));
title("q_3");legend("filter","actual");
xlabel('Timestep');ylabel('q')
hold off;

%omega
figure
subplot(2,2,1)
plot(gyro(1,:)-xfilt(5,:));
hold on;
plot(xtraj(5,:));
title("\omega_x");legend("filter","actual");
xlabel('Timestep');ylabel('\omega (rad/s)')
hold off;
subplot(2,2,2)
plot(gyro(2,:)-xfilt(6,:));
hold on;
plot(xtraj(6,:));
title("\omega_y");legend("filter","actual");
xlabel('Timestep');ylabel('\omega (rad/s)')
hold off;
subplot(2,2,3)
plot(gyro(3,:)-xfilt(7,:));
hold on;
plot(xtraj(7,:));
title("\omega_z");legend("filter","actual");
xlabel('Timestep');ylabel('\omega (rad/s)')
hold off;

%torque @@make sure this isn't going over 0.004 Nm of torque 
figure
subplot(2,2,1)
plot(u(1,:));
yline(0.004);yline(-0.004);
title('u_x');
xlabel('Timestep');ylabel('Torque (Nm)')
subplot(2,2,2)
plot(u(2,:));
yline(0.004);yline(-0.004);
title('u_y');
xlabel('Timestep');ylabel('Torque (Nm)')
subplot(2,2,3)
plot(u(3,:));
yline(0.004);yline(-0.004);
title('u_z');
xlabel('Timestep');ylabel('Torque (Nm)')

%angular momentum
figure
subplot(2,2,1)
plot(xtraj(8,:));
title('\rho_x');
subplot(2,2,2)
plot(xtraj(9,:));
title('\rho_y');
subplot(2,2,3)
plot(xtraj(10,:));
title('\rho_z');

%3-parameter tracking error
figure
subplot(2,2,1)
plot(error(2,:));
title('q_1 error');
xlabel('Timestep');ylabel('Error')
subplot(2,2,2)
plot(error(3,:));
title('q_2 error');
xlabel('Timestep');ylabel('Error')
subplot(2,2,3)
plot(error(4,:));
title('q_3 error');
xlabel('Timestep');ylabel('Error')

%orbit plot
figure
xlim([-7000,7000]);ylim([-7000,7000]);zlim([-7000,7000]);
hold on
plot3(xtraj(11,:)/1000,xtraj(12,:)/1000,xtraj(13,:)/1000);
plot3(0,0,0,'o','Color','b','MarkerSize',50,'MarkerFaceColor','#D9FFFF');
title("Orbit around Earth (size of earth not to scale)");xlabel('x (km)');ylabel('y (km)');zlabel('z (km)');

%range plot
figure
ranges=sqrt(xtraj(11,:).^2+xtraj(12,:).^2+xtraj(13,:).^2)/1000;
plot(ranges);
title("Range over time")
xlabel('Timestep');ylabel('Range (km)');

%drag torque plot
figure
plot(taumat(1:3,:)');
title("Drag torque over time")
xlabel('Timestep');ylabel('\tau (Nm)');
legend('x','y','z');

%GG torque plot
figure
plot(taumat(4:6,:)');
title("Gravity gradient torque over time")
xlabel('Timestep');ylabel('\tau (Nm)');
legend('x','y','z');