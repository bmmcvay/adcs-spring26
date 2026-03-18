%% parameters
clear;close all;
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
%function jail
function mat=hat(v)
    mat=[0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0];
end
function mat=L(q)
    mat=[q(1) -q(2:4)';q(2:4) q(1)*eye(3)+hat(q(2:4))];
end
function mat=R(q)
    mat=[q(1) -q(2:4)';q(2:4) q(1)*eye(3)-hat(q(2:4))];
end
function mat=G(q,H)
    mat=L(q)*H;
end
function mat=Q(q)
    H=[0 0 0;1 0 0;0 1 0;0 0 1];
    mat=H'*L(q)*R(q)'*H;
end
%% 1.1 safe mode
[V,D]=eig(J);
d=randn(3,1);
vv=randn(3,1);
Dp=D*(eye(3)+diag(hat(d)));
Vp=V*expm(hat(vv));
Jp=Vp*Dp*Vp';
% 1.2 angular velocity
omega_sp=[0;0;10*(2*pi/60)]+[0.5;0.5;0.5]; %rad/s
% 1.3 rotor momentum
omega_s=norm(omega_sp);
J_s=(omega_sp/omega_s)'*Jp*(omega_sp/omega_s);
rho_s=(1.2*max(max(Jp))-J_s)*omega_s;
omegahat=hat(omega_sp);
rho0=[omega_sp';omegahat]\[rho_s*omega_s;-omegahat*Jp*omega_sp];
%% 1.4 rotor implementation
function x=calcs2(q,omega,rho,J)
    rhodot=[0;0;0];
    q=q/norm(q);
    qdot=0.5*G(q,[0 0 0;1 0 0;0 1 0;0 0 1])*omega;
    omegadot=-J\(rhodot+hat(omega)*(J*omega+rho));
    x=[qdot;omegadot;rhodot];
end
function xnew=rk42(q,omega,rho,h,J)
    h1=calcs2(q,omega,rho,J);
    h2=calcs2(q+h/2*h1(1:4),omega+h/2*h1(5:7),rho+h/2*h1(8:10),J);
    h3=calcs2(q+h/2*h2(1:4),omega+h/2*h2(5:7),rho+h/2*h2(8:10),J);
    h4=calcs2(q+h*h3(1:4),omega+h*h3(5:7),rho+h*h3(8:10),J);
    xnew=[q;omega;rho]+h/6*(h1+2*h2+2*h3+h4);
    xnew(1:4)=xnew(1:4)/norm(xnew(1:4));
end
H=[0 0 0;1 0 0;0 1 0;0 0 1];
T=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
h=0.1;
n=300;
tf=n*h;

q0=[1;0;0;0];

xvec=zeros(10,n);
hvec=zeros(3,n);
tao=zeros(3,n);
xvec(:,1)=[q0;omega_sp;rho0];
for i=2:n
    xvec(:,i)=rk42(xvec(1:4,i-1),xvec(5:7,i-1),xvec(8:10,i-1),h,Jp);
    % poseplot(quaternion(xvec(1:4,i)')); %comment this out if you want the omega plot
    % drawnow %this too
end
for i=1:n
    hvec(:,i)=Q(xvec(1:4,i))*(J*xvec(5:7,i));
    omegadot=-Jp\(hat(xvec(5:7,i))*J*xvec(5:7,i));
    tao(:,i)=Jp*omegadot+cross(xvec(5:7,i),J*xvec(5:7,i));
end
figure
plot(xvec(5:7,:)');title('Angular velocity');xlabel('Timestep');ylabel('Omega');legend('x','y','z');

%% 2 spacecraft dynamics
[r0,v0]=orb2rv(a,e,deg2rad(i),deg2rad(o),deg2rad(w),deg2rad(v));
tf=100; %final time
n=200; %numsteps

q0=[1;0;0;0];
x0=[r0;v0;q0;omega_sp;rho0];
function xdot=blah2(t,x,J) %create xdot from x
    R=sqrt(x(1)^2+x(2)^2+x(3)^2); %km
    mu=398600.5; %km^3/s^2
    xdot(1:3,1)=x(4:6); %velocity  km/s
    xdot(4:6,1)=-mu*x(1:3)/R^3; %accleration  km
    xdot(7:10,1)=0.5*G(x(7:10),[0 0 0;1 0 0;0 1 0;0 0 1])*x(11:13); %qdot
    xdot(14:16,1)=[0;0;0]; %rhodot
    xdot(11:13,1)=-J\(xdot(14:16,1)+hat(x(11:13))*(J*x(11:13)+x(14:16))); %omegadot
end
t=linspace(0,tf,n);
options=odeset('RelTol',1e-6,'AbsTol',1e-4);
[tmat,xmat]=ode45(@(t,x) blah2(t,x,Jp),t,x0,options);

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
%pointing plot
figure
pointingerror=zeros(n,1);
for k=1:n
    qd=quatmultiply(xmat(k,7:10),quatinv([0 1 0 0]));
    pointingerror(k)=atan2d(sqrt(qd(2)^2+qd(3)^2+qd(4)^2),qd(1));
end
plot(tmat,pointingerror);xlabel("Timestep");ylabel("Pointing error (deg)");title("Pointing error over time");subtitle(['Time = ', num2str(tf), ' s']);
%% 3 covariance
P=[0.00002909 0 0;0 0.00002909 0;0 0 0.00019393];
qtrue=expq(0.5*randn(3,1));
Qtrue=Q(qtrue);
n=3; %number of observations
rn=randn(3,n); %position in ECI
for k=1:n
    rn(:,k)=rn(:,k)/norm(rn(:,k));
    rb(:,k)=expm(hat(chol(P)*randn(3,1)))*Qtrue'*rn(:,k); %position in body with noise
    rb(:,k)=rb(:,k)/norm(rb(:,k));
end
function s=residuals(rn,Q,q,rb)
    s=rn-Q(q)*rb;
end
function mat=expq(phi)
    mat=[cos(norm(phi));phi*sinc(norm(phi)/pi)];
end
cov(rb);
%% 4 Wahba's problem
%q-method
H=[0 0 0;1 0 0;0 1 0;0 0 1];
w=[1;0.4;0.5]; %weights - pick ur favorite numbers
D=zeros(4,4);
for k=1:n
    D=D+w(k)*L(H*rn(:,k))'*R(H*rb(:,k));
end
[v,d]=eig(D);
[womp whoomp]=max(max(abs(d)));
q_daven=-v(:,whoomp)/norm(v(:,whoomp));
Q_daven=Q(q_daven);

%SVD
B=zeros(3,3);
for k=1:length(rn)
    B=B+w(k)*rb(:,k)*rn(:,k)';
end
[U,S,V]=svd(B);
Q_svd=V*U';

%comparison
percent_daven=100*abs((Qtrue-Q_daven)./Qtrue)
percent_svd=100*abs((Qtrue-Q_svd)./Qtrue)