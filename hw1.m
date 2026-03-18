% odehybrid
%% parameters
clear;close
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                               HW 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4 dynamics
%EDIT THESE
[r0,v0]=orb2rv(a,e,deg2rad(i),deg2rad(o),deg2rad(w),deg2rad(v));
tf=6000; %final time
n=500; %numsteps

x0=[r0;v0];
function xdot=blah(t,x)
    R=sqrt(x(1)^2+x(2)^2+x(3)^2);
    G=6.6743e-20; %using km^3 not m^3
    M=5.9722e24;
    xdot(1:3,1)=x(4:6);
    xdot(4:6,1)=-G*M*x(1:3)/R^3;
end
t=linspace(0,tf,n);
options=odeset('RelTol',1e-6,'AbsTol',1e-4);
[tmat,xmat]=ode45(@blah,t,x0,options);

%plot
figure
xlim([-7000,7000]);ylim([-7000,7000]);zlim([-7000,7000]);
hold on
plot3(xmat(:,1),xmat(:,2),xmat(:,3));
plot3(0,0,0,'o','Color','b','MarkerSize',50,'MarkerFaceColor','#D9FFFF');
title("Orbit around Earth (size of earth not to scale)");xlabel('x');ylabel('y');zlabel('z');

%% 5.1 eulers
q=[1;0;0;0];
omega=[0;0;1];

omegadot=-J\(hat(omega)*J*omega);
tao=J*omegadot+cross(omega,J*omega);

%% 5.2 stability
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
function x=calcs(q,omega,J)
    qdot=0.5*G(q,[0 0 0;1 0 0;0 1 0;0 0 1])*omega;
    omegadot=-J\(hat(omega)*J*omega);
    x=[qdot;omegadot];
end
function xnew=rk4(q,omega,h,J)
    h1=calcs(q,omega,J);
    h2=calcs(q+h/2*h1(1:4),omega+h/2*h1(5:7),J);
    h3=calcs(q+h/2*h2(1:4),omega+h/2*h2(5:7),J);
    h4=calcs(q+h*h3(1:4),omega+h*h3(5:7),J);
    xnew=[q;omega]+h/6*(h1+2*h2+2*h3+h4);
    xnew(1:4)=xnew(1:4)/norm(xnew(1:4));
end
H=[0 0 0;1 0 0;0 1 0;0 0 1]; %@@make this a function
T=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
h=0.1;
n=3000;
tf=n*h;

q0=[1;0;0;0];
omega0=[10;0;0]+[0.5;0.5;0.5]; %plus perturbation?

xvec=zeros(7,n);
hvec=zeros(3,n);
tao=zeros(3,n);
xvec(:,1)=[q0;omega0];
figure
for i=2:n
    xvec(:,i)=rk4(xvec(1:4,i-1),xvec(5:7,i-1),h,J);
    % poseplot(quaternion(xvec(1:4,i)')); %comment this out if you want the omega plot
    % drawnow %this too
end
for i=1:n
    hvec(:,i)=Q(xvec(1:4,i))*(J*xvec(5:7,i));
    omegadot=-J\(hat(xvec(5:7,i))*J*xvec(5:7,i));
    tao(:,i)=J*omegadot+cross(xvec(5:7,i),J*xvec(5:7,i));
end
plot(xvec(5:7,:)');title('Angular velocity');xlabel('Timestep');ylabel('Omega');legend('x','y','z');

%% 5.3 momentum sphere
hsphere=zeros(3,n);
for i=1:n
    hsphere(:,i)=(J*xvec(5:7,i))/norm(J*xvec(5:7,i));
end
[x,y,z]=sphere;
figure
surf(x,y,z,'EdgeColor','none','FaceAlpha',0.5);
hold on
plot3([1;0;-1;0],[0;0;0;0],[0;1;0;-1],'o','Color','g','MarkerSize',5,'MarkerFaceColor','g');
plot3([0;0],[1;-1],[0;0],'o','Color','r','MarkerSize',5,'MarkerFaceColor','r');
plot3(hsphere(1,:)',hsphere(2,:)',hsphere(3,:)');
title('Momentum sphere');xlabel('x');ylabel('y');zlabel('z');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                               HW 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
%% 1.1 safe mode
[V,D]=eig(J);
d=randn(3,1);
vv=randn(3,1);
Dp=D*(eye(3)+diag(hat(d)));
Vp=V*expm(hat(vv));
Jp=Vp*Dp*Vp';
% 1.2 angular velocity
omega_sp=[0;0;10*(2*pi/60)]; %rpm
% 1.3 rotor momentum
omega_s=norm(omega_sp);
J_s=(omega_sp/omega_s)'*Jp*(omega_sp/omega_s);
rho_s=(1.2*max(max(Jp))-J_s)*omega_s;
omegahat=hat(omega_sp);
rho0=[omega_sp';omegahat]\[rho_s*omega_s;-omegahat*J*omega_sp];
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
figure
for i=2:n
    xvec(:,i)=rk42(xvec(1:4,i-1),xvec(5:7,i-1),xvec(8:10,i-1),h,Jp);
    poseplot(quaternion(xvec(1:4,i)')); %comment this out if you want the omega plot
    drawnow %this too
end
for i=1:n
    hvec(:,i)=Q(xvec(1:4,i))*(J*xvec(5:7,i));
    omegadot=-Jp\(hat(xvec(5:7,i))*J*xvec(5:7,i));
    tao(:,i)=Jp*omegadot+cross(xvec(5:7,i),J*xvec(5:7,i));
end
%plot(xvec(5:7,:)');title('Angular velocity');xlabel('Timestep');ylabel('Omega');legend('x','y','z');

%% 2 spacecraft dynamics
[r0,v0]=orb2rv(a,e,deg2rad(i),deg2rad(o),deg2rad(w),deg2rad(v));
tf=6000; %final time
n=500; %numsteps

q0=[1;0;0;0];
x0=[r0;v0;q0;omega_sp;rho0];
function xdot=blah2(t,x)
    J=diag([0.00834167083;0.041666666;0.041666666666]);
    [V,D]=eig(J);
    d=randn(3,1);
    vv=randn(3,1);
    Dp=D*(eye(3)+diag(hat(d)));
    Vp=V*expm(hat(vv));
    Jp=Vp*Dp*Vp';
    R=sqrt(x(1)^2+x(2)^2+x(3)^2);
    mu=398600.5; %km?@@
    xdot(1:3,1)=x(4:6); %velocity 
    xdot(4:6,1)=-mu*x(1:3)/R^3; %accleration
    xdot(7:10,1)=0.5*G(x(7:10),[0 0 0;1 0 0;0 1 0;0 0 1])*x(11:13); %qdot
    xdot(14:16,1)=[0;0;0]; %rhodot
    xdot(11:13,1)=-J\(xdot(14:16,1)+hat(x(11:13))*(J*x(11:13)+x(14:16))); %omegadot
end
t=linspace(0,tf,n);
options=odeset('RelTol',1e-6,'AbsTol',1e-4);
[tmat,xmat]=ode45(@blah2,t,x0,options);

%plot
figure
xlim([-7000,7000]);ylim([-7000,7000]);zlim([-7000,7000]);
hold on
plot3(xmat(:,1),xmat(:,2),xmat(:,3));
plot3(0,0,0,'o','Color','b','MarkerSize',50,'MarkerFaceColor','#D9FFFF');
title("Orbit around Earth (size of earth not to scale)");xlabel('x');ylabel('y');zlabel('z');
plot(tmat,xmat(:,1:4));legend('q1','q2','q3','q4');
%% 3 covariance
P=[0.00000277777778 0 0;0 0.00000277777778 0;0 0 0.000123456790123];
qtrue=expq(0.5*randn(3,1));
Qtrue=Q(qtrue);
n=3; %number of observations
rn=randn(3,n); %position in ECI
for k=1:n
    rn(:,k)=rn(:,k)/norm(rn(:,k));
    rb=expm(hat(chol(P)*randn(3,1)))*Qtrue'*rn(:,k); %position in body with noise
    rb(:,k)=rb(:,k)/norm(rb(:,k));
end
function s=residuals(rn,Q,q ,rb)
    s=rn-Q(q)*rb;
end
function mat=expq(phi)
    mat=[cos(norm(phi));phi*sinc(norm(phi)/pi)];
end
%% 4 Wahba's problem
%q-method
H=[0 0 0;1 0 0;0 1 0;0 0 1];
w=[1;0.4;0.5]; %weights - pick ur favorite numbers
for k=1:n
    D=D+w(k)*transpose(L(H*rn(:,k)))*R(H*rb(:,k));
end
[v,d]=eig(D);
[womp whoomp]=max(max(d));
q=D(:,whoomp);

%SVD
for k=1:length(rn)
    B=B+w(k)*rb(:,k)*rn(:,k)';
end
[U,S,V]=svd(B);
Qmat=V*U';