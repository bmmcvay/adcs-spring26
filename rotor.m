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

%% 1.1 safe mode
[V,D]=eig(J);
d=randn(3,1);
vv=randn(3,1);
Dp=D*(eye(3)+diag(hat(d)));
Vp=V*exp(hat(vv));
Jp=Vp*Dp*Vp';
% 1.2 angular velocity
omega_sp=[0;0;10];
% 1.3 rotor momentum
omega_s=norm(omega_sp);
J_s=(omega_sp/omega_s)'*Jp*(omega_sp/omega_s);
rho_s=(1.2*max(max(Jp))-J_s)*omega_s;
omegahat=hat(omega_sp);
rho0=[omega_sp';omegahat]\[rho_s*omega_s;-omegahat*J*omega_sp];
%% 1.4 rotor implementation
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
function mat=Q(q,H)
    mat=H'*L(q)*R(q)'*H;
end
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
n=6000;
tf=n*h;

q0=[1;0;0;0];

xvec=zeros(10,n);
hvec=zeros(3,n);
tao=zeros(3,n);
xvec(:,1)=[q0;omega_sp;rho0];
figure
for i=2:n
    xvec(:,i)=rk42(xvec(1:4,i-1),xvec(5:7,i-1),xvec(8:10,i-1),h,Jp);
    % poseplot(quaternion(xvec(1:4,i)')); %comment this out if you want the omega plot
    % drawnow %this too
end
for i=1:n
    hvec(:,i)=Q(xvec(1:4,i),H)*(J*xvec(5:7,i));
    omegadot=-Jp\(hat(xvec(5:7,i))*J*xvec(5:7,i));
    tao(:,i)=Jp*omegadot+cross(xvec(5:7,i),J*xvec(5:7,i));
end
plot(xvec(5:7,:)');title('Angular velocity');xlabel('Timestep');ylabel('Omega');legend('x','y','z');