%% FUNCTION JAIL
function [tauB uwheel controlpath]=controller(t,x,J) %using zac's paper
    %u is magnetic dipole moment
    q=x(1:4); %B to N
    omega=x(5:7); %B (rad/s)
    rho=x(8:10); %B (rad/s)
    r=x(11:13); %N (m)
    v=x(14:16); %N (m/s)
    range=norm(r); %(m)

    %mag field
    BvecN=50e-6*[sin(t/1000);sin(t/1000+2);sin(t/1000+6)]; %teslas @@use nigrf
    BvecB=Q(q)'*BvecN; %convert to body frame
    %Bdot=cross(omega,BvecB); %omega from gyro (Bcross)

    %control
    %omega=gyro+bias;
    k=0.0001; %wheel gain
    s=[-0.3497;-0.9330;-0.0852]; %unit sun vector in inertial frame (this is currently random) @@plot and make sure this is reasonable
    sB=Q(q)'*s; %sun vector in body frame
    hmax=0.015; %max angular momentum (kgm/s)
    hd=hmax*sB/norm(sB); %desired momentum (kgm/s)
    h=J*omega+rho; %angular momentum (kgm/s)
    umax=0.2; %max magnetic dipole (A/m^2) @@(my documentation only lists the average, what's a standard value for the max?)
    tolss=0.25; %error tolerance for spin stabilization
    tolp=0.15; %error tolerance for pointing
    uwheel=[0;0;0];
    [Jval ind]=max(max(J));
    Jmax=[0;0;1];
    if norm(Jmax-h/norm(hd))>tolss %@@why norm(hd) here and not norm(h)?
        % u=hat(BvecB)*(Jmax-h/norm(hd));
        u=hat(BvecB)*(Jmax*norm(hd)-h);
        u=umax*(u/norm(u));
        controlpath=1;
    elseif norm(sB-h/norm(hd))>tolp
        % u=hat(BvecB)*(sB-h/norm(h));
        u=hat(BvecB)*(sB*norm(hd)-h);
        u=umax*(u/norm(u));
        controlpath=2;
    else
        u=[0;0;0];
        controlpath=3;
        uwheel=-k*[0;0;omega(3)]; %proportional control to get z-axis spin to 0
    end
    tauB=cross(u,BvecB); %@@is this necessary?
end

function xdot=orbitaldynamicsBB(x,J,tau_B,uwheel)
    q=x(1:4); 
    omega=x(5:7); %(rad/s)
    rho=x(8:10);
    r=x(11:13); %(m)
    v=x(14:16); %(m/s)
    mu=3.986*10^14; %(m^3/s^2)
    range=norm(r); %(m)
    q=q/norm(q);

    qdot=0.5*G(q)*omega;
    omegadot=-J\(hat(omega)*(J*omega+rho)-tau_B-uwheel);
    rhodot=-uwheel;
    rdot=v;
    vdot=-mu*r/range^3;
    xdot=[qdot;omegadot;rhodot;rdot;vdot];
end

function [xnew,tauB,controlpath]=rk4(x,h,J,t)
    [tauB,uwheel,controlpath]=controller(t,x,J);
    h1=orbitaldynamicsBB(x,J,tauB,uwheel);
    h2=orbitaldynamicsBB(x+h/2*h1,J,tauB,uwheel);
    h3=orbitaldynamicsBB(x+h/2*h2,J,tauB,uwheel);
    h4=orbitaldynamicsBB(x+h*h3,J,tauB,uwheel);
    xnew=x+h/6*(h1+2*h2+2*h3+h4);
    xnew(1:4)=xnew(1:4)/norm(xnew(1:4));
end
%% BDOT CONTROLLER
clear;
a=417+6378.137; %(km)
e=0.002243;
i=51.6448; %(deg)
o=36.2327; %(deg)
w=27.8710; %(deg)
v=332.2560; %(deg)
m=5; %(kg)
Lx=0.3; %(m)
Ly=0.1001; %(m)
Lz=0.1; %(m)
Ix=m/12*(Ly^2+Lz^2); %(kg/m^2)
Iy=m/12*(Lx^2+Lz^2); %(kg/m^2)
Iz=m/12*(Lx^2+Ly^2); %(kg/m^2)
sortedI=sort([Ix Iy Iz]);
J=diag(sortedI);
[r0,v0]=orb2rv(a,e,deg2rad(i),deg2rad(o),deg2rad(w),deg2rad(v)); %(km, km/s)
r0=r0*1000; %convert to (m)
v0=v0*1000; %convert to (m)

%model parameters
h=0.1; %10 Hz is 0.1
n=550000;
tf=n*h; %one period is like 5500 s
tspan=linspace(0,tf,n);

%initial conditions
q0=[1;0;0;0];
omega0=0.1*randn(3,1);
rho0=zeros(3,1); %wheels - none present rn
x0=[q0;omega0;rho0;r0;v0];

%loop
xmat=zeros(16,n);
xmat(:,1)=x0;
tauBmat=zeros(3,n);
controlpathvec=zeros(1,n);
for k=1:n-1
    [xmat(:,k+1) tauBmat(:,k+1) controlpathvec(:,k+1)]=rk4(xmat(:,k),h,J,tspan(k));
end

%% PLOTTING
%quaternion
figure
subplot(2,2,1)
plot(xmat(1,:));
title('q_0');
xlabel('Timestep');ylabel('q')
subplot(2,2,2)
plot(xmat(2,:));
title('q_1');
xlabel('Timestep');ylabel('q')
subplot(2,2,3)
plot(xmat(3,:));
title('q_2');
xlabel('Timestep');ylabel('q')
subplot(2,2,4)
plot(xmat(4,:));
title('q_3');
xlabel('Timestep');ylabel('q')

%omega
figure
subplot(2,2,1)
plot(xmat(5,:));
title("\omega_x");
xlabel('Timestep');ylabel('\omega (rad/s)')
subplot(2,2,2)
plot(xmat(6,:));
title("\omega_y");
xlabel('Timestep');ylabel('\omega (rad/s)')
subplot(2,2,3)
plot(xmat(7,:));
title("\omega_z");
xlabel('Timestep');ylabel('\omega (rad/s)')

%angular momentum (body frame)
figure
subplot(2,2,1)
plot(J(1,1)*xmat(5,:));
title("h_xB");
xlabel('Timestep');ylabel('h (kgm/s)')
subplot(2,2,2)
plot(J(2,2)*xmat(6,:));
title("h_yB");
xlabel('Timestep');ylabel('h (kgm/s)')
subplot(2,2,3)
plot(J(3,3)*xmat(7,:));
title("h_zB");
xlabel('Timestep');ylabel('h (kgm/s)')

%angular momentum (inertial frame)
hNmat=zeros(3,n);
for i=1:n
    hNmat(:,i)=Q(xmat(1:4,i))*(J*xmat(5:7,i));
end
figure
subplot(2,2,1)
plot(hNmat(1,:)');
yline(-0.0052);
title("h_xN");
xlabel('Timestep');ylabel('h (kgm/s)')
subplot(2,2,2)
plot(hNmat(2,:)');
yline(-0.0140);
title("h_yN");
xlabel('Timestep');ylabel('h (kgm/s)')
subplot(2,2,3)
plot(hNmat(3,:)');
yline(-0.0013);
title("h_zN");
xlabel('Timestep');ylabel('h (kgm/s)')

%magnetorquer torque plot
figure
plot(tauBmat');
title("Torque rod activity over time")
xlabel('Timestep');ylabel('\tau (Nm)');
legend('x','y','z');

%control path plot
figure
plot(controlpathvec');
title("Controller path")
xlabel("Timestep")
ylabel("Path number")

%rho plot
figure
plot(xmat(10,:)');
title("\rho_z");
xlabel('Timestep');ylabel('h (kgm/s)')


%orbit plot
figure
xlim([-7000,7000]);ylim([-7000,7000]);zlim([-7000,7000]);
hold on
plot3(xmat(11,:)/1000,xmat(12,:)/1000,xmat(13,:)/1000);
plot3(0,0,0,'o','Color','b','MarkerSize',50,'MarkerFaceColor','#D9FFFF');
title("Orbit around Earth (size of earth not to scale)");xlabel('x (km)');ylabel('y (km)');zlabel('z (km)');