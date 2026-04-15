%% FUNCTION JAIL
function mat=hat(v)
    mat=[0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0];
end
function mat=L(q)
    mat=[q(1) -q(2:4)';q(2:4) q(1)*eye(3)+hat(q(2:4))];
end
function mat=R(q)
    mat=[q(1) -q(2:4)';q(2:4) q(1)*eye(3)-hat(q(2:4))];
end
function mat=G(q)
    H=[0 0 0;1 0 0;0 1 0;0 0 1];
    mat=L(q)*H;
end
function mat=Q(q)
    H=[0 0 0;1 0 0;0 1 0;0 0 1];
    mat=H'*L(q)*R(q)'*H;
end
function x=dynamics(x,J)
    q=x(1:4);
    omega=x(5:7);
    q=q/norm(q);
    qdot=0.5*G(q)*omega;
    omegadot=-J\(hat(omega)*J*omega);
    x=[qdot;omegadot];
end
function xnew=rk4(x,h,J)
    h1=dynamics(x,J);
    h2=dynamics(x+h/2*h1,J);
    h3=dynamics(x+h/2*h2,J);
    h4=dynamics(x+h*h3,J);
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
    H=[0 0 0;1 0 0;0 1 0;0 0 1];
    q=x(1:4);
    C=zeros(3*m,6);
    yk=reshape(y,4,m);
    for k=1:m
        C((k-1)*3+1:(k-1)*3+3,1:3)=H'*R(yk(:,k))*T()*G(q);
    end
end
function mat=expq(phi)
    mat=[cos(norm(phi));phi*sinc(norm(phi)/pi)];
end

%% SIMULATION 
close all
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

%parameters
h=0.1;
n=600;
tf=n*h;

%sample random initial attitude & angular velocity
q0=randn(4,1);
q0=q0/norm(q0);
omega0=0.1*randn(3,1);
x0=[q0;omega0];
     
% Noise
m=2; %number of star trackers
W=0.001*eye(3*m,3*m); %measurement noise
V=0.001*eye(6,6); %process noise
     
%Simulate n time steps
xtraj=zeros(7,n);
xtraj(:,1)=x0;
for k=1:(n-1)
    xtraj(:,k+1)=rk4(xtraj(:,k),h,J);
end
     
%generate noisy gyro measurements
gyro=zeros(3,n);
bias=zeros(3,n);
bias(:,1)=0.3*randn(3,1); %initial bias
M=0.001*eye(3,3)+eye(3,3); %mounting error
M=expm(hat(0.01*randn(3,1)))*M*(expm(hat(0.01*randn(3,1))))';
for k = 1:n
    gyro(:,k)=M*xtraj(5:7,k)+bias(:,k)+sqrt(V(4:6,4:6))*randn(3,1);
    bias(:,k+1)=bias(:,k)+0.01*randn(3,1); %random walk
end
bias(:,n+1)=[]; %remove the last column because i'm good at computer science
 
%Generate noisy star-tracker measurements
ytraj=zeros(4*m,n);
for k = 1:n
    qk = xtraj(1:4,k);
    yk = zeros(4,m);
    wk = reshape(sqrt(W)*randn(3*m,1),3,m);
    for i=1:m
        yk(:,i)=L(qk)*expq(wk(:,i)); %noise
    end
    yk(3,:)=yk(3,:)+0.0001*randn(1,m); %about-boresight error
    yk(1:2,:)=yk(1:2,:)+0.00001*randn(2,m); %cross-boresight error
    ytraj(:,k)=yk(:);
end

%% FILTER
% Filter Initialization
xfilt=zeros(7,n);
xfilt(1:4,1)=q0+0.01*randn(4,1);
xfilt(1:4,1)=xfilt(1:4,1)/norm(xfilt(1:4,1));

P=zeros(6,6,n);
P(:,:,1)=0.5*eye(6,6); %initial covariance guess
     
for k=1:(n-1)
    %Prediction
    xpred=state_prediction(xfilt(:,k),gyro(:,k),h);
    A=state_prediction_deriv(xfilt(:,k),gyro(:,k),h);
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

%% VECTOR SENSORS (SVD) WITH AFFINE ERROR MODEL
%initialization
qtrue=expq(0.5*randn(3,1));
Qtrue=Q(qtrue);
n=3; %number of observations
b=0.1*randn(3,1);
w=0.1*randn(3,1); %additive noise

%vector sim
rn=randn(3,n); %position in ECI
rb=zeros(3,n);
for k=1:n
    rn(:,k)=rn(:,k)/norm(rn(:,k));
    rb(:,k)=Qtrue'*(M*rn(:,k)+b+w); %position in body with noise
    rb(:,k)=rb(:,k)/norm(rb(:,k)); %re-normalize
end

%SVD
w=[1;0.4;0.5]; %weights - pick ur favorite numbers
B=zeros(3,3);
for k=1:length(rn)
    B=B+w(k)*rb(:,k)*rn(:,k)';
end
[U,S,V]=svd(B);
Q_svd=V*U';
errorsvd=mean(mean(100*abs((Qtrue-Q_svd)./Qtrue)))

%% ERROR ANALYSIS
%error & percent error calcs
quaterror=abs(xfilt(1:4,:)-xtraj(1:4,:));
biaserror=abs(xfilt(5:7,:)-bias);
quatpercent=100*abs((xtraj(1:4,:)-xfilt(1:4,:))./xtraj(1:4,:));
biaspercent=100*abs((bias-xfilt(5:7,:))./bias);
avgquaterror=mean(mean(quatpercent))
avgbiaserror=mean(mean(biaspercent))

%matrices for finding the mean of the error across multiple trials
%run the first time only
    % matqu=[avgquaterror];
    % matbi=[avgbiaserror];
    % matsvd=[errorsvd];
%run on subsequent trials & take the mean of these once done doing trials
    matqu=[matqu avgquaterror];
    matbi=[matbi avgbiaserror];
    matsvd=[matsvd errorsvd];
    
%% PLOTTING
%quaternions
figure
subplot(2,2,1)
plot(xfilt(1,:));
hold on;
plot(xtraj(1,:));
title("q0");legend("filter","actual");
hold off;
subplot(2,2,2)
plot(xfilt(2,:));
hold on;
plot(xtraj(2,:));
title("q1");legend("filter","actual");
hold off;
subplot(2,2,3)
plot(xfilt(3,:));
hold on;
plot(xtraj(3,:));
title("q2");legend("filter","actual");
hold off;
subplot(2,2,4)
plot(xfilt(4,:));
hold on;
plot(xtraj(4,:));
title("q3");legend("filter","actual");
hold off;

%biases
figure
subplot(2,2,1)
plot(xfilt(5,:));
hold on;
plot(bias(1,:));
title("x bias");legend("filter","actual");
hold off;
subplot(2,2,2)
plot(xfilt(6,:));
hold on;
plot(bias(2,:));
title("y bias");legend("filter","actual");
hold off;
subplot(2,2,3)
plot(xfilt(7,:));
hold on;
plot(bias(3,:));
title("z bias");legend("filter","actual");
hold off;

%errors
figure
subplot(1,2,1)
plot(quaterror');
title('quaternion error');legend("q0","q1","q2","q3");
subplot(1,2,2)
plot(biaserror');
title('bias error');legend('x','y','z');

%percent errors
figure
subplot(1,2,1)
plot(quatpercent');
title('quaternion percent error');legend("q0","q1","q2","q3");
subplot(1,2,2)
plot(biaspercent');
title('bias percent error');legend('x','y','z');