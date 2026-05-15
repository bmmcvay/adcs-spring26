function tau=Bdotcontroller(t,x,wgs84) %can also have mc as output
    q=x(1:4); %B to N
    omega=x(5:7); %B (rad/s)
    rho=x(8:10); %B
    r=x(11:13); %N (m)
    v=x(14:16); %N (m/s)
    range=norm(r); %(m)

    %mag field
    % jd=juliandate([2026,5,10,0,0,t]);
    % [gmst,gast]=siderealTime(jd); %get GST
    % Q_ECI2ECEF=[cos(gast) sin(gast) 0;-sin(gast) cos(gast) 0;0 0 1]; %converts from ECI to ECEF
    % R_ECEF=Q_ECI2ECEF*r; %convert range vector from ECI to ECEF
    % [lat,lon,h]=ecef2geodetic(wgs84,R_ECEF(1),R_ECEF(2),R_ECEF(3));
    % [Bvec,H,D,I,F]=igrfmagm(h,lat,lon,decyear(2026,5,10,0,0,t)); %Bvec is in nanoteslas in north-east-down frame, Bdot is in nanoteslas/year
    % [Bvecx Bvecy Bvecz]=ned2ecefv(Bvec(1),Bvec(2),Bvec(3),lat,lon); %convert Bvec from NED to ECEF
    % Bvec=Q_ECI2ECEF'*[Bvecx;Bvecy;Bvecz]; %...and from ECEF to ECI
    % BvecB=Q(q)'*Bvec*10e-9; %...and from ECI to body. also converting nanoteslas to teslas
    BvecN=50e-6*[sin(t/1000);sin(t/1000+5);sin(t/1000+17)]; 
    BvecB=Q(q)'*BvecN;
    Bdot=cross(omega,BvecB); %omega from gyro (Bcross)

    k=1000000; %gain
    mcmax=0.02; %Amp/m^2
    mc=-k*Bdot; %magnetic moment control law
    % mc=max(min(mc,mcmax),-mcmax); %saturate torquer dipole
    tau=cross(mc,BvecB);
end