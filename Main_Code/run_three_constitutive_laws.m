function [ Output_Vec,CFL_law ] = run_three_constitutive_laws( Input_Vec,lambzE )
% Main code for running each case. The code uses the case and input
% parameters from Input_Vec and simulates the problem using different
% constitutive laws. The first is the law that is presented in the paper,
% the second and third is the common non-linear constitutive relation but
% using different compliance estimates

H1=Input_Vec(1);
Ri=Input_Vec(2);
omega=Input_Vec(3);
lambzres=Input_Vec(4);
c=Input_Vec(5);
k1=Input_Vec(6);
k2=Input_Vec(7);
alpha=Input_Vec(8);
beta=Input_Vec(9);

%% Constitutive Parameters (to be varied for constitutive comparison)
kappa=360./(360.-omega);
lambz=lambzE*lambzres;

totalH = H1;

flow_scale_factor=1;

n_law=3; % number of constitutive laws used

for con_law = 1:n_law % law 1 = new law, law 2 = non-linear elastic, law 3 = power law, law 4 = linear elastic
    %% Load input file
    run('CAROTID.m');
    
    %% Element size and number of elements in each vessel
    dx=1;
    ld=Length/dx;
    vel=ceil(Length/dx);
    mvel=mod(vel,2);
    vel=vel+mvel;
    Vnelem=vel;
    
    %% Temporal mesh
    dt=1e-3; % 1e-3 is recommended for optimum accuracy/speed, but higher is fine for testing
    maxt=1.1;
    t=0:dt:maxt-dt;
    k(1)= 3/(2*dt);
    k(2)= -2/dt;
    k(3)= 1/(2*dt);
    
    %% find A0 and Beta
    A0_start=pi.*r0_start.^2;
    A0_end=pi.*r0_end.^2;
    E=700e3;
    h_start=.3;
    h_end=.3;
    Beta_start=4/3*sqrt(pi)*E*h_start;
    Beta_end=4/3*sqrt(pi)*E*h_end;
    
    %% Parameter estimation for non-linear elastic, and power laws
    if con_law==2
        A0_start=Ad;
        A0_end=Ad;
        Beta_start=2*c0sq1*rho*(Ad)^(1/4);
        Beta_end=2*c0sq1*rho*(Ad)^(1/4);
    elseif con_law==3
        A0_start=Ad;
        A0_end=Ad;
        Beta_start=2*c0sq1m*rho*(Ad)^(1/4);
        Beta_end=2*c0sq1m*rho*(Ad)^(1/4);
    end
    
    %% Input file and spacial mesh generation
    [VBeta,VA0,npnode,nqnode,nelem,L]=mesh(Vessel,Vnelem,Beta_start,Beta_end,A0_start,A0_end,Length); % Preliminary 1D mesh creation
    [gnode,gelem] = nodeposition( Vnelem,Vessel ); % Find first and last node (gnode), and first and last element (gelem) of each vessel for position in global 'stiffness' matrix
    [Avec,Pvec,inode,Qvec,upw,dnw]=globalnode(Vessel,Vnelem+1,gelem,gnode,nelem,npnode,nqnode,Daughter1,Daughter2,Daughter3,Valve); % Find local element position in global 'stiffness' matrix
    [A0,Beta,L] = lineardistribution( VA0,VBeta,Vnelem,Vessel,L',gnode,nqnode,gelem ); % Linearly distribute A0, Beta, L
    [ Pwvec,Awvec,welem,wnode,bc_nodw,bc_nodw1 ] = windkesselposition( Terminal,gnode,Pvec,Avec ); % Windkessel mesh creation
    bcw=[bc_nodw,bc_nodw1]; % Boundary Conditions of Windkessel
    [ Beta,A0,L,dbdx,da0dx,inode1,nvec,nvecq,upw,dnw ] = allelem( Beta,A0,nelem,L,inode,Avec,Qvec,Pvec,upw,dnw ); % Vectorisation function
    
    %% determine size of solution vector (after taking BC's into account)
    pn=(1:Pwvec(end))';
    pn(bcw)=[];
    
    %% Position vector for each 1D node (for fast global 'stiffness' matrix constuction)
    i=1:nelem;
    I(:,1,i)=[nvec(1,1,i);nvec(1,2,:);nvec(1,1,:);nvec(1,2,:)];
    J(:,1,i)=[nvec(1,1,:);nvec(1,1,:);nvec(1,2,:);nvec(1,2,:)];
    H(:,1,i)=[nvec(1,1,i);nvec(1,2,:)];
    
    %% Windkessel Model (excluding characteristic impedance)
    ila=1:1:size(wnode,1)*3;
    ilb=1:3:size(wnode,1)*3;
    ilc=2:3:size(wnode,1)*3;
    ild=3:3:size(wnode,1)*3;
    ile=1:1:(size(wnode,1)*3/3);
    a=1;
    for j=1:size(wnode,1)
        WKvec(1,:,a)=[wnode(j,2)  wnode(j,3)];
        WKvec(1,:,a+1)=[wnode(j,2)  wnode(j,4)];
        WKvec(1,:,a+2)=[wnode(j,1)  wnode(j,2)];
        a=a+3;
    end
    
    %% Windkessel Parameters (converting to cgs units)
    R=R/100000;
    Z=Z/100000;
    C=C*100000;
    R(R<0)=[];
    Z(Z<0)=[];
    C(C<0)=[];
    
    %% Position vector for each 0D node (for fast global 'stiffness' matrix constuction)
    Iwk(:,1,ila)=[WKvec(1,1,ila);WKvec(1,2,ila);WKvec(1,1,ila);WKvec(1,2,ila)];
    Jwk(:,1,ila)=[WKvec(1,1,ila);WKvec(1,1,ila);WKvec(1,2,ila);WKvec(1,2,ila)];
    Hwk(:,1,ile)=[WKvec(1,1,ilb);WKvec(1,2,ilb)];
    
    %% Windkessel elements
    % Arterial Compliance Elements
    Kwk(1,1,ilb)= k(1)*C(ile);
    Kwk(2,1,ilb)=-k(1)*C(ile);
    Kwk(3,1,ilb)=-k(1)*C(ile);
    Kwk(4,1,ilb)= k(1)*C(ile);
    
    % Vascular Bed Resistance
    Kwk(1,1,ilc)= 1./R(ile);
    Kwk(2,1,ilc)=-1./R(ile);
    Kwk(3,1,ilc)=-1./R(ile);
    Kwk(4,1,ilc)= 1./R(ile);
    
    % Venous Compliance Elements
    Kwk(1,1,ild)= 1./Z(ile);
    Kwk(2,1,ild)=-1./Z(ile);
    Kwk(3,1,ild)=-1./Z(ile);
    Kwk(4,1,ild)= 1./Z(ile);
    C1(1,1,ile)=C(ile);
    
    %% Free up memory
    clearvars Vessel Vnelem Beta_start Beta_end A0_start A0_end npnode nqnode Valve ...
        h_start h_end dx vel BP Fraction VBeta VA0 il dbdx da0dx inode1
    
    %% Material properties
    rho = 1.06;
    mu = .04;
    
    %% Boundary and Initial Conditions
    p0=10.933e4;
    
    if con_law~=1
        p0=Pd;
    end
    
    %% Pressure or flow rate / velocity wave input using spline interpolation from digitised pressure wave
    T=maxt;
    t=t+.055;
    qbc=(6.5+3.294*sin(2*pi*t/T-0.023974)+1.9262*sin(4*pi*t/T-1.1801)-1.4219*sin(6*pi*t/T+0.92701)-0.66627*sin(8*pi*t/T-0.24118)-0.33933*sin(10*pi*t/T-0.27471)-0.37914*sin(12*pi*t/T-1.0557)+0.22396*sin(14*pi*t/T+1.22)+0.1507*sin(16*pi*t/T+1.0984)+0.18735*sin(18*pi*t/T+0.067483)+0.038625*sin(20*pi*t/T+0.22262)+0.012643*sin(22*pi*t/T-0.10093)-0.0042453*sin(24*pi*t/T-1.1044)-0.012781*sin(26*pi*t/T-1.3739)+0.014805*sin(28*pi*t/T+1.2797)+0.012249*sin(30*pi*t/T+0.80827)+0.0076502*sin(32*pi*t/T+0.40757)+0.0030692*sin(34*pi*t/T+0.195)-0.0012271*sin(36*pi*t/T-1.1371)-0.0042581*sin(38*pi*t/T-0.92102)-0.0069785*sin(40*pi*t/T-1.2364)+0.0085652*sin(42*pi*t/T+1.4539)+0.0081881*sin(44*pi*t/T+0.89599)+0.0056549*sin(46*pi*t/T+0.17623)+0.0026358*sin(48*pi*t/T-1.3003)-0.0050868*sin(50*pi*t/T-0.011056)-0.0085829*sin(52*pi*t/T-0.86463));
    qbc=qbc*flow_scale_factor;
    t=t-.055;
    
    %% Initialisation of variables
    % pressure
    P = p0*ones(Pwvec(end),1); % pressure solution at t_(n+1)
    P(bcw) = 0;    % extravascular and venous pressure (boundary condition)
    Pm1 = P;                 % pressure solution at t_n
    Pm2 = P;                % pressure solution at t_(n-1)
    p_icycle = max(P);
    p_icycle_min = min(P);
    pext=0;         % external pressure
    Pn1=P;
    
    % flow
    Q = zeros(Qvec(end),1);
    Qm1 = Q;                 % flow solution at t_n
    Qm2 = Q;                % flow solution at t_n
    qext = zeros(Pwvec(end),1);           % external (prescribed) flows
    Qn1=Q;
    
    if con_law==1
        A=A0*0+0.25;
        Alm=A(:);
        
        %% Finding initial area and compliance
        epsilon=1e-8;
        P_n=P(nvec);
        P_n=P_n(:);
        An=Alm;
        deltap = P_n - Integrate_Pressure(Alm,Ri,totalH,k1,k2,c,beta,alpha,kappa,lambz);
        deltaA=deltap;
        nus=0;
        restrict=1;
        
        %% Secant method to find Area from pressure integral
        while norm(deltaA)>1e-8
            nus=nus+1;
            p_plus=Integrate_Pressure(Alm+epsilon,Ri,totalH,k1,k2,c,beta,alpha,kappa,lambz);
            p_minus=Integrate_Pressure(Alm-epsilon,Ri,totalH,k1,k2,c,beta,alpha,kappa,lambz);
            diff_p=(p_plus-p_minus)./2/epsilon;
            Alm = An + restrict*(P_n - Integrate_Pressure(Alm,Ri,totalH,k1,k2,c,beta,alpha,kappa,lambz))./diff_p;
            deltaA=Alm-An;
            An=Alm;
            deltap = P_n - Integrate_Pressure(Alm,Ri,totalH,k1,k2,c,beta,alpha,kappa,lambz);
            if nus==100
                restrict=0.01;
                disp('check area convergence in time loop')
                %pause
            end
        end
        
        A=An(nvecq);
        Alm=A(:);
        
        %% Compliance calculation
        p_plus=Integrate_Pressure(Alm+epsilon,Ri,totalH,k1,k2,c,beta,alpha,kappa,lambz);
        p_minus=Integrate_Pressure(Alm-epsilon,Ri,totalH,k1,k2,c,beta,alpha,kappa,lambz);
        Ca1(:,1)=2*epsilon./(p_plus-p_minus);
        Ca=Ca1(nvecq);
        
        c0sq1_initial=rho*A./Ca;
        Wave_Speed=sqrt(A./rho./Ca);
        
    elseif con_law==2
        A=A0;
        Alm=A(:);
        
        % Compliance calculation
        Ca = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*Beta )./Beta./Beta;   % initial Compliance
        Ca_p = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*(Beta+epsilon) )./(Beta+epsilon)./(Beta+epsilon);
        Ca_m = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*(Beta-epsilon) )./(Beta-epsilon)./(Beta-epsilon);
        Diff_Ca=(Ca_p-Ca_m)/2/epsilon;
        if norm(Diff_Ca(:))<1e-10
            Ca_p = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*(Beta+epsilon*100) )./(Beta+epsilon*100)./(Beta+epsilon*100);
            Ca_m = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*(Beta-epsilon*100) )./(Beta-epsilon*100)./(Beta-epsilon*100);
            Diff_Ca=(Ca_p-Ca_m)/2/epsilon/100;
        end
        Beta=Beta + (Cd-Ca)./Diff_Ca;
        Ca = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*Beta )./Beta./Beta;
        delta_Ca=Cd-Ca;
        np=norm(delta_Ca(:));
        while norm(delta_Ca(:))>1e-20
            Ca_p = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*(Beta+epsilon) )./(Beta+epsilon)./(Beta+epsilon);
            Ca_m = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*(Beta-epsilon) )./(Beta-epsilon)./(Beta-epsilon);
            Diff_Ca=(Ca_p-Ca_m)/2/epsilon;
            if norm(Diff_Ca(:))<1e-10
                Ca_p = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*(Beta+epsilon*100) )./(Beta+epsilon*100)./(Beta+epsilon*100);
                Ca_m = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*(Beta-epsilon*100) )./(Beta-epsilon*100)./(Beta-epsilon*100);
                Diff_Ca=(Ca_p-Ca_m)/2/epsilon/100;
            end
            Beta=Beta + (Cd-Ca)./Diff_Ca;
            Ca = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*Beta )./Beta./Beta;
            delta_Ca=Cd-Ca;
            np=norm(delta_Ca(:));
        end
        
        c0sq2=A(1)./rho/Ca(1);
        
        Wave_Speed=sqrt(A./rho./Ca);
        %% Secant method to determine initial Beta
        
        
    elseif con_law==3
        A=A0;
        Alm=A(:);
        
        % Compliance calculation
        Ca = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*Beta )./Beta./Beta;   % initial Compliance
        Ca_p = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*(Beta+epsilon) )./(Beta+epsilon)./(Beta+epsilon);
        Ca_m = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*(Beta-epsilon) )./(Beta-epsilon)./(Beta-epsilon);
        Diff_Ca=(Ca_p-Ca_m)/2/epsilon;
        if norm(Diff_Ca(:))<1e-10
            Ca_p = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*(Beta+epsilon*100) )./(Beta+epsilon*100)./(Beta+epsilon*100);
            Ca_m = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*(Beta-epsilon*100) )./(Beta-epsilon*100)./(Beta-epsilon*100);
            Diff_Ca=(Ca_p-Ca_m)/2/epsilon/100;
        end
        Beta=Beta + (Cdm-Ca)./Diff_Ca;
        Ca = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*Beta )./Beta./Beta;
        delta_Ca=Cdm-Ca;
        np=norm(delta_Ca(:));
        while norm(delta_Ca(:))>1e-20
            Ca_p = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*(Beta+epsilon) )./(Beta+epsilon)./(Beta+epsilon);
            Ca_m = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*(Beta-epsilon) )./(Beta-epsilon)./(Beta-epsilon);
            Diff_Ca=(Ca_p-Ca_m)/2/epsilon;
            if norm(Diff_Ca(:))<1e-10
                Ca_p = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*(Beta+epsilon*100) )./(Beta+epsilon*100)./(Beta+epsilon*100);
                Ca_m = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*(Beta-epsilon*100) )./(Beta-epsilon*100)./(Beta-epsilon*100);
                Diff_Ca=(Ca_p-Ca_m)/2/epsilon/100;
            end
            Beta=Beta + (Cdm-Ca)./Diff_Ca;
            Ca = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*Beta )./Beta./Beta;
            delta_Ca=Cdm-Ca;
            np=norm(delta_Ca(:));
        end
        
        c0sq2=A(1)./rho/Ca(1);
        Wave_Speed=sqrt(A./rho./Ca);
    end
    
    %% Adapt to CFL close to 1
    CFL=max(max(Wave_Speed*dt./[L,L]));
    dt=dt;
    
    t=0:dt:maxt-dt;
    k(1)= 3/(2*dt);
    k(2)= -2/dt;
    k(3)= 1/(2*dt);
    
    t=t+.055;
    qbc=(6.5+3.294*sin(2*pi*t/T-0.023974)+1.9262*sin(4*pi*t/T-1.1801)-1.4219*sin(6*pi*t/T+0.92701)-0.66627*sin(8*pi*t/T-0.24118)-0.33933*sin(10*pi*t/T-0.27471)-0.37914*sin(12*pi*t/T-1.0557)+0.22396*sin(14*pi*t/T+1.22)+0.1507*sin(16*pi*t/T+1.0984)+0.18735*sin(18*pi*t/T+0.067483)+0.038625*sin(20*pi*t/T+0.22262)+0.012643*sin(22*pi*t/T-0.10093)-0.0042453*sin(24*pi*t/T-1.1044)-0.012781*sin(26*pi*t/T-1.3739)+0.014805*sin(28*pi*t/T+1.2797)+0.012249*sin(30*pi*t/T+0.80827)+0.0076502*sin(32*pi*t/T+0.40757)+0.0030692*sin(34*pi*t/T+0.195)-0.0012271*sin(36*pi*t/T-1.1371)-0.0042581*sin(38*pi*t/T-0.92102)-0.0069785*sin(40*pi*t/T-1.2364)+0.0085652*sin(42*pi*t/T+1.4539)+0.0081881*sin(44*pi*t/T+0.89599)+0.0056549*sin(46*pi*t/T+0.17623)+0.0026358*sin(48*pi*t/T-1.3003)-0.0050868*sin(50*pi*t/T-0.011056)-0.0085829*sin(52*pi*t/T-0.86463));
    qbc=qbc*flow_scale_factor;
    t=t-.055;
    
    Kwk(1,1,ilb)= k(1)*C(ile);
    Kwk(2,1,ilb)=-k(1)*C(ile);
    Kwk(3,1,ilb)=-k(1)*C(ile);
    Kwk(4,1,ilb)= k(1)*C(ile);
    
    %% Points which are monitored
    % For Q
    monptsq=[2*gelem(:,1)-1;(round((2*gelem(:,1)-1+2*gelem(:,2))/2));2*gelem(:,2)];  % Midpoint Q
    % For A
    monptsa=(round((gelem(:,1)+gelem(:,2))/2));   % Midpoint A
    % For P
    monptsp=[gnode(:,1);round((gnode(:,1)+gnode(:,2))/2);gnode(:,2)];   % Midpoint P
    
    %% Other monitoring points based on vessel length
    dmonitor=1;
    ld=Length/dmonitor;
    vel=ceil(Length/dmonitor);
    mvel1=mod(vel,2);
    vel=vel+mvel1;
    Vnelem=vel;
    vm=0;
    for vmonitor=1:numel(Length)
        nmonitor=Vnelem(vmonitor);
        vm=vm(end)+1:vm(end)+nmonitor+1;
        gh=(numel(gelem(vmonitor,1):gelem(vmonitor,2)))/nmonitor;
        nmonitorp=round(gnode(vmonitor,1):(numel(gelem(vmonitor,1):gelem(vmonitor,2)))/nmonitor:gnode(vmonitor,2));
        nmonitorq=round([(2*gelem(vmonitor,1)-1:(numel(2*gelem(vmonitor,1)-1:2*gelem(vmonitor,2)))/nmonitor:2*gelem(vmonitor,2)),2*gelem(vmonitor,2)]);
        MonitorQ(vm,1)=nmonitorq;
        MonitorP(vm,1)=nmonitorp;
        VesselMonitor(vm,1)=vmonitor;
        mid_p=ceil((gnode(vmonitor,1)+gnode(vmonitor,2))/2);
        mid_q=ceil((2*gelem(vmonitor,1)-1+2*gelem(vmonitor,2))/2);
    end
    ncycle=3;
    
    %% Initialise solution space
    sol_outputq = zeros(numel(monptsq),numel(t),ncycle);
    sol_outputp = zeros(numel(monptsp),numel(t),ncycle);
    sol_outputa = zeros(numel(monptsq),numel(t),ncycle);
    sol_outputCa= zeros(numel(monptsq),numel(t),ncycle);
    differenceP = zeros(1,numel(t),ncycle);
    if con_law==1;
        All_Pressure=zeros(numel(monptsq),numel(t),n_law);
        All_Flow=zeros(numel(monptsq),numel(t),n_law);
        All_Area=zeros(numel(monptsq),numel(t),n_law);
        All_Compliance=zeros(numel(monptsq),numel(t),n_law);
        All_Wave_Speed=zeros(numel(monptsq),numel(t),n_law);
    end
    
    P_f = zeros(max(max(gelem))*2,numel(t),1);
    P_b = zeros(max(max(gelem))*2,numel(t),1);
    U_f = zeros(max(max(gelem))*2,numel(t),1);
    U_b = zeros(max(max(gelem))*2,numel(t),1);
    Wi_f=zeros(max(max(gelem))*2,numel(t),1);
    Wi_b=zeros(max(max(gelem))*2,numel(t),1);
    
    sol_Q=zeros(numel(MonitorQ),numel(t));
    sol_P=zeros(numel(MonitorP),numel(t));
    
    DiffQ=1;
    DiffP=1;
    
    %% Start of main code
    eps = 1;
    eps1=1;
    eps2=1;
    eps_q=1;
    cnt = 0;
    max_cycle=20; % number of cardiac cycles
    max_q = 0;
    max_iter=1;
    maxcfl=0;
    
    %% Convergence
    tol = 1e-6;      % convergence tolerance for periodicity of cardiac cycle (percentage)
    q_icycle=0;
    
    % Loop over cardiac cycles (several cycles required to converge to periodic solution)
    while (eps > tol || eps1 > tol || eps_q > tol || eps2 > tol) & cnt<max_cycle
        
        cnt = cnt + 1;
        maxpnorm = 0;
        %tic
        
        %% loop over timesteps
        for j=1:numel(t);
            %% Inflow
            qext(1)=qbc(j);
            iter=0;
            DiffQ=1;
            DiffP=1;
            while (DiffQ > tol || DiffP > tol) & iter<max_iter
                iter=iter+1;
                %% Calculates all (q^2/a) at all nodes
                dqdx=(Qn1(nvecq).^2)./(A);
                
                %% Calculates d((q^2/a))/dx within element (upwind for 2nd node in element)
                dqdx=(dqdx(1,2,:)-dqdx(1,1,:))./L(1,1,:);
                
                %% Element Stiffness matices
                % upwind dqdx
                %[Ke,fe] = Elementstiff_Full_Newton( L, Qm1(nvecq), Qm2(nvecq), Pm1(nvec), Pm2(nvec), A, Ca, nvec, dqdx(upw),Qn1(nvecq),k,rho,mu );
                
                % central difference dqdx
                [Ke,fe] = Elementstiff_Full_Newton( L, Qm1(nvecq), Qm2(nvecq), Pm1(nvec), Pm2(nvec), A, Ca, nvec, (dqdx(upw)+dqdx(dnw))/2,Qn1(nvecq),k,rho,mu );
                
                %% Windkessel and Vessel Junctions force vectors
                [ fewk ] = windkesselb( C1,Pm1(WKvec),Pm2(WKvec),ilb,ilc,ild,ile );
                
                %% Assemble Globals Matrices and Vectors do after lagrange and windkessel
                K=sparse([I(:);Iwk(:)],[J(:);Jwk(:)],[Ke(:);Kwk(:)],max(Iwk(:)),max(Iwk(:)));
                f=sparse([H(:);Hwk(:)],1,[fe(:);fewk(:)],max(Iwk(:)),1);
                
                %% Implement bc for Windkessel
                f = f - K(:,bcw)*P(bcw);
                f(bcw) = [];
                K(bcw,:) = [];
                K(:,bcw) = [];
                
                %% Implement flow rate boundary Conditions
                f(1)=f(1)+qext(1);
                
                %% Solve the system
                Pp = K\f;
                P(pn)=Pp;
                
                %% Calculate flow rates Q
                Pnode=P(nvec(1,:,i));
                Q(nvecq)=[Ke(1,1,i).*Pnode(1,1,i)+Ke(1,2,i).*Pnode(1,2,i)-fe(1,:,i),-(Ke(2,1,i).*Pnode(1,1,i)+Ke(2,2,i).*Pnode(1,2,i)-fe(2,:,i))];
                
                if con_law==1
                    P_n=P(nvec);
                    P_n=P_n(:);
                    An=Alm;
                    deltap = P_n - Integrate_Pressure(Alm,Ri,totalH,k1,k2,c,beta,alpha,kappa,lambz);
                    deltaA=deltap;
                    nus=0;
                    restrict=1;
                    while norm(deltaA)>1e-8
                        nus=nus+1;
                        p_plus=Integrate_Pressure(Alm+epsilon,Ri,totalH,k1,k2,c,beta,alpha,kappa,lambz);
                        p_minus=Integrate_Pressure(Alm-epsilon,Ri,totalH,k1,k2,c,beta,alpha,kappa,lambz);
                        diff_p=(p_plus-p_minus)./2/epsilon;
                        Alm = An + restrict*(P_n - Integrate_Pressure(Alm,Ri,totalH,k1,k2,c,beta,alpha,kappa,lambz))./diff_p;
                        deltaA=Alm-An;
                        An=Alm;
                        deltap = P_n - Integrate_Pressure(Alm,Ri,totalH,k1,k2,c,beta,alpha,kappa,lambz);
                        if nus==100
                            restrict=0.01;
                            disp('check area convergence in time loop')
                            %pause
                        end
                    end
                    A=An(nvecq);
                    Alm=A(:);
                    
                    p_plus=Integrate_Pressure(Alm+epsilon,Ri,totalH,k1,k2,c,beta,alpha,kappa,lambz);
                    p_minus=Integrate_Pressure(Alm-epsilon,Ri,totalH,k1,k2,c,beta,alpha,kappa,lambz);
                    Ca1(:,1)=2*epsilon./(p_plus-p_minus);
                    Ca=Ca1(nvecq);
                    Calm=Ca(:);
                    
                    Wave_Speed=sqrt(A./rho./Ca);
                    
                elseif con_law==2
                    A=(A0.*(P(nvec)-p0-pext)./Beta+sqrt(A0)).^2;    % initial Area
                    Ca = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*Beta )./Beta./Beta;
                    Alm=A(:);
                    Calm=Ca(:);
                    Wave_Speed=sqrt(A./rho./Ca);
                elseif con_law==3
                    A=(A0.*(P(nvec)-p0-pext)./Beta+sqrt(A0)).^2;    % initial Area
                    Ca = 2*A0.*( A0.*(P(nvec)-pext-p0) + sqrt(A0).*Beta )./Beta./Beta;
                    Alm=A(:);
                    Calm=Ca(:);
                    Wave_Speed=sqrt(A./rho./Ca);
                end
                DiffQ=norm(Qn1-Q);
                DiffP=norm(Pn1(pn)-P(pn));
                
                Qn1=Q;
                Pn1=P;
                
            end
            
            %% Update old solutions
            Am1=A;
            Pm2 = Pm1;
            Pm1 = P;
            Qm2 = Qm1;
            Qm1 = Q;
            
            %% Monitor value of zeta (determines effective area at aortic valve)
            if mod(cnt,ncycle)==0;
                pnt=ncycle;% point to plot
                pnt1=pnt-1;
            else
                pnt=mod(cnt,ncycle);
                pnt1=pnt-1;
            end
            Wave_Speed_vec=Wave_Speed(:);
            
            %% monitor points and draw
            sol_outputq(:,j,pnt) = Q(monptsq);
            sol_outputp(:,j,pnt) = P(monptsp);
            sol_outputa(:,j,pnt) = Alm(monptsq);
            sol_outputCa(:,j,pnt) = Calm(monptsq);
            differenceP(:,j,pnt)=P(gnode(1,1))-P(gnode(1,2));
            All_Pressure(:,j,con_law)=P(monptsp)/1333.22;
            All_Flow(:,j,con_law)=Q(monptsq);
            All_Area(:,j,con_law)=Alm(monptsq);
            All_Compliance(:,j,con_law)=Calm(monptsq);
            All_Wave_Speed(:,j,con_law)=Wave_Speed_vec(monptsq);
            
            wavespeed=Q(nvecq)./A+Wave_Speed;
            CFL=wavespeed*dt./[L,L];
            maxcfl=max(maxcfl,max(max(CFL)));
            
            
            %% Monitor for tree graph
            sol_Q(:,j)=Q(MonitorQ);
            sol_P(:,j)=P(MonitorP);
            
        end
        
        if mod(cnt-1,ncycle)==0;
            pnt1=3;
        end
        %% Calculate the error compared to the previous time step for convergence criteria
        eps=abs((max(sol_outputp(2,:,pnt)-p_icycle))./max(sol_outputp(2,:,pnt)));
        eps1=abs((min(sol_outputp(2,:,pnt)-p_icycle_min))./min(sol_outputp(2,:,pnt)));
        eps2=abs(((sol_outputp(2,end,pnt)-(sol_outputp(2,end,pnt1))))./(sol_outputp(2,end,pnt)));
        p_icycle=max(sol_outputp(2,:,pnt));
        p_icycle_min=min(sol_outputp(2,:,pnt));
        eps_q=(max(sol_outputq(2,:,pnt)-q_icycle))./max(sol_outputq(2,:,pnt));
        q_icycle=max(sol_outputq(2,:,pnt));
        
    end
    
    [PWV_pos,Pos_pwv1]=min(sol_outputp(1,:,pnt));
    [PWV_pos2,Pos_pwv2]=min(sol_outputp(3,:,pnt));
    DT=t(Pos_pwv2)-t(Pos_pwv1);
    DX=sum(L);
    PWV=DX/DT;
    
    m_1=sol_outputp(1,3:numel(t),pnt)-sol_outputp(1,1:numel(t)-2,pnt);
    m_2=sol_outputp(3,3:numel(t),pnt)-sol_outputp(3,1:numel(t)-2,pnt);
    [max1,max1i]=max(m_1);
    [max2,max2i]=max(m_2);
    c1=sol_outputp(1,max1i,pnt)-max1*t(max1i)/(2*dt);
    c1a=sol_outputp(1,Pos_pwv1,pnt);
    intersecting1=(c1a - c1)/(max1/(2*dt));
    
    c2=sol_outputp(3,max2i,pnt)-max2*t(max2i)/(2*dt);
    c2a=sol_outputp(3,Pos_pwv2,pnt);
    intersecting2=(c2a - c2)/(max2/(2*dt));
    
    DT1=t(max2i)-t(max1i);
    DT2=intersecting2-intersecting1;
    DX=sum(L);
    PWVb=DX/DT2;
    
    if con_law==1;
        Pd=min(sol_outputp(1,:,pnt));
        Ad=min(sol_outputa(1,:,pnt));
        pos=find(sol_outputa(1,:,pnt)==min(sol_outputa(1,:,pnt)));
        Cd=min(sol_outputCa(1,pos,pnt));
        Cdm=mean(mean(sol_outputCa(:,:,pnt)));
        c0sq1=rho*Ad./Cd;
        c0sq1m=rho*Ad./Cdm;
    end
    
    mean_pwv(1,con_law)=PWV;
    mean_pwvb(1,con_law)=PWVb;
    mean_ws(:,con_law)=mean(All_Wave_Speed(2,:,con_law));
    max_ws(:,con_law)=max(All_Wave_Speed(2,:,con_law));
    min_ws(:,con_law)=min(All_Wave_Speed(2,:,con_law));
    mean_p(:,con_law)=mean(All_Pressure(2,:,con_law));
    max_p(:,con_law)=max(All_Pressure(2,:,con_law));
    min_p(:,con_law)=min(All_Pressure(2,:,con_law));
    mean_a(:,con_law)=mean(All_Area(2,:,con_law));
    max_a(:,con_law)=max(All_Area(2,:,con_law));
    min_a(:,con_law)=min(All_Area(2,:,con_law));
    mean_ca(:,con_law)=mean(All_Compliance(2,:,con_law));
    max_ca(:,con_law)=max(All_Compliance(2,:,con_law));
    min_ca(:,con_law)=min(All_Compliance(2,:,con_law));
    delta_p(:,con_law)=mean(differenceP(:,:,pnt));
    CFL_law(con_law,1)=maxcfl;
    
end

Output_Vec=[mean_pwv',mean_ws',mean_p',max_p',min_p',mean_a',max_a',min_a',delta_p',mean_ca',max_ca',min_ca',mean_pwvb'];

end

