function ig = integrate_pressure(A,Ri,totalH,k1,k2,c,beta,alpha,kappa,lambz)


    zeta=[0 1/4 1/2 3/4 1]; % Split into 5 equally spaced nodes in [0,1] for Simpson's Rule with 2 integrals
    
    ri = ones(1,size(zeta,2)).*sqrt(A/pi);
    R = Ri+totalH*zeta;
    r = sqrt((R.*R-Ri*Ri)/(kappa*lambz)+ri.*ri);
    lambt=kappa*r./R;
    
    %
    f0 = 1.0./lambt/lambz;
    f1 = lambt;
    f2 = lambz;
    
    %
    beta_rad = beta*pi/180.0;
    I4=f1.*f1*cos(beta_rad)*cos(beta_rad)+f2*f2*sin(beta_rad)*sin(beta_rad);
    I1=f0.*f0+f1.*f1+f2*f2;
    Q=k2*(alpha*(I4-1).*(I4-1)+(1.0-alpha)*(I1-3).*(I1-3));
    
    ela=c/4*(f1.*f1-f0.*f0);
    fib1=(f1.*f1-f0.*f0)*k1*(1-alpha).*(I1-3).*exp(Q)+k1*alpha*f1.*f1*cos(beta_rad)*cos(beta_rad).*(I4-1.0).*exp(Q);
    fib2=fib1;
    
    cl = 4.0*(ela+fib1+fib2);
    
    %
    ig = cl*totalH./(f1*f2.*r);
    ig= (ig(:,1)*1 + ig(:,2)*4 + ig(:,3)*2 +ig(:,4)*4 + +ig(:,5)*1)/12; %% Simpson's Rule for two neighbouring integrals

end