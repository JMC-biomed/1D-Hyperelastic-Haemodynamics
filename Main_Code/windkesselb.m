function [ fe ] = windkessel( Ca,P1,P2,ilb,ilc,ild,ile  )
%WINDKESSEL finds K of windkessel models
global k

fe=zeros(1,2,size(P1,3)/3);

fe(1,1,ile)=Ca(1,1,ile).*(-k(2)*P1(1,1,ilb)-k(3)*P2(1,1,ilb))-Ca(1,1,ile).*(-k(2)*P1(1,2,ilb)-k(3)*P2(1,2,ilb));
fe(1,2,ile)=-Ca(1,1,ile).*(-k(2)*P1(1,1,ilb)-k(3)*P2(1,1,ilb))+Ca(1,1,ile).*(-k(2)*P1(1,2,ilb)-k(3)*P2(1,2,ilb));

end

