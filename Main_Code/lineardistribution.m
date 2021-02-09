function [ Ao,Betao,L ] = lineardistribution( A0,Beta,Vnelem,Vessel,Length,gnode,nqnode,gelem )
% LINEARDISTRIBUTION assigns an area, beta value to all nodes. 
% If a vessel tapers toward the end, the vessel will taper linearly 
% Lengths are assigned for each element

Ao=zeros(nqnode,1);
Betao=zeros(nqnode,1);

for i=1:numel(Vessel)
    %% Determine if A0 and Beta change in vessel
    if A0(i,1)==A0(i,2)
    Ao(gnode(i,1):gnode(i,2))=A0(i,1);
    else 
    Ao(gnode(i,1):gnode(i,2))=A0(i,1):(A0(i,2)-A0(i,1))/Vnelem(i):A0(i,2);    
    end
    
    if Beta(i,1)==Beta(i,2)
    Betao(gnode(i,1):gnode(i,2))=Beta(i,1);
    else 
    Betao(gnode(i,1):gnode(i,2))=Beta(i,1):(Beta(i,2)-Beta(i,1))/Vnelem(i):Beta(i,2);    
    end
    L(gelem(i,1):gelem(i,2),1)=Length(i);
end
end

