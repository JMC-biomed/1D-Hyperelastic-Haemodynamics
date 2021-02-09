function [ Pwvec,Awvec,welem,wnode,bc,bc1 ] = windkesselposition( Terminal,gnode,Pvec,Avec )
%WINDKESSELPOSITION takes note of nodal positions for pressure and
%flow/area for windkessel elements in global domain

Tv=find(Terminal==1);
Tn=gnode(Tv,2);
welem=numel(Tn);
Awvec=zeros(6*numel(Tn),1);
Pwvec=zeros(6*numel(Tn),1);
nP=Pvec(end)+1;
nA=Avec(end)+1;
node=1:6;
nodea=nA:nA+5;

for k=1:numel(Tn)
Pwvec(node(1))=Pvec(Tn(k)); % Windkessel first node connects to the last vessel
Pwvec(node(2):node(4))=nP;
Pwvec(node(5))=nP+1;
Pwvec(node(6))=nP+2;

wnode(k,:)=unique(Pwvec(node));
bc(k)=Pwvec(node(6));
bc1(k)=Pwvec(node(5));
Awvec(node)=nodea;
node=node+6;
nodea=nodea+6;
nP=nP+3;

end
end

