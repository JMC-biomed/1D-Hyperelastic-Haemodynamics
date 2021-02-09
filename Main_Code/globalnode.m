function [ Avec,Pvec,inode,Qvec,upw,dnw ]=globalnode(Vessel,Vnodes,gelem,gnode,nelem,npnode,nqnode,d1,d2,d3,valve);
%GLOBALNODE finds vectors Pvec, Qvec and Avec, which encode positions and 
%connectivities of each node (for pressure, flow, and area respectively)
% Also findd the node upwind and downwind of the current node (at the ends
% of the domain a ghost node is created and assumed to have the same value
% as the end node

Pvec=zeros(nqnode,1);
upw=zeros(nqnode,1);
Avec=[1:nqnode]';
Qvec=[1:nelem*2]';

%% initialise vector for first vessel
nodes=[gnode(1,1):gnode(1,2)]';
nel=[gelem(1,1):gelem(1,2)]';
d=[d1,d2,d3];
for V=1:numel(Vessel)
    Parent=find(d1==V|d2==V|d3==V);
    s=(d(Parent,:)>0);
    if isempty(Parent)==1
        Pvec(gnode(V,1):gnode(V,2))=nodes;
        upw(gnode(V,1):gnode(V,2))=[1;nel];
        dnw(gnode(V,1):gnode(V,2))=[nel;nel(end)];
    end
end

%% Finds first node of each element
inode=zeros(nelem,1);
p=1;
for o=1:nelem
    inode(o)=p;
    p=p+1;
    if isempty(gelem(gelem(:,2)==o))==0
        p=p+1;
    end
end

end