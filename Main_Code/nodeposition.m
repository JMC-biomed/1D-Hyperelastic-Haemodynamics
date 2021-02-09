function [ gnode,gelem ] = nodeposition( Vnelem,Vessel )
%NODEPOSITION finds the first and last, node and element, in each vessel
gnode=zeros(numel(Vnelem),2);k=0;j=0;
Vnode=Vnelem+1; vsum=sum(Vnode);
for elem=1:numel(Vnelem);
    gnode(elem,1)=k+1;
    gnode(elem,2)=k+Vnode(elem);
    gelem(elem,1)=j+1;
    gelem(elem,2)=j+Vnelem(elem);
    k=k+Vnode(elem);
    j=j+Vnelem(elem);
end
end

