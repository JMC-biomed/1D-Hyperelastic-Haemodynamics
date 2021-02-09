function [ Beta,A0,npnode,nqnode,nelem,L ] = mesh( Vessel,Vnelem,Beta_start,Beta_end,A0_start,A0_end,Length )
% MESH function creates and sets up the computational mesh
nelem=sum(Vnelem);
npnode=nelem+1; % number of pressure nodes
nqnode=sum(Vnelem+1); % number of flow nodes
A0=[A0_start,A0_end]; % Area at proximal and distal end
Beta=[Beta_start,Beta_end]; % Beta at proximal and distal end
for i=1:numel(Vessel)
    L(i)=Length(i)./Vnelem(i); % Length of an element in each vessel (uniform)
end
end

