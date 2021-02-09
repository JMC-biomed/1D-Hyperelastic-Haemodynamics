function [ B,A,Le,dbdx,da0dx,inode1,Avec1,Qvec1,upwind,downwind ] = allelem( Beta,A0,nelem,L,inode,Avec,Qvec,Pvec,upw,dnw )
%ALLELEM This function vectorises array indices, this allows efficient
%vectorisation to be utilised
k=1:2;
for i=1:nelem;
    Le(1,1,i)=L(i);
    inode1(1,1,i)=inode(i);
    upwind(1,:,i)=[upw(inode(i)) upw(inode(i)+1)];
    downwind(1,:,i)=[dnw(inode(i)) dnw(inode(i)+1)];
    Qvec1(1,:,i)=k;
    k=k+2;
end
Avec1=Avec([inode1 inode1+1]);
B=Beta(Avec1);
A=A0(Avec1);
dbdx=(B(1,2,:)-B(1,1,:))./Le;
da0dx=(A(1,2,:)-A(1,1,:))./Le;

end

