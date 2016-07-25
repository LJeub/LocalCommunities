function C=div_0(A,B)
% DIV_0 pointwise division such that 0/0=0

C=A;
ind=find(A);
C(ind)=A(ind)./B(ind);

end