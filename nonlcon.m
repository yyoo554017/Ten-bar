function [g,geq]= nonlcon(r)
[sigma,Q]=sol_TenBarTruss(r(1),r(2));
g(1)=-(min(sigma)+2.5*10^8);
g(2)=max(sigma)-2.5*10^8;
g(3)=sqrt(Q(3)^2+Q(4)^2)-0.02;


geq=[];