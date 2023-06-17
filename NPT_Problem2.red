rlset reals;
f0:=(z^2 + (1 - y)^2 - z*(1 - y) = p1)
    and (x^2 + (1 - z)^2 - x*(1 - z) = q1) 
    and (y^2 + (1 - x)^2 - y*(1 - x) = r1);
f1:=sub({x=x1-y1-z1+1/2, y=x1+2*z1+1/2, z=x1+y1-z1+1/2, 
    p1=(-3*p2+q2+2*r2)/6, q1=(3*p2+q2+2*r2)/6, r1=(r2-q2)/3}, f0);
f2:=ex({x1,y1,z1},f1);
f3:=rlgsn f2; f4:=rlqe f3;
f5:=sub({p2=-p^2+q^2, q2=p^2+q^2-2*r^2, r2=p^2+q^2+r^2}, f4);
f6:=f5 and p>=0 and q>=0 and r>=0;
g1:=all({p,q,r}, f6 impl p+q+r>=3/2);
% Output for Problem 1 a):
rlcad g1;
% Output for Problem 1 b):
g2:=ex({p,q,r},f6 and m=p+q+r);
g3:=rlqe g2; g4:=rltab g3;
;end;
