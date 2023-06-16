% Turn off fancy display in Reduce to avoid display errors
off fancy; 
% a, b, c angles of ABC; a1, b1, c1 angles of PQR
c:=pi-a-b; c1:=pi-a1-b1; 
% We develop computations related to (ii) in proof of Theorem 5
% d, d1, d2, d3, d4, d5, d6 are the delta_i angles
d1:=d; d2:=pi-b-d; d3:=b-a1+d; d4:=a+a1-d; d5:=-a+c1+d; d6:=pi-c1-d;
% Introduce H_0(delta)
h0:=sin(b1)*(sin(d1)+sin(d2))/sin(b)+sin(c1)*(sin(d3)
    +sin(d4))/sin(c)+sin(a1)*(sin(d5)+sin(d6))/sin(a); 
% Simplify H_0 and get its numerator and denominator	
sh0:=trigsimp(h0); nh0:=num(sh0); dh0:=den(sh0); 
% Get numerators of A_0 and B_0 in H_0
ha0:=lcof(nh0,cos(d)); hb0:=lcof(nh0, sin(d));
% Get A_0 and B_0, represented by aa0, bb0: 
aa0:=ha0/dh0; bb0:=hb0/dh0; 
% Show equality (6.13), that is, 
% A_0=(sin(c1)*sin(a+a1)*(sin(a)+sin(b)+sin(c)))/(sin(a)*sin(c))
trigsimp(aa0-(sin(c1)*sin(a+a1)*(sin(a)+sin(b)+sin(c)))
/(sin(a)*sin(c)));
% Get tau_0 squared, which we represent here by tt0
tt0:=trigsimp((sin(a)+sin(b)+sin(c))^2/(aa0^2+bb0^2));
% Equation of the boundary of S_0 with p1=p^2, q1=q^2 and
% r1=r^2
s0:=sin(c)^2*sin(a)^2*q1^2 
    + sin(a)^2*sin(b)^2*r1^2 
    + sin(b)^2*sin(c)^2*p1^2 
    - 2*sin(a)*sin(c)*sin(b)^2*cos(b)*p1*r1 
    - 2*sin(b)*sin(c)*sin(a)^2*cos(a)*q1*r1
    - 2*sin(a)*sin(b)*sin(c)^2*cos(c)*p1*q1 
    - 2*sin(c)^3*sin(a)^3*sin(b)^2*cos(b)*q1 
    - 2*sin(a)^3*sin(b)^3*sin(c)^2*cos(c)*r1 
    - 2*sin(b)^3*sin(c)^3*sin(a)^2*cos(a)*p1 
    + sin(a)^4*sin(b)^4*sin(c)^4; 
% tau_0(sin(a1), sin(b1), sin(c1)) is in the boundary of S_0:
ss0:=sub({p1=tt0*sin(a1)^2, q1=tt0*sin(b1)^2, r1=tt0*sin(c1)^2},s0);
trigsimp(ss0);
% We repeat similar computations for (iii) in the proof of Theorem 5
% Here a1, b1, c1 are changed to -a1, -b1, -c1
d1:=d; d2:=pi-b-d; d3:=b+a1+d; d4:=a-a1-d; d5:=-a-c1+d; d6:=pi+c1-d;
h1:=sin(-b1)*(sin(d1)+sin(d2))/sin(b)+sin(-c1)*(sin(d3)
    +sin(d4))/sin(c)+sin(-a1)*(sin(d5)+sin(d6))/sin(a);
% Now everything runs as in (ii):
sh1:=trigsimp(h1); nh1:=num(sh1); dh1:=den(sh1);
ha1:=lcof(nh1,cos(d)); hb1:=lcof(nh1, sin(d));
aa1:=ha1/dh1; bb1:=hb1/dh1;
tt1:=trigsimp((sin(a)+sin(b)+sin(c))^2/(aa1^2+bb1^2));
ss1:=sub({p1=tt1*sin(a1)^2, q1=tt1*sin(b1)^2, r1=tt1*sin(c1)^2},s0);
trigsimp(ss1);
% Here we check the value of discriminant D 
% and equalities T_0=tau_0^2, T_1=tau_1^2
% Here tt represents t^2 in the text
st0:=sub({p1=tt*sin(a1)^2, q1=tt*sin(b1)^2, r1=tt*sin(c1)^2},s0); 
% We get coefficients F_2, F_1 and F_0
f2:=lcof(st0,tt); f1:=lcof(st0-f2*tt^2,tt); f0:=st0-f2*tt^2-f1*tt;
% We compute discriminant D (disc) and check its simplified value 
disc:=trigsimp(f1^2-4*f2*f0); 
trigsimp(disc-16*sin(a)^6*sin(b)^6*sin(c)^6
*sin(a1)^2*sin(b1)^2*sin(c1)^2);
% Set sqdisc to square root of discriminant D
sqdisc:=4*sin(a)^3*sin(b)^3*sin(c)^3*sin(a1)*sin(b1)*sin(c1);
% Declare roots T_0 and T_1 of the quadratic equation
root0:=(-f1-sqdisc)/(2*f2); root1:=(-f1+sqdisc)/(2*f2);
% Check that T_0=tau_0^2, T_1=tau_1^2 (here T_0=root0, T_1=root1)
trigsimp(root0-tt0); trigsimp(root1-tt1);
% Check that root0 is in fact -f0/f1 when a1=a, b1=b, c1=c:
a1:=a; b1:=b; sr:=trigsimp((2*f0)/(-f1+sqdisc)+f0/f1);
;end;