% Turn off fancy display in Reduce
off fancy; 
% a, b, c angles of ABC, a1, b1, c1 angles of PQR
c:=pi-a-b; a1:=pi; b1:=0; c1:=pi-a1-b1;  
% d1, d2, d3, d4, d5, d6 are the delta_i angles
d1:=d; d2:=pi-b-d; d3:=b-a1+d; d4:=a+a1-d; d5:=-a+c1+d; d6:=pi-c1-d;
% Now we compute H_2(delta)
h2:=e*(sin(d1)+sin(d2))/sin(b)+(1-e)*(sin(d3)+sin(d4))/sin(c)
    	+1*(sin(d5)+sin(d6))/sin(a); 
% Simplify H_2 and get its numerator and denominator
sh2:=trigsimp(h2); nh2:=num(sh2); dh2:=den(sh2); 
% Get numerators of A_2 and B_2 in H_2
ha2:=lcof(nh2,cos(d)); hb2:=lcof(nh2, sin(d)); 
% Get A_2 and B_2
aa2:=ha2/dh2; bb2:=hb2/dh2;
% Get tau_2 squared (represented with tt2 here)
tt2:=trigsimp((sin(a)+sin(b)+sin(c))^2/(aa2^2+bb2^2)); 
% Equation of the boundary of S_0 (renaming p^2,q^2,r^2 to p1,q1,r1)
s0:=sin(c)^2*sin(a)^2*q1^2 
    + sin(a)^2*sin(b)^2*r1^2 
    + sin(b)^2*sin(c)^2*p1^2 
    - 2*sin(a)*sin(c)*sin(b)^2*cos(b)*p1*r1 
    - 2*sin(b)*sin(c)*sin(a)^2*cos(a)*q1*r1 
    - 2*sin(a)*sin(b)*sin(c)^2*cos(c)*p1*q1
    - 2*sin(c)^3*sin(a)^3*sin(b)^2*cos(b)*q1 
    - 2*sin(a)^3*sin(b)^3*sin(c)^2*cos(c)*r1 
    - 2*sin(b)^3*sin(c)^3*sin(a)^2*cos(a)*p1 
    + sin(a)^4*sin(b)^4*sin(c)^4; %  
% We check that tau_2(1,e,1-e) is in the boundary of S_0
ss2:=sub({p1=tt2, q1=e^2*tt2, r1=(1-e)^2*tt2},s0);
trigsimp(ss2);
;end;