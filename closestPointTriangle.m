function [clp,u,v,w] = closestPointTriangle(p, a, b, c, u, v, w)

ab = b-a;
ac = c-a;
ap = p-a;
d1 = dot(ab,ap);
d2 = dot(ac,ap);
if (d1 <= 0.0) && (d2 <= 0.0)
    u = 1.0;
    v = 0.0;
    w = 0.0;
    clp = a;
end;

bp = p-b;
d3 = dot(ab,bp);
d4 = dot(ac,bp);
if (d3 >= 0.0) && (d4 <= d3) 
    u = 0.0;
    v = 1.0;
    w = 0.0;
    clp = b;
end;
   
vc = d1*d4 - d3*d2;
if (vc <= 0.0) && (d1 >= 0.0) && (d3 <= 0.0) 
    v = d1 / (d1 - d3);
    u = 1.0 - v;
    w = 0.0;
    clp = a + ab * v;
end;

cp = p - c;
d5 = dot(ab,cp);
d6 = dot(ac,cp);
if (d6 >= 0.0) && (d5 <= d6)
    u = 0.0;
    v = 0.0;
    w = 1.0;
    clp = c;
end;

vb = d5*d2 - d1*d6;
if (vb <= 0.0) && (d2 >= 0.0) && (d6 <= 0.0) 
    w = d2 / (d2 - d6);
    u = 1.0 - w;
    v = 0.0;	
    clp = a + ac * w;
end;

va = d3*d6 - d5*d4;
if (va <= 0.0) && ((d4 - d3) >= 0.0) && ((d5 - d6) >= 0.0) 
    w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    u = 0.0;
    v = 1.0 - w;
    clp = b + (c - b) * w;
end;

denom = 1.0 / (va + vb + vc);
v = vb * denom;
w = vc * denom;
u = 1.0 - v - w;
clp = a + ab * v + ac * w;

