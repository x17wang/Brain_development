function [l1,l2,l3] = EV(X)

c1 = X(1,1)*X(2,2) + X(1,1)*X(3,3) + X(2,2)*X(3,3) - X(1,2)*X(1,2) - X(2,3)*X(2,3) - X(1,3)*X(1,3);
c0 = X(3,3)*X(1,2)*X(1,2) + X(1,1)*X(2,3)*X(2,3) + X(2,2)*X(1,3)*X(1,3) - X(1,1)*X(2,2)*X(3,3) - 2.0*X(1,3)*X(1,2)*X(2,3);
p = trace(X)*trace(X) - 3.0*c1;
q = trace(X)*(p - 3.0/2.0*c1) - 27.0/2.0*c0;

phi = 27.0 * (0.25*c1*c1*(p-c1) + c0*(q + 27.0/4.0*c0));
phi = 1.0/3.0 * atan2(sqrt(abs(phi)), q);
t = sqrt(abs(p))*cos(phi);
s = 1.0/sqrt(3.0)*sqrt(abs(p))*sin(phi);

l3 = 1.0/3.0*(trace(X) - t) - s;
l2 = l3 + 2.0*s;
l1 = l3 + t + s;	
	