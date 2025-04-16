function dfdx = diff6(f,x)
a1 = 150/256;
a2 = -25/256;
a3 = 3/256;
b1 = 1/16;
b2 = 1/24;
c1 = 75/64;
c2 = -25/384;
c3 = 3/640;
d1 = 1.5;
d2 = -0.3;
d3 = 1/30;

nx = length(x);

for i=3:nx-3
    Q(i) = a1*(f(i+1)+f(i))+a2*(f(i+2)+f(i-1))+a3*(f(i+3)+f(i-2));
end

Q0     =b1* (35 *f(1) - 35 *f(2)   + 21 *f(3)   - 5 *f(4));
Q(1)   =b1* (5  *f(1) + 15 *f(2)   - 5  *f(3)   +    f(4));
Q(2)   =b1* (-   f(1) + 9  *f(2)   + 9  *f(3)   -    f(4));
Q(nx-2)=b1* (-   f(nx)+ 9  *f(nx-1)+ 9  *f(nx-2)-    f(nx-3));
Q(nx-1)=b1* (5  *f(nx)+ 15 *f(nx-1)- 5  *f(nx-2)+    f(nx-3));
Q(nx)  =b1* (35 *f(nx)- 35 *f(nx-1)+ 21 *f(nx-2)- 5 *f(nx-3));

for i=4:nx-3
    fx(i) = (c1*(Q(i)-Q(i-1))+c2*(Q(i+1)-Q(i-2))+c3*(Q(i+2)-Q(i-3)));
end

fx(1)   = b2* (-23  *Q0     +21*Q(1)   +3 *Q(2)   -  Q(3));
fx(2)   = b2* (-22  *Q(1)   +17*Q(2)   +9 *Q(3)   -5*Q(4)   +Q(5));
fx(3)   = b2* (      Q(1)   -27*Q(2)   +27*Q(3)   -  Q(4));
fx(nx-2)=-b2* (      Q(nx-1)-27*Q(nx-2)+27*Q(nx-3)-  Q(nx-4));
fx(nx-1)=-b2* (-22  *Q(nx-1)+17*Q(nx-2)+9 *Q(nx-3)-5*Q(nx-4)+Q(nx-5));
fx(nx)  =-b2* (-23  *Q(nx)  +21*Q(nx-1)+3 *Q(nx-2)-  Q(nx-3));

for i=3:nx-3
    QQ(i) = a1*(x(i+1)+x(i))+a2*(x(i+2)+x(i-1))+a3*(x(i+3)+x(i-2));
end

QQ0     =b1* (35 *x(1) - 35 *x(2)   + 21 *x(3)   - 5 *x(4));
QQ(1)   =b1* (5  *x(1) + 15 *x(2)   - 5  *x(3)   +    x(4));
QQ(2)   =b1* (-   x(1) + 9  *x(2)   + 9  *x(3)   -    x(4));
QQ(nx-2)=b1* (-   x(nx)+ 9  *x(nx-1)+ 9  *x(nx-2)-    x(nx-3));
QQ(nx-1)=b1* (5  *x(nx)+ 15 *x(nx-1)- 5  *x(nx-2)+    x(nx-3));
QQ(nx)  =b1* (35 *x(nx)- 35 *x(nx-1)+ 21 *x(nx-2)- 5 *x(nx-3));

%for i=4:nx-3
%    sx(i) = (c1*(QQ(i)-QQ(i-1))+c2*(QQ(i+1)-QQ(i-2))+c3*(QQ(i+2)-QQ(i-3)));
%end
%
%sx(1)   = b2* (-23  *QQ0     +21*QQ(1)   +3 *QQ(2)   -  QQ(3));
%sx(2)   = b2* (-22  *QQ(1)   +17*QQ(2)   +9 *QQ(3)   -5*QQ(4)   +QQ(5));
%sx(3)   = b2* (      QQ(1)   -27*QQ(2)   +27*QQ(3)   -  QQ(4));
%sx(nx-2)=-b2* (      QQ(nx-1)-27*QQ(nx-2)+27*QQ(nx-3)-  QQ(nx-4));
%sx(nx-1)=-b2* (-22  *QQ(nx-1)+17*QQ(nx-2)+9 *QQ(nx-3)-5*QQ(nx-4)+QQ(nx-5));
%sx(nx)  =-b2* (-23  *QQ(nx)  +21*QQ(nx-1)+3 *QQ(nx-2)-  QQ(nx-3));

for i=3:nx-1
    sx(i) =(d1*(QQ(i)-QQ(i-1))+d2*(x(i+1)-x(i-1))+d3*(QQ(i+1)-QQ(i-2)));
end

sx(1) = b2* (-23 *QQ0    +21*QQ(1)   +3*QQ(2)    -QQ(3));
sx(2) = (d1*(QQ(2)-QQ(1))+d2*(x(3)-x(1))+d3*(QQ(3)-QQ0));
sx(nx)=-b2* (-23 *QQ(nx) +21*QQ(nx-1)+3*QQ(nx-2) -QQ(nx-3));

dfdx=fx./sx;
%% exception for the derivation at the boundary
dfdx(1)=(f(2)-f(1))/(x(2)-x(1));
dfdx(nx)=(f(nx)-f(nx-1))/(x(nx)-x(nx-1));



dfdx=dfdx';
end
