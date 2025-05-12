function [x,y,dydx]=NDerivacaoDFCENT(f,a,b,h,y)
x=a:h:b;
n=length(x);
if nargin==4
    y=f(x);
end;
dydx=zeros(1,n);
dydx(1)=(y(3)-y(2))/h;
for i=2:n-1
    dydx(i)=(y(i+1)-y(i-1))/(2*h);
end;
dydx(n)=(y(n-2)-y(n-1))/h;
