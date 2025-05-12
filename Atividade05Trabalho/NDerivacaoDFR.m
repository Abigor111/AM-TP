function [x,y,dydx] = NDerivacaoDFR(~, f,a,b,h,y)
    x = a:h:b;
    n = length(x);
    if nargin == 4
        y = f(x);
    end
    dydx = zeros(1,n);
    dydx(1) = (y(2) - y(1)) / h;           % Ainda é progressiva no primeiro ponto!
    for i = 2:n
        dydx(i) = (y(i) - y(i-1)) / h;     % Diferença regressiva
    end
end
