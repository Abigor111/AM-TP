function [t,y] = ODEfun(f,a, b, y0)
%ODE45 

% metodo de resolução numerica de Equações diferencias de 1a ordem
%   [a,b] - intervalo de valores da variável independente t
%   y0 - aproximação inicial y(a)=y0
%   y'=f(t,y), t=[a,b], y(a)=y0

%INPUT:
%f - função introduzida pelo utilizador
% a- menor valor do intervalo
%b - valor superior do intervalo
% y0 - valor de iteração inicial

%output: 
%dydt - intervalo de derivação em cada unidade de tempo dentro do intervalo
%t;

tspan = [a,b ];
[t, y] = ode45(f, tspan, y0);



end

