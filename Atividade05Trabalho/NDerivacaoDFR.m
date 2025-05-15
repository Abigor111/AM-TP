%---> NDerivacaoDFR Deriva��o num�rica usando diferen�as finitas regressivas
%
%INPUT:
% f  - fun��o a derivar (ou vetor de valores se y for omitido)
% [a,b] - intervalo onde se pretende calcular a derivada
% h  - passo de discretiza��o
% y  - (opcional) vetor de valores j� calculados da fun��o f(x)
%
%OUTPUT:
% x  - vetor de pontos no intervalo [a,b] com passo h
% y  - vetor dos valores de f(x) nos pontos x(i)
% dydx - vetor das derivadas aproximadas em cada x(i)
%
%AUTORES:
% Igor Carvalheira    a2024128677@isec.pt
% Lucas Pantarotto    a2024143625@isec.pt
% Rafael Carvalho     a2024143302@isec.pt
function [x,y,dydx] = NDerivacaoDFR(~,f,a,b,h,y)
    x = a:h:b;              % Cria vetor de pontos entre a e b com passo h
    n = length(x);          % N�mero total de pontos
    if nargin == 5
        y = f(x);           % Calcula os valores da fun��o se y n�o for dado
    end
    dydx = zeros(1,n);      % Inicializa o vetor da derivada com zeros
    dydx(1) = (y(2)-y(1))/h; % Estimativa inicial usando diferen�a progressiva
    for k = 2:n
        dydx(k) = (y(k)-y(k-1))/h; % Diferen�a regressiva para os restantes pontos
    end
end