%---> Método de Adams-Bashforth de ordem 2 para resolução de SED/PVI
%
%INPUT:
% f - função da 1ª EDO do sistema: u' = f(t,u,v)
% g - função da 2ª EDO do sistema: v' = g(t,u,v)
% [a,b] - intervalo de valores da variável independente t
% n - número de subintervalos ou iterações do método
% u0 - aproximação inicial da 1ª variável (u(a) = u0)
% v0 - aproximação inicial da 2ª variável (v(a) = v0)
%
%OUTPUT:
% t - vetor do intervalo [a,b] discretizado
% u - vetor das soluções aproximadas da 1ª equação nos pontos t(i)
% v - vetor das soluções aproximadas da 2ª equação nos pontos t(i)
%     com t(i+1) = t(i) + h;
%
%AUTORES:
% Igor Carvalheira    a2024128677@isec.pt
% Lucas Pantarotto    a2024143625@isec.pt
% Rafael Carvalho     a2024143302@isec.pt

function [t,u,v] = AdamsSED(~,f,g,a,b,n,u0,v0)
    h = (b-a)/n;                          % Calcula o passo h
    t = a:h:b;                            % Cria o vetor de tempo com n+1 pontos
    u = zeros(1,n+1);                     % Inicializa vetor da 1ª variável
    v = zeros(1,n+1);                     % Inicializa vetor da 2ª variável
    u(1) = u0;                            % Condição inicial u(a) = u0
    v(1) = v0;                            % Condição inicial v(a) = v0

    % Passo 1: usar método de Euler para calcular o segundo ponto
    u(2) = u(1) + h * f(t(1), u(1), v(1));
    v(2) = v(1) + h * g(t(1), u(1), v(1));

    % Passos seguintes: método de Adams-Bashforth de 2ª ordem
    for i = 2:n
        u(i+1) = u(i) + h/2 * (3*f(t(i),u(i),v(i)) - f(t(i-1),u(i-1),v(i-1)));
        v(i+1) = v(i) + h/2 * (3*g(t(i),u(i),v(i)) - g(t(i-1),u(i-1),v(i-1)));
    end
end
