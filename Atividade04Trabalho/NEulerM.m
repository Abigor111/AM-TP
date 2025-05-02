% --> NEULER Método de Euler Melhorado (Heun) para resolução numérica de EDO/PVI
%INPUT:
%   f - função da EDO y'=f(t,y)
%   [a,b] - intervalo de valores da variável independente t
%   n - número de subintervalos ou iterações do método
%   y0 - aproximação inicial y(a)=y0
%OUTPUT:
%   t - vetor do intervalo [a,b] discretizado
%   y - vetor das soluções aproximadas do PVI em cada um dos t(i)
%   t(i+1) = t(i)+h;
%AUTORES
% Igor Carvalheira a2024128677@isec.pt
% Lucas Pantarotto a2024143625@isec.pt
% Rafael Carvalho a2024143302@isec.pt
function [t, y] = NEulerM(f, a, b, n, y0)
    h = (b - a) / n;         % Passo
    t = a:h:b;               % Vetor de tempo
    y = zeros(1, n+1);       % Inicializa vetor das soluções
    y(1) = y0;               % Condição inicial

    for i = 1:n
        k1 = f(t(i), y(i));                          % Pendente no início
        y_pred = y(i) + h * k1;                      % Previsão com Euler
        k2 = f(t(i+1), y_pred);                      % Pendente na previsão
        y(i+1) = y(i) + (h/2) * (k1 + k2);            % Correção (média das pendentes)
    end
end