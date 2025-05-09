% ---> NEULER Método de Euler para resolução numérica de EDO/PVI
%INPUT:
% f - função da EDO y'=f(t,y)
% [a,b] - intervalo de valores da variável independente t
% n - número de subintervalos ou iterações do método
% y0 - aproximação inicial y(a)=y0
%OUTPUT:
% t - vetor do intervalo [a,b] discretizado
% y - vetor das soluções aproximadas do PVI em cada um dos t(i)
%     com t(i+1)=t(i)+h;
%AUTORES:
% Igor Carvalheira a2024128677@isec.pt
% Lucas Pantarotto a2024143625@isec.pt
% Rafael Carvalho a2024143302@isec.pt
function [t,y] = NEuler(f,a,b,n,y0)
    h = (b-a)/n;            % Calcula o tamanho do passo (h) com base no intervalo [a,b] e número de subintervalos n
    t = a:h:b;              % Cria o vetor t com os pontos discretizados do intervalo [a,b]
    y = zeros(1,n+1);        % Inicializa o vetor y com n elementos (devia ser n+1 para ser geral)
    y(1) = y0;              % Define a condição inicial: y(a) = y0
    for i = 1:n
        y(i+1) = y(i)+h*f(t(i),y(i));  % Aplica o método de Euler: aproxima o valor de y(i+1) com base na inclinação local
    end
    end