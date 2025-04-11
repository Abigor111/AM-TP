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
    h=(b-a)/n;                         %Calcula o passo h com base no intervalo e no número de subintervalos
    t(1)=a;                            %Define o primeiro valor de t como o início do intervalo
    y(1)=y0;                           %Define a condição inicial y(a)=y0
    for i=1:n
        y(i+1)=y(i)+h*f(t(i),y(i));    %Aplica a fórmula de Euler para calcular y(i+1)
        t(i+1)=t(i)+h;                 %Atualiza o valor de t para o próximo passo
    end
    end