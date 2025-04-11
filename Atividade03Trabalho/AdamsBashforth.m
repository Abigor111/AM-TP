% ---> Método de Adams-Bashforth de ordem 2
%INPUT:
% f - função da EDO y'=f(t,y)
% [a,b] - intervalo de valores da variável independente t
% n - número de subintervalos ou iterações do método
% y0 - aproximação inicial y(a)=y0
%OUTPUT:
% t - vetor do intervalo [a,b] discretizado
% y - vetor das soluções aproximadas do PVI em cada t(i)
%     com t(i+1)=t(i)+h;
%AUTORES:
% Igor Carvalheira a2024128677@isec.pt
% Lucas Pantarotto a2024143625@isec.pt
% Rafael Carvalho a2024143302@isec.pt
function [t,y] = AdamsBashforth(f,a,b,n,y0)
    h=(b-a)/n;                        %Calcula o passo com base no intervalo e no número de subintervalos
    t=a:h:b;                          %Cria o vetor de pontos no intervalo [a,b]
    y=zeros(1,n+1);                   %Inicializa o vetor y com zeros
    y(1)=y0;                          %Define a condição inicial y(a)=y0
    y(2)=y(1)+h*f(t(1),y(1));         %Calcula o segundo valor usando o método de Euler
    for i=2:n
        y(i+1)=y(i)+h/2*(3*f(t(i),y(i))-f(t(i-1),y(i-1)));  %Fórmula de Adams-Bashforth de 2ª ordem
    end
    
    end
    