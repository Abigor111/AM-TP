% ---> NRK2 Método de Runge-Kutta de ordem 2 para resolução numérica de EDO/PVI
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
function [t,y] = NRK2(f,a,b,n,y0)

    h=(b-a)/n;                         %Calcula o passo h com base no intervalo e no número de subintervalos
    t=a:h:b;                           %Cria o vetor de pontos t no intervalo [a,b]
    y=zeros(1,n+1);                    %Inicializa o vetor y com zeros
    y(1)=y0;                           %Define a condição inicial y(a)=y0

    for i=1:n
        k1=h*f(t(i),y(i));             %Calcula k1 com base em t(i) e y(i)
        k2=h*f(t(i+1),y(i)+k1);        %Calcula k2 usando t(i+1) e y(i)+k1
        y(i+1)=y(i)+(k1+k2)/2;         %Atualiza y(i+1) pela média dos incrementos k1 e k2
    end
    end