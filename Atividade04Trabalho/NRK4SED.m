%---> NRK4 Método de Runge-Kutta de ordem 4 para resolução numérica de um sistema de EDO/PVI
%
%INPUT:
% f - função da primeira EDO u'=f(t,u,v)
% g - função da segunda EDO v'=g(t,u,v)
% [a,b] - intervalo de valores da variável independente t
% n - número de subintervalos ou iterações do método
% u0 - valor inicial u(a)=u0
% v0 - valor inicial v(a)=v0
%
%OUTPUT:
% t - vetor do intervalo [a,b] discretizado
% u - vetor das soluções aproximadas de u(t)
% v - vetor das soluções aproximadas de v(t)
%     com t(i+1)=t(i)+h;
%
%AUTORES:
% Igor Carvalheira a2024128677@isec.pt
% Lucas Pantarotto a2024143625@isec.pt
% Rafael Carvalho a2024143302@isec.pt

function [t,u,v] = NRK4SED(~,f,g,a,b,n,u0,v0)
h = (b-a)/n;                                %Calcula o passo h
t = a:h:b;                                  %Vetor do tempo
u = zeros(1,n+1);                           %Inicializa vetor u
v = zeros(1,n+1);                           %Inicializa vetor v
u(1) = u0;                                  %Condição inicial para u
v(1) = v0;                                  %Condição inicial para v

for i = 1:n
    k1u = h*f(t(i),u(i),v(i));              %RK4: k1 para u
    k1v = h*g(t(i),u(i),v(i));              %RK4: k1 para v

    k2u = h*f(t(i)+h/2,u(i)+k1u/2,v(i)+k1v/2); %RK4: k2 para u
    k2v = h*g(t(i)+h/2,u(i)+k1u/2,v(i)+k1v/2); %RK4: k2 para v

    k3u = h*f(t(i)+h/2,u(i)+k2u/2,v(i)+k2v/2); %RK4: k3 para u
    k3v = h*g(t(i)+h/2,u(i)+k2u/2,v(i)+k2v/2); %RK4: k3 para v

    k4u = h*f(t(i)+h,u(i)+k3u,v(i)+k3v);       %RK4: k4 para u
    k4v = h*g(t(i)+h,u(i)+k3u,v(i)+k3v);       %RK4: k4 para v

    u(i+1) = u(i)+(1/6)*(k1u+2*k2u+2*k3u+k4u); %Atualiza u(i+1)
    v(i+1) = v(i)+(1/6)*(k1v+2*k2v+2*k3v+k4v); %Atualiza v(i+1)
end

end
