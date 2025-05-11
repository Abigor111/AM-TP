%---> NRK2 Método de Runge-Kutta de ordem 2 para resolução numérica de um sistema de EDO/PVI
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
function [t,u,v] = NRK2SED(~,f,g,a,b,n,u0,v0)
h=(b-a)/n;                                  %Calcula o passo h com base no intervalo e número de subintervalos
t=a:h:b;                                    %Cria o vetor de pontos t no intervalo [a,b]
u=zeros(1,n+1);                             %Inicializa vetor u com zeros
v=zeros(1,n+1);                             %Inicializa vetor v com zeros
u(1)=u0;                                    %Define a condição inicial u(a)=u0
v(1)=v0;                                    %Define a condição inicial v(a)=v0

for i=1:n
    k1u=h*f(t(i),u(i),v(i));                %Calcula k1 para u usando f
    k1v=h*g(t(i),u(i),v(i));                %Calcula k1 para v usando g
    k2u=h*f(t(i+1),u(i)+k1u,v(i)+k1v);      %Calcula k2 para u usando valores intermediários
    k2v=h*g(t(i+1),u(i)+k1u,v(i)+k1v);      %Calcula k2 para v usando valores intermediários
    u(i+1)=u(i)+(k1u+k2u)/2;                %Atualiza u usando a média ponderada de k1 e k2
    v(i+1)=v(i)+(k1v+k2v)/2;                %Atualiza v usando a média ponderada de k1 e k2
end

end
