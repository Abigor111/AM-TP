%---> NEULER Método de Euler para resolução numérica de um sistema de EDO/PVI
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
function [t,u,v] = NEulerSED(~,f,g,a,b,n,u0,v0)
h=(b-a)/n;                             %Calcula o passo h com base no intervalo e número de subintervalos
t=a:h:b;                               %Cria o vetor de pontos t no intervalo [a,b]
u=zeros(1,n+1);                        %Inicializa vetor u com zeros
v=zeros(1,n+1);                        %Inicializa vetor v com zeros
u(1)=u0;                               %Define a condição inicial u(a)=u0
v(1)=v0;                               %Define a condição inicial v(a)=v0

for i=1:n
    u(i+1)=u(i)+h*f(t(i),u(i),v(i));   %Atualiza u usando a EDO f com método de Euler
    v(i+1)=v(i)+h*g(t(i),u(i),v(i));   %Atualiza v usando a EDO g com método de Euler
end

end
