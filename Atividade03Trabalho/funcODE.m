%FUNCODE Resolução de EDO/PVI utilizando o método ode45 do MATLAB
%INPUT:
% f - função da EDO y'=f(t,y)
% [a,b] - intervalo de valores da variável independente t
% n - número de subintervalos (define a discretização do vetor t)
% y0 - aproximação inicial y(a)=y0
%OUTPUT:
% y - vetor transposto com as soluções aproximadas do PVI nos pontos t(i)
%AUTORES:
% Igor Carvalheira a2024128677@isec.pt
% Lucas Pantarotto a2024143625@isec.pt
% Rafael Carvalho a2024143302@isec.pt
function y = funcODE(f,a,b,n,y0)
    h=(b-a)/n;                          %Calcula o passo de discretização h
    t=a:h:b;                            %Cria o vetor t no intervalo [a,b]
    [~,yODE45]=ode45(f,t,y0);          %Resolve a EDO com ode45
    y=yODE45.';                         %Transpõe o resultado para manter consistência com os outros métodos
    end