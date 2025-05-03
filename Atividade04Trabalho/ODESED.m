%---> FUNCODE Resolução de um sistema de EDOs/PVI utilizando o método ode45 do MATLAB
%
%INPUT:
% f - função vetorial do sistema: y' = f(t, y), onde y = [u; v]
% [a,b] - intervalo de valores da variável independente t
% n - número de subintervalos (define a discretização do vetor t)
% u0 - valor inicial da primeira variável (u(a) = u0)
% v0 - valor inicial da segunda variável (v(a) = v0)
%
%OUTPUT:
% y - matriz transposta com as soluções aproximadas do sistema nos pontos t(i)
%     y(1,:) contém os valores de u(t), y(2,:) os valores de v(t)
%
%AUTORES:
% Igor Carvalheira a2024128677@isec.pt
% Lucas Pantarotto a2024143625@isec.pt
% Rafael Carvalho a2024143302@isec.pt
function y = ODESED(~,f,a,b,n,u0,v0)
    h = (b-a)/n;                          % Calcula o passo de discretização h
    t = a:h:b;                            % Cria o vetor t com n+1 pontos
    y0 = [u0; v0];                        % Condição inicial vetorial
    [~, yODE45] = ode45(f, t, y0);        % Chama ode45 para o sistema
    y = yODE45.';                         % Transpõe para manter y(1,:) = u, y(2,:) = v
end
