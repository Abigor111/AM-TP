%---> NEulerM Método de Euler Melhorado (Heun) para sistemas de EDO/PVI
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

function [t,u,v] = NEulerMSED(~,f,g,a,b,n,u0,v0)
    h = (b - a) / n;                      % Passo
    t = a:h:b;                            % Vetor do tempo
    u = zeros(1, n+1);                    % Inicializa vetor u
    v = zeros(1, n+1);                    % Inicializa vetor v
    u(1) = u0;                            % Condição inicial u(a) = u0
    v(1) = v0;                            % Condição inicial v(a) = v0

    for i = 1:n
        % --- k1 ---
        k1u = f(t(i), u(i), v(i));
        k1v = g(t(i), u(i), v(i));

        % Previsão (Euler simples)
        u_pred = u(i) + h * k1u;
        v_pred = v(i) + h * k1v;

        % --- k2 ---
        k2u = f(t(i+1), u_pred, v_pred);
        k2v = g(t(i+1), u_pred, v_pred);

        % Correção (Heun)
        u(i+1) = u(i) + (h/2) * (k1u + k2u);
        v(i+1) = v(i) + (h/2) * (k1v + k2v);
    end

end
