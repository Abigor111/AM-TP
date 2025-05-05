function [t, x, y, z] = NEulerM3(f1, f2, f3, a, b, n, x0, y0, z0)
% NEulerM3 - Euler Melhorado (Heun) para sistemas de 3 EDOs
% INPUT:
%   f1, f2, f3 - funções derivadas para x, y, z
%   [a,b] - intervalo de tempo
%   n - número de passos
%   x0, y0, z0 - condições iniciais
% OUTPUT:
%   t - vetor do tempo
%   x, y, z - soluções aproximadas

    h = (b - a) / n;       % Passo
    t = a:h:b;             % Vetor do tempo

    % Inicializar vetores solução
    x = zeros(1, n+1); y = zeros(1, n+1); z = zeros(1, n+1);

    % Condições iniciais
    x(1) = x0; y(1) = y0; z(1) = z0;

    for i = 1:n
        % --- k1 ---
        k1x = f1(t(i), x(i), y(i), z(i));
        k1y = f2(t(i), x(i), y(i), z(i));
        k1z = f3(t(i), x(i), y(i), z(i));

        % Previsão (Euler simples)
        x_pred = x(i) + h * k1x;
        y_pred = y(i) + h * k1y;
        z_pred = z(i) + h * k1z;

        % --- k2 ---
        k2x = f1(t(i+1), x_pred, y_pred, z_pred);
        k2y = f2(t(i+1), x_pred, y_pred, z_pred);
        k2z = f3(t(i+1), x_pred, y_pred, z_pred);

        % Correção (Heun)
        x(i+1) = x(i) + (h/2) * (k1x + k2x);
        y(i+1) = y(i) + (h/2) * (k1y + k2y);
        z(i+1) = z(i) + (h/2) * (k1z + k2z);
    end
end
