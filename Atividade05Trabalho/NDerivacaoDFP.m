% ---> NDERIVACAODFP Derivação Numérica por Diferenças Finitas Progressivas (DFP)
%
%INPUT:
% f   - função a derivar (caso y não seja fornecido)
% [a,b] - intervalo da variável independente x
% h   - passo da discretização
% y   - (opcional) vetor com os valores da função f(x) nos pontos x
%
%OUTPUT:
% x   - vetor do intervalo [a,b] discretizado
% y   - vetor dos valores da função (calculado se não fornecido)
% dydx - vetor das derivadas aproximadas de f em cada ponto x(i)
%        usando diferença finita progressiva: f'(x) ≈ (f(x+h)-f(x))/h
%
%AUTORES:
% Igor Carvalheira    a2024128677@isec.pt
% Lucas Pantarotto    a2024143625@isec.pt
% Rafael Carvalho     a2024143302@isec.pt

function [x,y,dydx] = NDerivacaoDFP(~,f,a,b,h,y)
    x = a:h:b;                           % Gera o vetor de pontos no intervalo [a,b]
    n = length(x);                       % Número de pontos
    if nargin == 5                       % Se y não for fornecido, calcula-se f(x)
        y = f(x);
    end
    dydx = zeros(1,n);                   % Inicializa vetor para derivadas
    for k = 1:n-1
        dydx(k) = (y(k+1)-y(k))/h;       % Diferença finita progressiva
    end
    dydx(n) = (y(n)-y(n-1))/h;           % Diferença finita regressiva no último ponto
end
