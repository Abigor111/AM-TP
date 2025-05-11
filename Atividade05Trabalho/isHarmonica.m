%---> isHarmonica Verifica se uma função é harmónica
%
%INPUT:
% f - função simbólica de duas variáveis (x, y)
%
%OUTPUT:
% YN - valor lógico (1 se a função é harmónica, 0 caso contrário)
%
%AUTORES:
% Igor Carvalheira    a2024128677@isec.pt
% Lucas Pantarotto    a2024143625@isec.pt
% Rafael Carvalho     a2024143302@isec.pt
function YN = isHarmonica(~,f)
    % Define as variáveis simbólicas x e y
    syms z(x,y)
    % Verifica se a soma das derivadas parciais de segunda ordem é nula
    if simplify(diff(f,x,2) + diff(f,y,2)) == 0
        YN = 1;  % É harmónica
    else
        YN = 0;  % Não é harmónica
    end           
end
