%Pentru a rezolva sistemul 
% 2x + y − z = 8
% 0.5y + 0.5z = 1
% −z = 1

% Definim matricea coeficienților și vectorul termenilor liberi
% A = [2 1 -1; -3 -1 2; -2 1 2];
% b = [8; -11; -3];

% Apelăm funcția pentru rezolvarea sistemului
% x = gauss_elimination(A, b);

% Afișăm rezultatul
% disp('Soluția sistemului Ax = b este:');
% disp(x);

%Eliminare Gauss
% function x = gauss_elimination(A, b);
%     % Dimensiunea sistemului
%     n = length(b);
% 
%     % Matricea extinsă
%     A = [A b];
% 
%     % Eliminare Gauss cu pivotare parțială
%     for k = 1:n-1;
%         % Pivotare parțială: alegem cel mai mare pivot din coloana curentă
%         [~, pivotRow] = max(abs(A(k:n, k)));
%         pivotRow = pivotRow + k - 1;
% 
%         % Schimbăm liniile dacă pivotul nu este deja pe poziția diagonală
%         if pivotRow ~= k
%             A([k, pivotRow], :) = A([pivotRow, k], :);
%         end
% 
%         % Eliminare pentru a obține forma triunghiulară
%         for i = k+1:n
%             m = A(i, k) / A(k, k); % Coeficient de eliminare
%             A(i, k:end) = A(i, k:end) - m * A(k, k:end);
%         end
%     end
% 
%     % Substituție inversă pentru obținerea soluției
%     x = zeros(n, 1);
%     for i = n:-1:1
%         x(i) = (A(i, end) - A(i, i+1:n) * x(i+1:n)) / A(i, i);
%     end
% end

%Task 1
%Pentru a rezolva sistemul 
% 3x + 2y − z = 1
% 2x − 2y + 4z = −2
% −x + 0.5y − z = 0

% Definim matricea coeficienților și vectorul termenilor liberi
A = [3 2 -1; 2 -2 4; -1 0.5 -1];
b = [1; -2; 0];

% Apelăm funcția pentru rezolvarea sistemului
x = gauss_elimination(A, b);

% Afișăm rezultatul
disp('Soluția sistemului Ax = b este:');
disp(x);
%Eliminare Gauss
function x = gauss_elimination(A, b);
   % Dimensiunea sistemului
   n = length(b);
   % Matricea extinsă
   A = [A b];
   % Eliminare Gauss cu pivotare parțială
   for k = 1:n-1;
       % Pivotare parțială: alegem cel mai mare pivot din coloana curentă
       [~, pivotRow] = max(abs(A(k:n, k)));
       pivotRow = pivotRow + k - 1;
       % Schimbăm liniile dacă pivotul nu este deja pe poziția diagonală
       if pivotRow ~= k
           A([k, pivotRow], :) = A([pivotRow, k], :);
       end
       % Eliminare pentru a obține forma triunghiulară
       for i = k+1:n
           m = A(i, k) / A(k, k); % Coeficient de eliminare
           A(i, k:end) = A(i, k:end) - m * A(k, k:end);
       end
   end
   % Substituție inversă pentru obținerea soluției
   x = zeros(n, 1);
   for i = n:-1:1
       x(i) = (A(i, end) - A(i, i+1:n) * x(i+1:n)) / A(i, i);
   end
end



