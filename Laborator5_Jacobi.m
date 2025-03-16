% % Definim matricea coeficienților și vectorul termenilor liberi
% A = [10 1 1; 1 10 -1; 1 -1 10];
% b = [6; 3; 6];
% 
% % Setăm toleranța și numărul maxim de iterații
% eps_sis = 1e-6;
% Nmax = 5000;
% 
% % Inițializăm soluția inițială
% X0 = zeros(length(b), 1);
% 
% % Apelăm funcția Jacobi
% [x, num_iter] = jacobi(A, b, eps_sis, Nmax, X0);
% 
% % Afișăm soluția
% fprintf('Soluția sistemului Ax = b este:\n');
% disp(x);
% 
% % Afișăm numărul de iterații necesar
% fprintf('Numărul de iterații necesar: %d\n', num_iter);
% 
% % Verificăm soluția: calculăm A * x și comparăm cu b
% b_calc = A * x;
% fprintf('Verificare: A * x trebuie să fie aproximativ egal cu b:\n');
% disp(b_calc);
% fprintf('b original:\n');
% disp(b);
% 
% function [x, num_iter] = jacobi(A, b, eps_sis, Nmax, X0)
%     % Funcție pentru rezolvarea sistemului Ax = b prin metoda iterativă Jacobi
%     % Intrare:
%     % A - matricea coeficienților (n x n)
%     % b - vectorul termenilor liberi (n x 1)
%     % eps_sis - toleranța pentru criteriul de oprire
%     % Nmax - numărul maxim de iterații admis
%     % X0 - soluția inițială (n x 1)
%     % Ieșire:
%     % x - soluția sistemului (n x 1)
%     % num_iter - numărul de iterații necesar
% 
%     n = length(b);
%     x = X0;
%     x_new = zeros(n, 1);
%     num_iter = 0;
% 
%     for k = 1:Nmax
%         for i = 1:n
%             sum_term = A(i, 1:i-1) * x(1:i-1) + A(i, i+1:n) * x(i+1:n);
%             x_new(i) = (b(i) - sum_term) / A(i, i);
%         end
% 
%         % Verificăm criteriul de oprire
%         if norm(x_new - x, inf) < eps_sis
%             num_iter = k;
%             x = x_new;
%             return;
%         end
% 
%         x = x_new;
%     end
% 
%     num_iter = Nmax;
% end

% % Definim matricea coeficienților și vectorul termenilor liberi
% A = [10 1 1; 1 10 -1; 1 -1 10];
% b = [6; 3; 6];
% 
% % Setăm toleranța și numărul maxim de iterații
% eps_sis = 1e-6;
% Nmax = 5000;
% 
% % Generăm date de test
% nlinii = 3;
% ncoloane = 3;
% [A, X, b] = genereaza_date_test(nlinii, ncoloane);
% 
% % Inițializăm soluția inițială
% X0 = zeros(length(b), 1);
% 
% % Apelăm funcția Jacobi
% [x, num_iter] = jacobi(A, b, eps_sis, Nmax, X0);
% 
% % Afișăm soluția
% fprintf('Soluția sistemului Ax = b este:\n');
% disp(x);
% 
% % Afișăm numărul de iterații necesar
% fprintf('Numărul de iterații necesar: %d\n', num_iter);
% 
% % Verificăm soluția: calculăm A * x și comparăm cu b
% b_calc = A * x;
% fprintf('Verificare: A * x trebuie să fie aproximativ egal cu b:\n');
% disp(b_calc);
% fprintf('b original:\n');
% disp(b);
% 
% function [x, num_iter] = jacobi(A, b, eps_sis, Nmax, X0)
%     n = length(b);
%     x = X0;
%     x_new = zeros(n, 1);
%     num_iter = 0;
% 
%     for k = 1:Nmax
%         for i = 1:n
%             sum_term = A(i, 1:i-1) * x(1:i-1) + A(i, i+1:n) * x(i+1:n);
%             x_new(i) = (b(i) - sum_term) / A(i, i);
%         end
% 
%         if norm(x_new - x, inf) < eps_sis
%             num_iter = k;
%             x = x_new;
%             return;
%         end
% 
%         x = x_new;
%     end
% 
%     num_iter = Nmax;
% end
% 
% function [A, X, b] = genereaza_date_test(nlinii, ncoloane)
%     A = rand(nlinii, ncoloane) * 10;
%     for i = 1:nlinii
%         A(i, i) = sum(abs(A(i, :))) + rand() * 10;
%     end
%     X = rand(ncoloane, 1) * 10;
%     b = A * X;
% end

% % Definim parametrii
% nlinii = 3;
% ncoloane = 3;
% eps_sis = 1e-6;
% Nmax = 5000;
% 
% % Generăm și testăm 3 perechi (A, B)
% for test_idx = 1:3
%     [A, X_exact, b] = genereaza_date_test(nlinii, ncoloane);
%     X0 = zeros(length(b), 1);
% 
%     % Apelăm funcția Jacobi
%     [X_implementat, num_iter] = jacobi(A, b, eps_sis, Nmax, X0);
% 
%     % Calculăm eroarea relativă
%     eroare_relativa = norm(X_implementat - X_exact);
% 
%     % Afișăm rezultatele
%     fprintf('Test %d:\n', test_idx);
%     fprintf('Soluția exactă:\n');
%     disp(X_exact);
%     fprintf('Soluția obținută prin Jacobi:\n');
%     disp(X_implementat);
%     fprintf('Numărul de iterații necesar: %d\n', num_iter);
%     fprintf('Eroare relativă: %.10f\n', eroare_relativa);
%     fprintf('Comparare cu eps_sis: %.10f\n\n', eps_sis);
% end
% 
% function [x, num_iter] = jacobi(A, b, eps_sis, Nmax, X0)
%     n = length(b);
%     x = X0;
%     x_new = zeros(n, 1);
%     num_iter = 0;
% 
%     for k = 1:Nmax
%         for i = 1:n
%             sum_term = A(i, 1:i-1) * x(1:i-1) + A(i, i+1:n) * x(i+1:n);
%             x_new(i) = (b(i) - sum_term) / A(i, i);
%         end
% 
%         if norm(x_new - x, inf) < eps_sis
%             num_iter = k;
%             x = x_new;
%             return;
%         end
% 
%         x = x_new;
%     end
% 
%     num_iter = Nmax;
% end
% 
% function [A, X, b] = genereaza_date_test(nlinii, ncoloane)
%     A = rand(nlinii, ncoloane) * 10;
%     for i = 1:nlinii
%         A(i, i) = sum(abs(A(i, :))) + rand() * 10;
%     end
%     X = rand(ncoloane, 1) * 10;
%     b = A * X;
% end

% Definim parametrii
nlinii = 3;
ncoloane = 3;
Nmax = 5000;

% Generăm matricea și vectorul b
[A, X_exact, b] = genereaza_date_test(nlinii, ncoloane);
X0 = zeros(length(b), 1);

% Definim un vector de precizii eps_sis
precizii = logspace(-1, -16, 16);
numar_iteratii = zeros(size(precizii));

% Rezolvăm sistemul pentru fiecare precizie
for i = 1:length(precizii)
    eps_sis = precizii(i);
    [~, num_iter] = jacobi(A, b, eps_sis, Nmax, X0);
    numar_iteratii(i) = num_iter;
end

% Reprezentăm grafic variația numărului de iterații în funcție de precizie
figure;
semilogx(precizii, numar_iteratii, '-o', 'LineWidth', 2);
grid on;
xlabel('Precizie (eps\_sis)');
ylabel('Numărul de iterații');
title('Variatia numărului de iterații în funcție de precizie');

function [x, num_iter] = jacobi(A, b, eps_sis, Nmax, X0)
    n = length(b);
    x = X0;
    x_new = zeros(n, 1);
    num_iter = 0;

    for k = 1:Nmax
        for i = 1:n
            sum_term = A(i, 1:i-1) * x(1:i-1) + A(i, i+1:n) * x(i+1:n);
            x_new(i) = (b(i) - sum_term) / A(i, i);
        end
        
        if norm(x_new - x, inf) < eps_sis
            num_iter = k;
            x = x_new;
            return;
        end
        
        x = x_new;
    end
    
    num_iter = Nmax;
end

function [A, X, b] = genereaza_date_test(nlinii, ncoloane)
    A = rand(nlinii, ncoloane) * 10;
    for i = 1:nlinii
        A(i, i) = sum(abs(A(i, :))) + rand() * 10;
    end
    X = rand(ncoloane, 1) * 10;
    b = A * X;
end
