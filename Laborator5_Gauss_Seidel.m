clc; clear;

% % Definim matricea coeficienților și vectorul termenilor liberi
% A = [4 1 1; 2 5 2; 1 2 3];
% b = [7; 3; 5];
% 
% % Setăm toleranța și numărul maxim de iterații
% tol = 1e-10;
% max_iter = 1000;
% 
% % Apelăm funcția Gauss-Seidel
% x = gauss_seidel(A, b, tol, max_iter);
% 
% % Afișăm soluția
% disp('Soluția sistemului Ax = b este:')
% disp(x)
% 
% % Verificăm soluția: calculăm A * x și comparăm cu b
% b_calc = A * x;
% disp('Verificare: A * x trebuie să fie aproximativ egal cu b:')
% disp(b_calc)
% disp('b original:')
% disp(b)
% 
% % Calculăm eroarea
% error = norm(b - b_calc);
% disp(['Eroare între Ax și b: ', num2str(error)])
% 
% function x = gauss_seidel(A, b, tol, max_iter)
% % Funcție pentru rezolvarea sistemului Ax = b prin metoda Gauss-Seidel
% % Intrare:
% % A - matricea coeficienților (n x n)
% % b - vectorul termenilor liberi (n x 1)
% % tol - toleranța pentru criteriul de oprire
% % max_iter - numărul maxim de iterații
% % Ieșire:
% % x - soluția sistemului (n x 1)
% 
% n = length(b);
% x = zeros(n, 1); % Inițializare
% x_old = x;
% 
% for k = 1:max_iter
%     for i = 1:n
%         % Calculăm suma folosind noile valori pe loc
%         sum1 = A(i, 1:i-1) * x(1:i-1);
%         sum2 = A(i, i+1:n) * x_old(i+1:n);
%         x(i) = (b(i) - sum1 - sum2) / A(i, i);
%     end
% 
%     % Verificăm criteriul de oprire
%     if norm(x - x_old, inf) < tol
%         break;
%     end
% 
%     x_old = x; % Actualizăm valorile anterioare
% end
% end

% % Definim matricea coeficienților și vectorul termenilor liberi
% A = [4 1 1; 2 5 2; 1 2 3];
% b = [7; 3; 5];
% 
% % Setăm precizia și numărul maxim de iterații
% eps_sis = 1e-10;
% Nmax = 5000;
% X0 = zeros(length(b), 1);
% 
% % Apelăm funcția Gauss-Seidel
% [x, num_iter] = gauss_seidel(A, b, eps_sis, Nmax, X0);
% 
% % Afișăm soluția
% disp('Soluția sistemului Ax = b este:')
% disp(x)
% 
% % Verificăm soluția: calculăm A * x și comparăm cu b
% b_calc = A * x;
% disp('Verificare: A * x trebuie să fie aproximativ egal cu b:')
% disp(b_calc)
% disp('b original:')
% disp(b)
% 
% % Calculăm eroarea
% error = norm(b - b_calc);
% disp(['Eroare între Ax și b: ', num2str(error)])
% disp(['Număr de iterații necesare: ', num2str(num_iter)])
% 
% function [x, num_iter] = gauss_seidel(A, b, eps_sis, Nmax, X0)
% % Funcție pentru rezolvarea sistemului Ax = b prin metoda Gauss-Seidel
% % Intrare:
% % A - matricea coeficienților (n x n)
% % b - vectorul termenilor liberi (n x 1)
% % eps_sis - precizia pentru criteriul de oprire
% % Nmax - numărul maxim de iterații
% % X0 - iteratia inițială
% % Ieșire:
% % x - soluția sistemului (n x 1)
% % num_iter - numărul de iterații necesar
% 
% n = length(b);
% x = X0;
% num_iter = 0;
% 
% for k = 1:Nmax
%     x_old = x;
%     for i = 1:n
%         sum1 = A(i, 1:i-1) * x(1:i-1);
%         sum2 = A(i, i+1:n) * x_old(i+1:n);
%         x(i) = (b(i) - sum1 - sum2) / A(i, i);
%     end
% 
%     if norm(x - x_old, inf) < eps_sis
%         num_iter = k;
%         return;
%     end
% end
% 
% num_iter = Nmax;
% end

% % Generăm date de test
% [nlinii, ncoloane] = deal(3, 3);
% [A, X_exact, b] = genereaza_date_test(nlinii, ncoloane);
% 
% % Setăm precizia și numărul maxim de iterații
% eps_sis = 1e-10;
% Nmax = 5000;
% X0 = zeros(length(b), 1);
% 
% % Apelăm funcția Gauss-Seidel
% [x, num_iter] = gauss_seidel(A, b, eps_sis, Nmax, X0);
% 
% % Afișăm soluția
% disp('Soluția sistemului Ax = b este:')
% disp(x)
% 
% % Verificăm soluția: calculăm A * x și comparăm cu b
% b_calc = A * x;
% disp('Verificare: A * x trebuie să fie aproximativ egal cu b:')
% disp(b_calc)
% disp('b original:')
% disp(b)
% 
% % Calculăm eroarea
% error = norm(b - b_calc);
% disp(['Eroare între Ax și b: ', num2str(error)])
% disp(['Număr de iterații necesare: ', num2str(num_iter)])
% 
% function [x, num_iter] = gauss_seidel(A, b, eps_sis, Nmax, X0)
% % Funcție pentru rezolvarea sistemului Ax = b prin metoda Gauss-Seidel
% % Intrare:
% % A - matricea coeficienților (n x n)
% % b - vectorul termenilor liberi (n x 1)
% % eps_sis - precizia pentru criteriul de oprire
% % Nmax - numărul maxim de iterații
% % X0 - iteratia inițială
% % Ieșire:
% % x - soluția sistemului (n x 1)
% % num_iter - numărul de iterații necesar
% 
% n = length(b);
% x = X0;
% num_iter = 0;
% 
% for k = 1:Nmax
%     x_old = x;
%     for i = 1:n
%         sum1 = A(i, 1:i-1) * x(1:i-1);
%         sum2 = A(i, i+1:n) * x_old(i+1:n);
%         x(i) = (b(i) - sum1 - sum2) / A(i, i);
%     end
% 
%     if norm(x - x_old, inf) < eps_sis
%         num_iter = k;
%         return;
%     end
% end
% 
% num_iter = Nmax;
% end
% 
% function [A, X, b] = genereaza_date_test(nlinii, ncoloane)
% % Funcție pentru generarea datelor de test
% % A - matrice cu diagonala dominantă
% % X - vector soluție
% % b - vectorul termenilor liberi
% 
% A = rand(nlinii, ncoloane) * 10;
% for i = 1:nlinii
%     A(i, i) = sum(abs(A(i, :))) + rand() * 10; % Asigură diagonală dominantă
% end
% X = rand(ncoloane, 1) * 10;
% b = A * X;
% end

% % Generăm și testăm cel puțin 3 perechi (A, B)
% n_tests = 3;
% eps_sis = 1e-10;
% Nmax = 5000;
% 
% for test = 1:n_tests
%     [A, X_exact, b] = genereaza_date_test(3, 3);
%     X0 = zeros(length(b), 1);
% 
%     % Apelăm funcția Gauss-Seidel
%     [X_implementat, num_iter] = gauss_seidel(A, b, eps_sis, Nmax, X0);
% 
%     % Calculăm eroarea relativă
%     eroare_relativa = norm(X_implementat - X_exact) / norm(X_exact);
% 
%     % Afișăm rezultatele
%     fprintf('Test %d:\n', test);
%     disp('Soluția implementată:');
%     disp(X_implementat);
%     disp('Soluția exactă:');
%     disp(X_exact);
%     disp(['Eroare relativă: ', num2str(eroare_relativa)]);
%     disp(['Număr de iterații necesare: ', num2str(num_iter)]);
% 
%     % Comparăm rezultatul obținut cu eps_sis
%     if eroare_relativa < eps_sis
%         disp('Eroarea relativă este sub pragul de precizie.');
%     else
%         disp('Eroarea relativă este peste pragul de precizie.');
%     end
%     disp('--------------------------------');
% end
% 
% function [x, num_iter] = gauss_seidel(A, b, eps_sis, Nmax, X0)
% % Funcție pentru rezolvarea sistemului Ax = b prin metoda Gauss-Seidel
% % Intrare:
% % A - matricea coeficienților (n x n)
% % b - vectorul termenilor liberi (n x 1)
% % eps_sis - precizia pentru criteriul de oprire
% % Nmax - numărul maxim de iterații
% % X0 - iteratia inițială
% % Ieșire:
% % x - soluția sistemului (n x 1)
% % num_iter - numărul de iterații necesar
% 
% n = length(b);
% x = X0;
% num_iter = 0;
% 
% for k = 1:Nmax
%     x_old = x;
%     for i = 1:n
%         sum1 = A(i, 1:i-1) * x(1:i-1);
%         sum2 = A(i, i+1:n) * x_old(i+1:n);
%         x(i) = (b(i) - sum1 - sum2) / A(i, i);
%     end
% 
%     if norm(x - x_old, inf) < eps_sis
%         num_iter = k;
%         return;
%     end
% end
% 
% num_iter = Nmax;
% end
% 
% function [A, X, b] = genereaza_date_test(nlinii, ncoloane)
% % Funcție pentru generarea datelor de test
% % A - matrice cu diagonala dominantă
% % X - vector soluție
% % b - vectorul termenilor liberi
% 
% A = rand(nlinii, ncoloane) * 10;
% for i = 1:nlinii
%     A(i, i) = sum(abs(A(i, :))) + rand() * 10; % Asigură diagonală dominantă
% end
% X = rand(ncoloane, 1) * 10;
% b = A * X;
% end


% Generăm și testăm variația numărului de iterații în funcție de precizie
n_tests = 3;  % Numărul de teste
eps_sis_values = logspace(-1, -16, 16);  % Intervalul de valori pentru eps_sis (de la 10^-1 la 10^-16)
Nmax = 5000;  % Numărul maxim de iterații

% Stocăm numărul de iterații și precizia pentru fiecare test
num_iterations = zeros(length(eps_sis_values), 1);
precizii = zeros(length(eps_sis_values), 1);

% Pentru fiecare test, vom rezolva sistemul cu diverse valori ale lui eps_sis
for test = 1:n_tests
    [A, X_exact, b] = genereaza_date_test(3, 3);  % Generăm datele de test
    X0 = zeros(length(b), 1);  % Soluția inițială (0)

    % Pentru fiecare valoare de precizie (eps_sis), testăm Gauss-Seidel
    for i = 1:length(eps_sis_values)
        eps_sis = eps_sis_values(i);
        
        % Apelăm funcția Gauss-Seidel
        [X_implementat, num_iter] = gauss_seidel(A, b, eps_sis, Nmax, X0);
        
        % Stocăm numărul de iterații și precizia
        num_iterations(i) = num_iter;
        precizii(i) = eps_sis;
    end
    
end

% Generăm graficul
figure;
loglog(precizii, num_iterations, 'o-', 'LineWidth', 2);
xlabel('Precizia (\epsilon_{sis})');
ylabel('Număr de iterații');
title('Variația numărului de iterații în funcție de precizie');
grid on;

% Gauss-Seidel Method Function
function [x, num_iter] = gauss_seidel(A, b, eps_sis, Nmax, X0)
    % Funcție pentru rezolvarea sistemului Ax = b prin metoda Gauss-Seidel
    % Intrare:
    % A - matricea coeficienților (n x n)
    % b - vectorul termenilor liberi (n x 1)
    % eps_sis - precizia pentru criteriul de oprire
    % Nmax - numărul maxim de iterații
    % X0 - iteratia inițială
    % Ieșire:
    % x - soluția sistemului (n x 1)
    % num_iter - numărul de iterații necesar
    
    n = length(b);
    x = X0;  % Soluția inițială
    num_iter = 0;

    for k = 1:Nmax
        x_old = x;  % Salvăm soluția anterioară
        for i = 1:n
            % Calculăm suma pentru Gauss-Seidel
            sum1 = A(i, 1:i-1) * x(1:i-1);  % Partea stângă a ecuației
            sum2 = A(i, i+1:n) * x_old(i+1:n);  % Partea dreaptă (folosind soluția anterioară)
            
            % Actualizăm valoarea x(i)
            x(i) = (b(i) - sum1 - sum2) / A(i, i);
        end

        % Verificăm criteriul de oprire
        if norm(x - x_old, inf) < eps_sis
            num_iter = k;
            return;
        end
    end

    num_iter = Nmax;  % Dacă nu s-a atins precizia, returnăm numărul maxim de iterații
end

%Generare date test
function [A, X, b] = genereaza_date_test(nlinii, ncoloane)
    % Funcție pentru generarea datelor de test
    % A - matrice cu diagonala dominantă
    % X - vector soluție
    % b - vectorul termenilor liberi
    
    A = rand(nlinii, ncoloane) * 10;  % Generăm matricea A cu valori aleatorii
    for i = 1:nlinii
        A(i, i) = sum(abs(A(i, :))) + rand() * 10; % Asigură diagonală dominantă
    end
    
    X = rand(ncoloane, 1) * 10;  % Generăm vectorul soluției X aleatoriu
    b = A * X;  % Calculăm vectorul b
end
