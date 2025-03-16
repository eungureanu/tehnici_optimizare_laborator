%Exercițiu 1: Calculați erorile absolute și relative pentru următoarele exemple. x1 este valoarea exactă și a1 o aproximație a sa, unde i = ̅1̅,̅4̅:
% Pentru a1 =3.14; x1=3.141592
% Pentru a2 =999996; x2=1000000
% Pentru a3 =0.000009; x3=0.000012
% Pentru a4 =1.00345; x4=1.000145

% a) Pentru toate exemplele de mai sus să se calculeze Δai, |Δai| și δai
% b) Care dintre măsurătorile rezultate este mai bună?
% c) Când se realizează o aproximare prin lipsă și când prin adaos?

fprintf('Exercitiul 1');
format long %afisarea numerelor cu 15 zecimale

a_values = [3.14, 999996, 0.000009, 1.00345];
x_values = [3.141592, 1000000, 0.000012, 1.000145];

for i = 1:length(a_values)
    a = a_values(i);
    x = x_values(i);
    delta_a=x-a

    if delta_a == 0
        fprintf('delta_a este 0\n')
        continue;
    end

    eroare_absoluta=abs(x-a);
    eroare_relativa=eroare_absoluta/abs(a);

    % pentru afisarea numerelor cu un numar specificat de zecimale
    fprintf('eroare_absoluta = %.8f\n',eroare_absoluta)
    fprintf('eroare_relativa = %.8f\n',eroare_relativa)
end

% Exerciţiul 2:
fprintf('Exercitiul 2');
%completați programul pentru a afișa un mesaj prin care se menționează dacă a îl aproximează pe x prin lipsă sau prin adaos.
%utilizați funcția input() pentru citirea de la tastatură a valorilor lui a și x.

