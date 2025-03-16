%Task 3 Laborator: Compuneti functii cu un singur minim local (de ex: o functie care reprezinta o sfera, un paraboloid) 
%cu functii trigonometrice si analizati-le minimele.

% %raza si intervalul pe care vrem sa construim graficul
% r = 1;
% t = linspace(0, 2*pi, 50);
% p = linspace(0, pi, 50);
% 
% [t, p] = meshgrid(t, p);
% x = r * cos(t) .* sin(p);
% y = r * sin(t) .* sin(p);
% z = r * cos(p);
% 
% figure;
% surf(x,y,z);
% axis equal;
% xlabel('X'); ylabel('Y'); zlabel('Z');
% title('Sfera reprezentata trigonometric');
% 
% %ecuatia pentru o sfera de raza r este r.^2 = x.^2+y.^2+z.^2
% f = sqrt(r.^2-x.^2-y.^2);
% minime_locale = islocalmin(f,1) & islocalmin(f, 2);
% maxime_locale = islocalmax(f,1) & islocalmax(f, 2);
% 
% 
% %r este o constanta (l-am definit mai sus) deci nu are un minim si un maxim
% 
% x_min = x(minime_locale);
% y_min = y(minime_locale);
% % f de fapt este z
% z_min = z(minime_locale);
% f_min = f(minime_locale);
% 
% x_max = x(maxime_locale);
% y_max = y(maxime_locale);
% z_max = z(maxime_locale);
% f_max = f(maxime_locale);
% 
% hold on
% plot3(x_min, y_min, z_min, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
% plot3(x_max, y_max, z_max, 'go', 'MarkerSize', 10, 'LineWidth', 2);

%Functie obtinere functie
%1. Meshgrid
[x,y] = meshgrid(linspace(-5,5,100), linspace(-5,5,10));

%2. constanta A
A = 10;
B = 10;

%3. 
functie = x.^2./A.^2 + y.^2/B.^2;
%functie = y.^2./B.^2 - x.^2./A.^2;

%4. minime si maxime locale;
minime_locale = islocalmin(functie,1)& islocalmin(functie,2);
maxime_locale = islocalmax(functie,1) & islocalmax(functie,2);

%5. Extragem coordonatele pentru minime si maxime
x_min = x(minime_locale);
y_min = y(minime_locale);

f_min = functie(minime_locale);

x_max = x(maxime_locale);
y_max = y(maxime_locale);

f_max = functie(maxime_locale);

%6. Afisam graficul 3D
figure;
surf(x,y,functie,'EdgeColor','none');
colormap jet; shading interp; 
hold on;

%7. Marcam minime locale cu cercuri rosii
plot3(x_min, y_min, f_min, 'ro','MarkerSize',10,'LineWidth',2);

%8. Marcam maxime locale cu cercuri versi
plot3(x_max, y_max, f_max, 'go', 'MarkerSize', 10, 'LineWidt', 2);

%10. Setari pentru vizualizare
xlabel('X');
ylabel('Y');
legend('f(x,y)','Minime locale', 'Maxime locale', 'Location', 'best');
hold off;