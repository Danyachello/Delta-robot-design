function []= InputConstant()
%Задаём параметры механизма
global f F R_l R_r cos120 sin120 cos240 sin240 VM OQ
f = 110;    %Длина стороны треугольника платформы
F = 270;    %Длина стороны треугольника основания
R_l = 170;  %Длина рычагов
R_r = 320;  %Длина штанг
%Расчёт констант
cos120 = cosd(120);
sin120 = sind(120);
cos240 = cosd(240);
sin240 = sind(240);

VM = f*sqrt(3)/6;
OQ = F*sqrt(3)/6;
end