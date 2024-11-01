%Объявляем глобальные переменные характерных размеров робота и констант
global R_l R_r  VM OQ cos120 sin120 cos240 sin240
%Вычисляем константы
cos120 = cosd(120);
sin120 = sind(120);
cos240 = cosd(240);
sin240 = sind(240);
%Задаём размеры робота [мм]
F = 270; %Длина стороны треугольника основания
f = 80; %Длина стороны треугольника платформы
R_l = 170;  %Длина рычагов (по осям)
R_r = 320;  %Длина штанг (по осям)
%Вычисляем радиусы вписанных окружностей
OQ = F*sqrt(3)/6; %Радиус окружности осей шарниров
VM = f*sqrt(3)/6; %Радиус окружности осей рычагов
%Задаём желаемые координаты центра подвижной платформы [мм]
Theta1 = 115; %Координата по оси X
Theta2 = 115; %Координата по оси Y
Theta3 = 115; %Координата по оси Z
%Вычисляем углы, на которые нужно повернуть рычаги
[L1, L2, L3, V] = PZK(Theta1, Theta2, Theta3)