global vmax %Максимальная требуемая скорость исполнительного звена
global Jr mk %Массоинерционных характеристики
global amax %Максимальное требуемое ускорение рабочего органа
%Задаём начальное значение
%Добавляем глобальные переменные входных и выходных данных
global WZ_D WZ_H WZ_Z WZ_d WZ_h %Рабочая зона
global addSeg %Переменная, показывающая необходимость дополнительного
%сегмента в рабочей зоне
global dotDensity %Плотность точек

%Плотность точек анализа
dotDensity = 30;
%Параметры рабочей зоны
addSeg = 1; %Добавляемый к рабочей зоне сегмент. 0 - ничего; 1 - усечённый
%конус; %2 - часть сферы
WZ_D = 1000; %Диаметр цилиндра рабочей зоны
WZ_H = 400; %Высота цилиндра рабочей зоны
WZ_Z = -1150; %Координаты дна цилиндра рабочей зоны
WZ_a = WZ_D*sqrt(2)/2; %Сторона вписанного в основание цилиндра квадрата
WZ_d = 500; %Диаметр нижней окружности усечённого конуса
WZ_h = 200; %Высота добавляемой к РЗ части

%Максимальная требуемая скорость рабочего органа
vmax = 1560;
%Максимальное требуемое от рабочего органа ускорение
amax = 10400;
%Момент инерции рычага
Jr = 200000;
%Масса каретки
mk = 2;
%Размеры звеньев
F = 600; %Длина стороны треугольника основания
f = 200; %Длина стороны треугольника платформы
R_l = 630;  %Длина рычагов (по осям)
R_r = 1000;  %Длина штанг (по осям)%Константы
cos120 = cosd(120);
sin120 = sind(120);
cos240 = cosd(240);
sin240 = sind(240);
VM = f*sqrt(3)/6;
OQ = F*sqrt(3)/6;
%Генерируем массив точек для анализа
[TestPointsX, TestPointsY, TestPointsZ] = genTP();

%Задаём начальные значения
Mmax = 0;
Mmin = 1e10;

%Запускаем цикл проверки в каждой точке
for k = 1:size(TestPointsX, 2)
J = [                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          TestPointsX(k)/(R_l*(1 - ((OQ - VM + TestPointsY(k))^2 + R_l^2 - R_r^2 + TestPointsX(k)^2 + TestPointsZ(k)^2)^2/(4*R_l^2*((OQ - VM + TestPointsY(k))^2 + TestPointsZ(k)^2)))^(1/2)*((OQ - VM + TestPointsY(k))^2 + TestPointsZ(k)^2)^(1/2)),                                                                                                                                                                                                                                                                                                   ((2*OQ - 2*VM + 2*TestPointsY(k))/(2*R_l*((OQ - VM + TestPointsY(k))^2 + TestPointsZ(k)^2)^(1/2)) - ((2*OQ - 2*VM + 2*TestPointsY(k))*((OQ - VM + TestPointsY(k))^2 + R_l^2 - R_r^2 + TestPointsX(k)^2 + TestPointsZ(k)^2))/(4*R_l*((OQ - VM + TestPointsY(k))^2 + TestPointsZ(k)^2)^(3/2)))/(1 - ((OQ - VM + TestPointsY(k))^2 + R_l^2 - R_r^2 + TestPointsX(k)^2 + TestPointsZ(k)^2)^2/(4*R_l^2*((OQ - VM + TestPointsY(k))^2 + TestPointsZ(k)^2)))^(1/2) + (1/((OQ - VM + TestPointsY(k))^2 + TestPointsZ(k)^2)^(1/2) - ((2*OQ - 2*VM + 2*TestPointsY(k))*(OQ - VM + TestPointsY(k)))/(2*((OQ - VM + TestPointsY(k))^2 + TestPointsZ(k)^2)^(3/2)))/(1 - (OQ - VM + TestPointsY(k))^2/((OQ - VM + TestPointsY(k))^2 + TestPointsZ(k)^2))^(1/2),                                                                                                                                                                                                                                 (TestPointsZ(k)/(R_l*((OQ - VM + TestPointsY(k))^2 + TestPointsZ(k)^2)^(1/2)) - (TestPointsZ(k)*((OQ - VM + TestPointsY(k))^2 + R_l^2 - R_r^2 + TestPointsX(k)^2 + TestPointsZ(k)^2))/(2*R_l*((OQ - VM + TestPointsY(k))^2 + TestPointsZ(k)^2)^(3/2)))/(1 - ((OQ - VM + TestPointsY(k))^2 + R_l^2 - R_r^2 + TestPointsX(k)^2 + TestPointsZ(k)^2)^2/(4*R_l^2*((OQ - VM + TestPointsY(k))^2 + TestPointsZ(k)^2)))^(1/2) - (TestPointsZ(k)*(OQ - VM + TestPointsY(k)))/(((OQ - VM + TestPointsY(k))^2 + TestPointsZ(k)^2)^(3/2)*(1 - (OQ - VM + TestPointsY(k))^2/((OQ - VM + TestPointsY(k))^2 + TestPointsZ(k)^2))^(1/2));
    (3^(1/2)/(2*((OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)^(1/2)) - (3^(1/2)*(OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2)/(2*((OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)^(3/2)))/(1 - (OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2/((OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2))^(1/2) + ((TestPointsX(k)/2 + (3^(1/2)*TestPointsY(k))/2 + 3^(1/2)*(OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2))/(2*R_l*((OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)^(1/2)) - (3^(1/2)*(OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)*((TestPointsX(k)/2 + (3^(1/2)*TestPointsY(k))/2)^2 + (OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + R_l^2 - R_r^2 + TestPointsZ(k)^2))/(4*R_l*((OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)^(3/2)))/(1 - ((TestPointsX(k)/2 + (3^(1/2)*TestPointsY(k))/2)^2 + (OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + R_l^2 - R_r^2 + TestPointsZ(k)^2)^2/(4*R_l^2*((OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)))^(1/2), ((VM - OQ + TestPointsY(k)/2 - (3^(1/2)*TestPointsX(k))/2 + 3^(1/2)*(TestPointsX(k)/2 + (3^(1/2)*TestPointsY(k))/2))/(2*R_l*((OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)^(1/2)) + ((OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)*((TestPointsX(k)/2 + (3^(1/2)*TestPointsY(k))/2)^2 + (OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + R_l^2 - R_r^2 + TestPointsZ(k)^2))/(4*R_l*((OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)^(3/2)))/(1 - ((TestPointsX(k)/2 + (3^(1/2)*TestPointsY(k))/2)^2 + (OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + R_l^2 - R_r^2 + TestPointsZ(k)^2)^2/(4*R_l^2*((OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)))^(1/2) + ((OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2/(2*((OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)^(3/2)) - 1/(2*((OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)^(1/2)))/(1 - (OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2/((OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2))^(1/2), (TestPointsZ(k)/(R_l*((OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)^(1/2)) - (TestPointsZ(k)*((TestPointsX(k)/2 + (3^(1/2)*TestPointsY(k))/2)^2 + (OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + R_l^2 - R_r^2 + TestPointsZ(k)^2))/(2*R_l*((OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)^(3/2)))/(1 - ((TestPointsX(k)/2 + (3^(1/2)*TestPointsY(k))/2)^2 + (OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + R_l^2 - R_r^2 + TestPointsZ(k)^2)^2/(4*R_l^2*((OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)))^(1/2) - (TestPointsZ(k)*(OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2))/((1 - (OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2/((OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2))^(1/2)*((OQ - VM - TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)^(3/2));
    ((TestPointsX(k)/2 - (3^(1/2)*TestPointsY(k))/2 + 3^(1/2)*(VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2))/(2*R_l*((VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)^(1/2)) - (3^(1/2)*(VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)*((TestPointsX(k)/2 - (3^(1/2)*TestPointsY(k))/2)^2 + (VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + R_l^2 - R_r^2 + TestPointsZ(k)^2))/(4*R_l*((VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)^(3/2)))/(1 - ((TestPointsX(k)/2 - (3^(1/2)*TestPointsY(k))/2)^2 + (VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + R_l^2 - R_r^2 + TestPointsZ(k)^2)^2/(4*R_l^2*((VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)))^(1/2) - (3^(1/2)/(2*((VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)^(1/2)) - (3^(1/2)*(VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2)/(2*((VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)^(3/2)))/(1 - (VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2/((VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2))^(1/2), ((VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2 - 3^(1/2)*(TestPointsX(k)/2 - (3^(1/2)*TestPointsY(k))/2))/(2*R_l*((VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)^(1/2)) - ((VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)*((TestPointsX(k)/2 - (3^(1/2)*TestPointsY(k))/2)^2 + (VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + R_l^2 - R_r^2 + TestPointsZ(k)^2))/(4*R_l*((VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)^(3/2)))/(1 - ((TestPointsX(k)/2 - (3^(1/2)*TestPointsY(k))/2)^2 + (VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + R_l^2 - R_r^2 + TestPointsZ(k)^2)^2/(4*R_l^2*((VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)))^(1/2) + ((VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2/(2*((VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)^(3/2)) - 1/(2*((VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)^(1/2)))/(1 - (VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2/((VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2))^(1/2), (TestPointsZ(k)/(R_l*((VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)^(1/2)) - (TestPointsZ(k)*((TestPointsX(k)/2 - (3^(1/2)*TestPointsY(k))/2)^2 + (VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + R_l^2 - R_r^2 + TestPointsZ(k)^2))/(2*R_l*((VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)^(3/2)))/(1 - ((TestPointsX(k)/2 - (3^(1/2)*TestPointsY(k))/2)^2 + (VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + R_l^2 - R_r^2 + TestPointsZ(k)^2)^2/(4*R_l^2*((VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)))^(1/2) + (TestPointsZ(k)*(VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2))/((1 - (VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2/((VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2))^(1/2)*((VM - OQ + TestPointsY(k)/2 + (3^(1/2)*TestPointsX(k))/2)^2 + TestPointsZ(k)^2)^(3/2))];
%Находим максимальные скорости входных звеньев
w_max_TP = [norm(J(1, :))*vmax, norm(J(2, :))*vmax, norm(J(3, :))*vmax];
%Вычисляем приведённую массу в точке V
m_pr = mk + Jr*(w_max_TP(1)/vmax)^2 + Jr*(w_max_TP(2)/vmax)^2 + Jr*(w_max_TP(3)/vmax)^2;
%Вычисляем силу, действующую на механизм
F = m_pr*amax*(10^(-3));
J = transpose(inv(J));
M_max_TP = [norm(J(1, :))*F, norm(J(2, :))*F, norm(J(3, :))*F];
%Находим максимальный момент среди рычагов
Mmax = max([abs(M_max_TP), Mmax]);
Mmin = min([abs(M_max_TP), Mmin]);
end
%Переводим результат в H•м
maxM = Mmax*10^(-3)
%Вычисляем неравномерность максимального крутящего момента по рабочей зоне
uneven = Mmax/Mmin