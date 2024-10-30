%Генерирует массив точек для анализа
function [TestPointsX, TestPointsY, TestPointsZ] = genTP()
%Добавляем глобальные переменные входных и выходных данных
global WZ_D WZ_H WZ_Z WZ_d WZ_h %Рабочая зона
global addSeg %Переменная, показывающая необходимость дополнительного
%сегмента в рабочей зоне
global dotDensity %Плотность точек

%Для тестирования
dotDensity = 30;
%Задаём характеристики желаемой рабочей зоны
addSeg = 1; %Добавляемый к рабочей зоне сегмент. 0 - ничего; 1 - усечённый
%конус; %2 - часть сферы
WZ_D = 1200; %Диаметр цилиндра рабочей зоны
WZ_H = 400; %Высота цилиндра рабочей зоны
WZ_Z = -1150; %Координаты дна цилиндра рабочей зоны
WZ_a = WZ_D*sqrt(2)/2; %Сторона вписанного в основание цилиндра квадрата
WZ_d = 500; %Диаметр нижней окружности усечённого конуса
WZ_h = 200; %Высота добавляемой к РЗ части

%Задаём массив точек, заведомо покрывающий интересующий нас участок рабочей
%зоны
x = linspace(-WZ_D/2, 0, dotDensity);
y = linspace(-WZ_D/2, 0, dotDensity);
z = linspace(WZ_Z-WZ_h, WZ_Z+WZ_H, dotDensity);
[X, Y, Z] = meshgrid(x, y, z);
X = reshape(X, [1, size(X, 1)^3, 1]);
Y = reshape(Y, [1, size(Y, 1)^3, 1]);
Z = reshape(Z, [1, size(Z, 1)^3, 1]);

%Находим индексы точек, попадающих в цилиндрическую область
inCil = X.^2+Y.^2 <= (WZ_D/2)^2;
%И попадающих в сектр 60°
inSec = Y < tand(30)*X;
%Находим инексы точек попадающих
%1 Только в цилиндрическую РЗ
if addSeg == 0
    inDopSeg = Z > WZ_Z;
end
%В цилиндрическую РЗ и усёчённый конус
if addSeg == 1
    z0 = WZ_Z - WZ_h*(WZ_D/2)/((WZ_D-WZ_d)/2);
    a2_c2 = ((WZ_D/2)/(WZ_Z-z0))^2;
    inDopSeg = X.^2 + Y.^2 <= a2_c2*(Z-z0).^2;
end
%В цилиндрическую РЗ и в часть сферы
if addSeg == 2
    inDopSeg = (X.^2 + Y.^2 + (Z-(WZ_Z-WZ_h+WZ_d/2)).^2 <= (WZ_d/2)^2)|(Z > WZ_Z);
end
%Формируем индексы точек, удовлетворяющих всем условиям
inWZ = inCil & inDopSeg & inSec;
TestPointsX = X(inWZ);
TestPointsY = Y(inWZ);
TestPointsZ = Z(inWZ);

%Тестовая отрисовка точек
plot3(X(inWZ), Y(inWZ), Z(inWZ), '.r');
xlabel('x, мм');
ylabel('y, мм');
zlabel('z, мм');
grid on;
end