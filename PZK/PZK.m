%Находит решение прямой задачи кинематически для дельта-робота
%Возвращает координаты точек L1, L2, L3, необходимые для построения робота и
%координаты точки V
function [L1, L2, L3, V] = PZK(Tetta1, Tetta2, Tetta3)
global R_l R_r VM OQ
%Расчёт координат концов рычагов
x_L1 = 0;
y_L1 = -OQ-R_l*cosd(Tetta1-180);
z_L1 = -R_l*sind(Tetta1-180);
z_L2 = -R_l*sind(Tetta2-180);
z_L3 = -R_l*sind(Tetta3-180);
y_L2S = -OQ-R_l*cosd(Tetta2-180);
y_L3S = -OQ-R_l*cosd(Tetta3-180);
x_L2 = y_L2S*sind(120);
y_L2 = y_L2S*cosd(120);
x_L3 = y_L3S*sind(240);
y_L3 = y_L3S*cosd(240);
%Расчёт координат центров сфер (сдвинутых концов рычагов)
x_P1 = x_L1;
y_P1 = y_L1+VM;
z_P1 = z_L1;
x_P2 = x_L2+VM*cosd(30);
y_P2 = y_L2-VM*sind(30);
z_P2 = z_L2;
x_P3 = x_L3-VM*cosd(30);
y_P3 = y_L3-VM*sind(30);
z_P3 = z_L3;
%Расчёт точки пересечения сфер
w1 = R_r^2-x_P1^2-y_P1^2-z_P1^2;
w2 = R_r^2-x_P2^2-y_P2^2-z_P2^2;
w3 = R_r^2-x_P3^2-y_P3^2-z_P3^2;
a1 = x_P2-x_P1;
a2 = x_P3-x_P1;
b1 = y_P2-y_P1;
b2 = y_P3-y_P1;
c1 = z_P2-z_P1;
c2 = z_P3-z_P1;
d1 = (w1-w2)/2;
d2 = (w1-w3)/2;
e1 = ((b1*c2 - b2*c1)/(a1*b2 - a2*b1));
f1 = -(b1*d2 - b2*d1)/(a1*b2 - a2*b1);
e2 = -(a1*c2 - a2*c1)/(a1*b2 - a2*b1);
f2 = (a1*d2 - a2*d1)/(a1*b2 - a2*b1);
a_KU = (e1^2 + e2^2 + 1);
b_KU = (2*e1*(f1 - x_P1) - 2*z_P1 + 2*e2*(f2 - y_P1));
c_KU = z_P1^2 + (f1 - x_P1)^2 + (f2 - y_P1)^2 - R_r^2;
z_V = ((-b_KU-sqrt(b_KU^2-4*a_KU*c_KU)))/(2*a_KU);
x_V = e1*z_V+f1;
y_V = e2*z_V+f2;
L1 = [x_L1, y_L1, z_L1];
L2 = [x_L2, y_L2, z_L2];
L3 = [x_L3, y_L3, z_L3];
V = [x_V, y_V, z_V];
end