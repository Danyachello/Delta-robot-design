
%Задаём символьные переменные
syms Tetta_1 Tetta_2 Tetta_3 X_V Y_V Z_V
syms OQ VM R_l R_r F f X_0 Y_0 Z_0
syms omega_1(t) v_VX(t) v_VY(t) v_VZ(t)
syms epsilon_1 a_VX a_VY a_VZ
syms v_VX v_VY v_VZ
%Записываем решение обратной кинематической задачи в аналитическом виде
NL = sqrt(R_r^2-X_V^2);
y_M = -VM + Y_V;
y_Q = -OQ;
const_1 = y_M - y_Q;
NQ = sqrt(const_1^2 + Z_V^2);
Eq1 = Tetta_1 == 2*pi - acos((R_l^2 + NQ^2 - NL^2)/(2*R_l*NQ))...
    - acos(const_1/NQ);
%Формируем уравнение связи
F_1 = 2*pi - acos((R_l^2 + NQ^2 - NL^2)/(2*R_l*NQ))...
    - acos(const_1/NQ) - Tetta_1;
% F_1 = - (R_l^2 + NQ^2 - NL^2)/(2*R_l*NQ)...
%     - const_1/NQ+cos(2*pi - Tetta_1);
    %Проделываем те же операции для 2-го плеча
    X_V_120 = X_V*cos(2*pi/3) - Y_V*sin(2*pi/3);
    Y_V_120 = X_V*sin(2*pi/3) + Y_V*cos(2*pi/3);
    Z_V_120 = Z_V;
    NL = sqrt(R_r^2-X_V_120^2);
    y_M = -VM + Y_V_120;
    y_Q = -OQ;
    const_1 = y_M - y_Q;
    NQ = sqrt(const_1^2 + Z_V_120^2);
    Eq2 = Tetta_2 == 2*pi - acos((R_l^2 + NQ^2 - NL^2)/(2*R_l*NQ))...
        - acos(const_1/NQ);
    F_2 = 2*pi - acos((R_l^2 + NQ^2 - NL^2)/(2*R_l*NQ))...
        - acos(const_1/NQ) - Tetta_2;
%     F_2 = - (R_l^2 + NQ^2 - NL^2)/(2*R_l*NQ)...
%         - const_1/NQ + cos(2*pi - Tetta_2);
        %Для 3-го плеча
        X_V_240 = X_V*cos(4*pi/3) - Y_V*sin(4*pi/3);
        Y_V_240 = X_V*sin(4*pi/3) + Y_V*cos(4*pi/3);
        Z_V_240 = Z_V;
        NL = sqrt(R_r^2-X_V_240^2);
        y_M = -VM + Y_V_240;
        y_Q = -OQ;
        const_1 = y_M - y_Q;
        NQ = sqrt(const_1^2 + Z_V_240^2);
        Eq3 = Tetta_3 == 2*pi - acos((R_l^2 + NQ^2 - NL^2)/(2*R_l*NQ))...
            - acos(const_1/NQ);
        F_3 = 2*pi - acos((R_l^2 + NQ^2 - NL^2)/(2*R_l*NQ))...
            - acos(const_1/NQ) - Tetta_3;
%                 F_3 = - (R_l^2 + NQ^2 - NL^2)/(2*R_l*NQ)...
%             - const_1/NQ + cos(2*pi - Tetta_3);

%Составляем матрицы А и В
A = [diff(F_1, X_V), diff(F_1, Y_V), diff(F_1, Z_V);
     diff(F_2, X_V), diff(F_2, Y_V), diff(F_2, Z_V);
     diff(F_3, X_V), diff(F_3, Y_V), diff(F_3, Z_V)];
B = [diff(F_1, Tetta_1), 0, 0;
     0, diff(F_2, Tetta_2), 0;
     0, 0, diff(F_3, Tetta_3)];
J = -A*(inv(B))

%Подставляем конкретные значения (для проверки)
J_ch = double(subs(J, [X_V, Y_V, Z_V, OQ, VM, R_l, R_r, F, f],...
  [0, 0, -800, 114.34, 28.87, 630, 1000, 500, 100]));

%Вычисляем скорости в тестовой точке
v = [0; 0; 1000];
w = J_ch*v;
n = 30*w/pi
%Ответ [-3.1503; -3.3611; -4.4766]