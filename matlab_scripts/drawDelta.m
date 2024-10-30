%Строит схаматично дельта-робота на графике
function [] = drawDelta(dim, V, Theta)
f = dim(1);    %Длина стороны треугольника платформы
F = dim(2);    %Длина стороны треугольника основания
R_l = dim(3);  %Длина рычагов (по осям)
R_r = dim(4);  %Длина штанг (по осям)
X_V = V(1); Y_V = V(2); Z_V = V(3); %Координаты точки V
Theta1 = Theta(1); Theta2 = Theta(2); Theta3 = Theta(3); %Углы поворота рычагов

%Находим все точки для построения механизма
%Подготовительные действия
global cos120 sin120 cos240 sin240
VM = f*sqrt(3)/6;
OQ = F*sqrt(3)/6;
%Находим крайние точки рычагов и плеч
X_Q1 = 0;
Y_Q1 = -OQ;
Z_Q1 = 0;
X_L1 = 0;
Y_L1 = -(OQ + R_l*sind(Theta1-90));
Z_L1 =  R_l*cosd(Theta1-90);
X_M1 = X_V;
Y_M1 = Y_V - VM;
Z_M1 = Z_V;
    X_Q2 = -Y_Q1*sin120;
    Y_Q2 = Y_Q1*cos120;
    Z_Q2 = 0;
    X_L2 = (OQ + R_l*sind(Theta2-90))*sin120;
    Y_L2 = -(OQ + R_l*sind(Theta2-90))*cos120;
    Z_L2 =  R_l*cosd(Theta2-90);
    X_M2 = X_V - (-VM)*sin120;
    Y_M2 = Y_V + (-VM)*cos120;
    Z_M2 = Z_V;
        X_Q3 = -Y_Q1*sin240;
        Y_Q3 = Y_Q1*cos240;
        Z_Q3 = 0;
        X_L3 = (OQ + R_l*sind(Theta3-90))*sin240;
        Y_L3 = -(OQ + R_l*sind(Theta3-90))*cos240;
        Z_L3 =  R_l*cosd(Theta3-90);
        X_M3 = X_V - (-VM)*sin240;
        Y_M3 = Y_V + (-VM)*cos240;
        Z_M3 = Z_V;
%Находим точки треугольника основания
A = [-F/2, -OQ, 0];
B = [F/2, A(2), 0];
C = [0, F/sqrt(3), 0];
%Находим точки треугольника каретки
An = [X_V-(f/2), Y_V-VM, Z_V];
Bn = [X_V+(f/2), Y_V-VM, Z_V];
Cn = [X_V, Y_V+f/sqrt(3), Z_V];
%Формирование выходных векторов для построения
BaseTri = [A(1) B(1) C(1); A(2) B(2) C(2); A(3) B(3) C(3)];
CarTri = [An(1) Bn(1) Cn(1); An(2) Bn(2) Cn(2); An(3) Bn(3) Cn(3)];
Pl1 = [X_Q1 X_L1 X_M1; Y_Q1 Y_L1 Y_M1; Z_Q1 Z_L1 Z_M1];
Pl2 = [X_Q2 X_L2 X_M2; Y_Q2 Y_L2 Y_M2; Z_Q2 Z_L2 Z_M2];
Pl3 = [X_Q3 X_L3 X_M3; Y_Q3 Y_L3 Y_M3; Z_Q3 Z_L3 Z_M3];

%Построение робота
cla;
patch(BaseTri(1, :), BaseTri(2, :), BaseTri(3, :), 'c', 'FaceAlpha',0.7);
hold('on');
patch(CarTri(1, :), CarTri(2, :), CarTri(3, :), 'm', 'FaceAlpha',0.7);
plot3(Pl1(1, :), Pl1(2, :), Pl1(3, :), 'LineWidth', 1.5, 'Color', 'b');
plot3(Pl2(1, :), Pl2(2, :), Pl2(3, :), 'LineWidth', 1.5, 'Color', 'b');
plot3(Pl3(1, :), Pl3(2, :), Pl3(3, :), 'LineWidth', 1.5, 'Color', 'b');
xlabel('x, мм');
ylabel('y, мм');
zlabel('z, мм');
title('Вид робота');
view([ -30 , 20 ]);
axis('equal');
grid('on');
hold('off');
end