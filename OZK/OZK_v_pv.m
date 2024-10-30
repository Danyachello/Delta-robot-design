%Находит решение задачи о скоростях методом планов
function [w_1, w_2, w_3] = OZK_v_pv(X_V, Y_V, Z_V, v_VX, v_VY, v_VZ)
global R_l VM OQ cos120 sin120 cos240 sin240 %Размеры и константы
%Расчёт углов поворота рычагов в соответствующих системах координат
[vectTheta] = OZK(X_V, Y_V, Z_V);
Theta1 = vectTheta(1); Theta2 = vectTheta(2); Theta3 = vectTheta(3);
            %Вычисляем координаты точек в системах координат XOY,
            %X120Y120Z120 и X240Y240Z240
            Y_Q1 = -OQ;
            Z_Q1 = 0;
            X_L1 = 0;
            Y_L1 = -(OQ + R_l*sind(Theta1-90));
            Z_L1 =  R_l*cosd(Theta1-90);
            X_M1 = X_V;
            Y_M1 = Y_V - VM;
            Z_M1 = Z_V;
                X_V_120 = X_V*cos120 - Y_V*sin120;
                Y_V_120 = X_V*sin120 + Y_V*cos120;
                Z_V_120 = Z_V;
                Y_Q2_120 = -OQ;
                Z_Q2_120 = 0;
                X_L2_120 = 0;
                Y_L2_120 = -(OQ + R_l*sind(Theta2-90));
                Z_L2_120 =  R_l*cosd(Theta2-90);
                X_M2_120 = X_V_120;
                Y_M2_120 = Y_V_120 - VM;
                Z_M2_120 = Z_V_120;
                    X_V_240 = X_V*cos240 - Y_V*sin240;
                    Y_V_240 = X_V*sin240 + Y_V*cos240;
                    Z_V_240 = Z_V;
                    Y_Q3_240 = -OQ;
                    Z_Q3_240 = 0;
                    X_L3_240 = 0;
                    Y_L3_240 = -(OQ + R_l*sind(Theta3-90));
                    Z_L3_240 =  R_l*cosd(Theta3-90);
                    X_M3_240 = X_V_240;
                    Y_M3_240 = Y_V_240 - VM;
                    Z_M3_240 = Z_V_240;
                        v_VX_120 = v_VX*cos120 - v_VY*sin120;
                        v_VY_120 = v_VX*sin120 + v_VY*cos120;
                        v_VZ_120 = v_VZ;
                            v_VX_240 = v_VX*cos240 - v_VY*sin240;
                            v_VY_240 = v_VX*sin240 + v_VY*cos240;
                            v_VZ_240 = v_VZ;



x = (v_VX*(X_L1 - X_M1) + v_VY*(Y_L1 - Y_M1) + v_VZ*(Z_L1 - Z_M1))/(Y_M1 - Y_L1 + ((Y_L1 - Y_Q1)*(Z_L1 - Z_M1))/(Z_L1 - Z_Q1));
y = (Y_Q1-Y_L1)/(Z_Q1-Z_L1)*x;
v_L1 = sqrt(x^2+y^2);
if y>0
    znak = -1;
else
    znak = 1;
end
w_1 = znak*v_L1/R_l;
    x = (v_VX_120*(X_L2_120 - X_M2_120) + v_VY_120*(Y_L2_120 - Y_M2_120) + v_VZ_120*(Z_L2_120 - Z_M2_120))/(Y_M2_120 - Y_L2_120 + ((Y_L2_120 - Y_Q2_120)*(Z_L2_120 - Z_M2_120))/(Z_L2_120 - Z_Q2_120));
    y = (Y_Q2_120-Y_L2_120)/(Z_Q2_120-Z_L2_120)*x;
    v_L2_120 = sqrt(x^2+y^2);
    if y>0
        znak = -1;
    else
        znak = 1;
    end
    w_2 = znak*v_L2_120/R_l;
        x = (v_VX_240*(X_L3_240 - X_M3_240) + v_VY_240*(Y_L3_240 - Y_M3_240) + v_VZ_240*(Z_L3_240 - Z_M3_240))/(Y_M3_240 - Y_L3_240 + ((Y_L3_240 - Y_Q3_240)*(Z_L3_240 - Z_M3_240))/(Z_L3_240 - Z_Q3_240));
        y = (Y_Q3_240-Y_L3_240)/(Z_Q3_240-Z_L3_240)*x;
        v_L3_240 = sqrt(x^2+y^2);
        if y>0
            znak = -1;
        else
            znak = 1;
        end
        w_3 = znak*v_L3_240/R_l;
end