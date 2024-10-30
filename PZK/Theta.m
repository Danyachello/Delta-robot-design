%Находит углы поворота рычагов в соответствующих системах координат
function [Theta] = Theta(X_V, Y_V, Z_V)
global R_l R_r  VM OQ
y_M = -VM + Y_V;
y_Q = -OQ;
%Первый метод вычисления
sigma = R_l^2 - R_r^2 + X_V^2 + y_M^2 - y_Q^2 + Z_V^2;
a = (2*y_M - 2*y_Q)^2/(4*Z_V^2) + 1;
b = - 2*y_Q - ((2*y_M - 2*y_Q)*sigma)/(2*Z_V^2);
c = sigma^2/(4*Z_V^2) - R_l^2 + y_Q^2;
y_L = (-b-sqrt(b^2-4*a*c))/(2*a);
z_L = ( - 2*y_L*y_M  + 2*y_L*y_Q + sigma)/(2*Z_V);
Theta = 180 + atand((-z_L)./(y_Q-y_L));
%Второй метод вычисления
% NL = sqrt(R_r^2-X_V^2);
% const_1 = y_M - y_Q;
% NQ = sqrt(const_1^2 + Z_V^2);
% Theta = 360 - acosd((R_l^2 + NQ^2 - NL^2)/(2*R_l*NQ)) - acosd(const_1/NQ);
end