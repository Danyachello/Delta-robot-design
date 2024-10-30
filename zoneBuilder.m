%Функция построения реальной и желаемой рабочих областей
%дельта-робота
function [] = zoneBuilder()
global R_l R_r VM OQ cos120 sin120 cos240 sin240 %Размеры и константы
global minTheta QG phiMax varthetaMax %Ограничения
global WZ_D WZ_H WZ_Z WZ_a WZ_d WZ_h %Рабочая зона
global stepTheta %Шаг изменения углов
global minTheta_real maxTheta_real %Максимальный и минимальный углы рычагов
%при движении в реальной рабочей области
global minTheta_wish maxTheta_wish %Максимальный и минимальный углы рычагов
%при движении в желаемой рабочей области
global addSeg %Переменная, показывающая необходимость дополнительного
%сегмента в рабочей зоне
global status %Флаг свидетельствующий об удовлетворении условия нахождения
%желаемой рабочей зоны в реальной рабочей зоне
%поворота рычагов, которые могут быть с учёном ограничений
global Dg Hg hg %Габаритные размеры почившегося робота
maxTheta = 270; %Максимальный угол для перебора

%Задаём переменные максимальных и минимальных углов поворота рычагов,
%которые могут быть с учёном ограничений
minTheta_real = maxTheta;
maxTheta_real = minTheta;

minTheta_wish = maxTheta;
maxTheta_wish = minTheta;


%Подготовительные вычисления
gammaMin = asind(QG/R_l);
gammaMax = 180;
%Создаём вектор с одной нулевой точкой (потом её удалим)
V_RZ = [0, 0, 0];

%Создаём цикл перебора всех возможных комбинаций углов поворота рычагов
for Theta1 = minTheta:stepTheta:maxTheta
    for Theta2 = minTheta:stepTheta:maxTheta
        for Theta3 = minTheta:stepTheta:maxTheta
            %Вызываем функцию, решающую ПЗК
            [~,~,~,V] = PZK(Theta1, Theta2, Theta3);
            X_V = V(1); Y_V = V(2); Z_V = V(3);
            %Вычисляем координаты точек в системах координат XOY
            %X120Y120Z120 и X240Y240Z240
            X_Q1 = 0;
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
                X_Q2_120 = 0;
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
                    X_Q3_240 = 0;
                    Y_Q3_240 = -OQ;
                    Z_Q3_240 = 0;
                    X_L3_240 = 0;
                    Y_L3_240 = -(OQ + R_l*sind(Theta3-90));
                    Z_L3_240 =  R_l*cosd(Theta3-90);
                    X_M3_240 = X_V_240;
                    Y_M3_240 = Y_V_240 - VM;
                    Z_M3_240 = Z_V_240;
            %Собираем это в вектора
            Q1 = [X_Q1, Y_Q1, Z_Q1];
            L1 = [X_L1, Y_L1, Z_L1];
            M1 = [X_M1, Y_M1, Z_M1];
            V_120 = [X_V_120, Y_V_120, Z_V_120];
            Q2_120 = [X_Q2_120, Y_Q2_120, Z_Q2_120];
            L2_120 = [X_L2_120, Y_L2_120, Z_L2_120];
            M2_120 = [X_M2_120, Y_M2_120, Z_M2_120];
            V_240 = [X_V_240, Y_V_240, Z_V_240];
            Q3_240 = [X_Q3_240, Y_Q3_240, Z_Q3_240];
            L3_240 = [X_L3_240, Y_L3_240, Z_L3_240];
            M3_240 = [X_M3_240, Y_M3_240, Z_M3_240];
            %Вычисляем координаты точек N1 N2 N3
            N1 = [L1(1); M1(2); V(3)];
            N2_120 = [L2_120(1); M2_120(2); V_120(3)];
            N3_240 = [L3_240(1); M3_240(2); V_240(3)];
            %------Проверяем ограничения------%
            % 1 Проверяем отсутствие пересечений
            gamma1 = 180-acosd(-L1(3)/R_l)-asind((N1(2)-L1(2))/R_r);
            gamma2 = 180-acosd(-L2_120(3)/R_l)-asind((N2_120(2)-L2_120(2))/R_r);
            gamma3 = 180-acosd(-L3_240(3)/R_l)-asind((N3_240(2)-L3_240(2))/R_r);
            if (gamma1 < gammaMax) && (gamma1 > gammaMin) && (gamma2 < gammaMax) && (gamma2 > gammaMin) && (gamma3 < gammaMax) && (gamma3 > gammaMin)
                % 2 Проверяем, могут ли трёхподвижные шарниры позволить механизму принять
                % такое положение
                varphi1 = asind(abs(V(1))/R_r);
                varphi2 = asind(abs(V_120(1))/R_r);
                varphi3 = asind(abs(V_240(1))/R_r);
                if (varphi1<phiMax) && (varphi2<phiMax) && (varphi3<phiMax)
                    % 3 Проверяем, не превышают ли углы давления заданное
                    % максимальное значение
                    p1 = [M1(1)-L1(1); M1(2)-L1(2); M1(3)-L1(3)];
                    q1 = [0; L1(3)-Q1(3); Q1(2)-L1(2)];
                    vartheta1 = acosd(abs(p1(1)*q1(1)+p1(2)*q1(2)+p1(3)*q1(3))/(sqrt(p1(1)^2+p1(2)^2+p1(3)^2)*sqrt(q1(1)^2+q1(2)^2+q1(3)^2)));
                    p2 = [M2_120(1)-L2_120(1); M2_120(2)-L2_120(2); M2_120(3)-L2_120(3)];
                    q2 = [0; L2_120(3)-Q2_120(3); Q2_120(2)-L2_120(2)];
                    vartheta2 = acosd(abs(p2(1)*q2(1)+p2(2)*q2(2)+p2(3)*q2(3))/(sqrt(p2(1)^2+p2(2)^2+p2(3)^2)*sqrt(q2(1)^2+q2(2)^2+q2(3)^2)));
                    p3 = [M3_240(1)-L3_240(1); M3_240(2)-L3_240(2); M3_240(3)-L3_240(3)];
                    q3 = [0; L3_240(3)-Q3_240(3); Q3_240(2)-L3_240(2)];
                    vartheta3 = acosd(abs(p3(1)*q3(1)+p3(2)*q3(2)+p3(3)*q3(3))/(sqrt(p3(1)^2+p3(2)^2+p3(3)^2)*sqrt(q3(1)^2+q3(2)^2+q3(3)^2)));
                    if (vartheta1<varthetaMax) && (vartheta2<varthetaMax) && (vartheta3<varthetaMax)
                        %Если все условия выполнены - добавляем точку в
                        %вектор облака точек реальной РЗ
                        V_RZ = cat(1, V_RZ, V);
                        %Находим максимальный и минимальный углы поворота рычагов,
                        %которые могут быть с учёном ограничений
                        minTheta_real = min([minTheta_real, Theta1, Theta2, Theta3]);
                        maxTheta_real = max([maxTheta_real, Theta1, Theta2, Theta3]);
                    end
                end
            end
        end
    end
end

%Удаляем первую нулевую точку
V_RZ(1,:) = [];
%Находим какие точки в нашем облаке являются граничными и создаём
%группы по 3 точки, которые будем соединять треугольниками (для
%образования поверхности)
[K] = boundary(V_RZ);

%Находим точки для построения цилиндра
%Крышка
t = linspace(0, 2*pi, 50);
Xo1 = (WZ_D/2) * cos(t);
Yo1 = (WZ_D/2) * sin(t);
Zo1 = ones(1, 50)*(WZ_Z+WZ_H);
Zo2 = ones(1, 50)*WZ_Z;
%Цилиндрическая поверхность
[Xc,Yc,Zc] = cylinder((WZ_D/2), 50);
Zc = Zc*WZ_H+WZ_Z;

%Находим точки для построения параллелепипеда
Xp = [0 WZ_a WZ_a 0 0 WZ_a WZ_a 0];
Yp = [0 0 WZ_a WZ_a 0 0 WZ_a WZ_a] ;
Zp = Xp'*Xp*Yp'*Yp/WZ_a^2/WZ_a^2*WZ_H/2 + WZ_Z;
Xp = Xp - WZ_a/2;
Yp = Yp - WZ_a/2;

%Отрисовка рабочих зон

%Отображаем реальную РЗ
p1 = trisurf(K,V_RZ(:,1),V_RZ(:,2),V_RZ(:,3),'Facecolor','red', 'FaceAlpha',0.2);
hold('on');
%Отображаем крышку и боковую поверхность цилиндра
p2 = patch(Xo1, Yo1, Zo1, 'c', 'FaceAlpha',0.5);
p3 = surf(Xc,Yc,Zc, 'FaceColor', 'c',  'FaceAlpha',0.5);
%Отображаем параллелепипед
p4 = surf(Xp, Yp, Zp, 'FaceColor', 'g', 'FaceAlpha', 0.2);
p4.Visible = 0; %Скрываем
%Отображаем дополнительные сегменты, если нужно
if addSeg == 0
    %Дно цилиндра
    p5 = patch(Xo1, Yo1, Zo2, 'c', 'FaceAlpha',0.5);
end
if addSeg == 1
    %Дно усечённого конуса
    t = linspace(0, 2*pi, 50);
    Xo2 = (WZ_d/2) * cos(t);
    Yo2 = (WZ_d/2) * sin(t);
    p6 = patch(Xo2, Yo2, Zo2-WZ_h, 'c', 'FaceAlpha',0.5);
    %Коническая поверхность
    [Xco,Yco,Zco] = cylinder([(WZ_d/2), (WZ_D/2)], 50);
    Zco = Zco*WZ_h+WZ_Z-WZ_h;
    p7 = surf(Xco,Yco,Zco, 'FaceColor', 'c',  'FaceAlpha',0.5);
end
if addSeg == 2
    %Сферическая поверхность
    thetta = linspace(pi-asin((WZ_D/2)/(WZ_d/2)), pi, 50);
    phi = linspace(0, 2*pi, 50);
    [THETTA, PHI] = meshgrid(thetta, phi);
    Xsph = (WZ_d/2)*sin(THETTA).*cos(PHI);
    Ysph = (WZ_d/2)*sin(THETTA).*sin(PHI);
    Zsph = (WZ_Z+(WZ_d/2-WZ_h))+(WZ_d/2)*cos(THETTA);
    p8 = surf(Xsph,Ysph,Zsph, 'FaceColor', 'c',  'FaceAlpha',0.5);
end
%Настраиваем отображение
hold('off');
axis('equal');
xlabel('x, мм');
ylabel('y, мм');
zlabel('z, мм');
title('Рабочие зоны механизма');
%legend('Реальная РЗ', 'Требуемая РЗ');
view([-30 , 20]);
grid('on');

%Проверяем, подходит ли данная конфигурация робота под желаемую рабочую
%зону
%Формируем массив точек для проверки
%Точки крышки цилиндра
[t, r] = meshgrid(linspace(-pi/2, -pi/2-pi/3, 10), linspace(0, WZ_D/2,10));
[Xtc, Ytc] = pol2cart(t,r);
Ztc = ones(10, 10)*(WZ_Z+WZ_H);
%Точки боковой поверхности цилиндра
[t,r] = meshgrid(linspace(-pi/2, -pi/2-pi/3, 10), WZ_D/2);
[Xmv, Ymv] = pol2cart(t,r);
mv = linspace(WZ_Z, WZ_Z+WZ_H, 10);
for k = 2:9
    Xss((k-1)*10-9:(k-1)*10) = Xmv;
    Yss((k-1)*10-9:(k-1)*10) = Ymv;
    Zss((k-1)*10-9:(k-1)*10) = ones(1, 10)*mv(k);
end
if addSeg == 0
    %Точки дна цилиндра
    [t, r] = meshgrid(linspace(-pi/2, -pi/2-pi/3, 10), linspace(0, WZ_D/2,10));
    [Xbc, Ybc] = pol2cart(t,r);
    Zbc = ones(10, 10)*(WZ_Z);
    %Формируем итоговый вектор точек
    PointsToBeChecked_X = [reshape(Xtc',1,[]), reshape(Xss',1,[]), reshape(Xbc',1,[])];
    PointsToBeChecked_Y = [reshape(Ytc',1,[]), reshape(Yss',1,[]), reshape(Ybc',1,[])];
    PointsToBeChecked_Z = [reshape(Ztc',1,[]), reshape(Zss',1,[]), reshape(Zbc',1,[])];
end
if addSeg == 1
    %Точки дна конуса
    [t, r] = meshgrid(linspace(-pi/2, -pi/2-pi/3, 10), linspace(0, WZ_d/2,10));
    [Xbco, Ybco] = pol2cart(t,r);
    Zbco = ones(10, 10)*(WZ_Z-WZ_h);
    %Точки боковой поверхности конуса
    mv = linspace(WZ_Z-WZ_h, WZ_Z, 10);
    mvD = linspace(WZ_d, WZ_D, 10);
    for k = 2:10
        [t,r] = meshgrid(linspace(-pi/2, -pi/2-pi/3, 10), mvD(k)/2);
        [Xmv, Ymv] = pol2cart(t,r);
        Xssco((k-1)*10-9:(k-1)*10) = Xmv;
        Yssco((k-1)*10-9:(k-1)*10) = Ymv;
        Zssco((k-1)*10-9:(k-1)*10) = ones(1, 10)*mv(k);
    end
    %Формируем итоговый вектор точек
    PointsToBeChecked_X = [reshape(Xtc',1,[]), reshape(Xss',1,[]), reshape(Xbco',1,[]), reshape(Xssco',1,[])];
    PointsToBeChecked_Y = [reshape(Ytc',1,[]), reshape(Yss',1,[]), reshape(Ybco',1,[]), reshape(Yssco',1,[])];
    PointsToBeChecked_Z = [reshape(Ztc',1,[]), reshape(Zss',1,[]), reshape(Zbco',1,[]), reshape(Zssco',1,[])];
end
if addSeg == 2
    %Точки сферы
    thetta = linspace(pi-asin((WZ_D/2)/(WZ_d/2)), pi, 10);
    phi = linspace(-pi/2, -pi/2-pi/3, 10);
    [THETTA, PHI] = meshgrid(thetta, phi);
    Xsph = (WZ_d/2)*sin(THETTA).*cos(PHI);
    Ysph = (WZ_d/2)*sin(THETTA).*sin(PHI);
    Zsph = (WZ_Z+(WZ_d/2-WZ_h))+(WZ_d/2)*cos(THETTA);
    %Формируем итоговый вектор точек
    PointsToBeChecked_X = [reshape(Xtc',1,[]), reshape(Xss',1,[]), reshape(Xsph',1,[])];
    PointsToBeChecked_Y = [reshape(Ytc',1,[]), reshape(Yss',1,[]), reshape(Ysph',1,[])];
    PointsToBeChecked_Z = [reshape(Ztc',1,[]), reshape(Zss',1,[]), reshape(Zsph',1,[])];
end

% %Показываем точки для анализа
% hold('on');
% plot3(PointsToBeChecked_X,PointsToBeChecked_Y, PointsToBeChecked_Z, '.r');
% hold('off');
    
%Проверяем каждую точку на условие нахождения её в реальной рабочей зоне
%Задаём начальный вектор статуса (ни одна точка не подходит)
stt(1:size(PointsToBeChecked_X, 2)) = 0;
for k=1:size(PointsToBeChecked_X, 2)
    %Вызываем функцию, решающую ПЗК
    vectTheta = [0, 0, 0];
    [vectTheta(1),vectTheta(2), vectTheta(3)] = OZK(PointsToBeChecked_X(k), PointsToBeChecked_Y(k), PointsToBeChecked_Z(k));
    Theta1 = vectTheta(1); Theta2 = vectTheta(2); Theta3 = vectTheta(3);
    X_V = PointsToBeChecked_X(k); Y_V = PointsToBeChecked_Y(k); Z_V = PointsToBeChecked_Z(k);
    if imag(vectTheta) == [0, 0, 0]
        %Вычисляем координаты точек в системах координат XOY
        %X120Y120Z120 и X240Y240Z240
        X_Q1 = 0;
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
        X_Q2_120 = 0;
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
        X_Q3_240 = 0;
        Y_Q3_240 = -OQ;
        Z_Q3_240 = 0;
        X_L3_240 = 0;
        Y_L3_240 = -(OQ + R_l*sind(Theta3-90));
        Z_L3_240 =  R_l*cosd(Theta3-90);
        X_M3_240 = X_V_240;
        Y_M3_240 = Y_V_240 - VM;
        Z_M3_240 = Z_V_240;
        %Собираем это в вектора
        Q1 = [X_Q1, Y_Q1, Z_Q1];
        L1 = [X_L1, Y_L1, Z_L1];
        M1 = [X_M1, Y_M1, Z_M1];
        V_120 = [X_V_120, Y_V_120, Z_V_120];
        Q2_120 = [X_Q2_120, Y_Q2_120, Z_Q2_120];
        L2_120 = [X_L2_120, Y_L2_120, Z_L2_120];
        M2_120 = [X_M2_120, Y_M2_120, Z_M2_120];
        V_240 = [X_V_240, Y_V_240, Z_V_240];
        Q3_240 = [X_Q3_240, Y_Q3_240, Z_Q3_240];
        L3_240 = [X_L3_240, Y_L3_240, Z_L3_240];
        M3_240 = [X_M3_240, Y_M3_240, Z_M3_240];
        %Вычисляем координаты точек N1 N2 N3
        N1 = [L1(1); M1(2); V(3)];
        N2_120 = [L2_120(1); M2_120(2); V_120(3)];
        N3_240 = [L3_240(1); M3_240(2); V_240(3)];
        %------Проверяем ограничения------%
        % 1 Проверяем отсутствие пересечений
        gamma1 = 180-acosd(-L1(3)/R_l)-asind((N1(2)-L1(2))/R_r);
        gamma2 = 180-acosd(-L2_120(3)/R_l)-asind((N2_120(2)-L2_120(2))/R_r);
        gamma3 = 180-acosd(-L3_240(3)/R_l)-asind((N3_240(2)-L3_240(2))/R_r);
        if (gamma1 < gammaMax) && (gamma1 > gammaMin) && (gamma2 < gammaMax) && (gamma2 > gammaMin) && (gamma3 < gammaMax) && (gamma3 > gammaMin)
            % 2 Проверяем, могут ли трёхподвижные шарниры позволить механизму принять
            % такое положение
            varphi1 = asind(abs(V(1))/R_r);
            varphi2 = asind(abs(V_120(1))/R_r);
            varphi3 = asind(abs(V_240(1))/R_r);
            if (varphi1<phiMax) && (varphi2<phiMax) && (varphi3<phiMax)
                % 3 Проверяем, не превышают ли углы давления заданное
                % максимальное значение
                p1 = [M1(1)-L1(1); M1(2)-L1(2); M1(3)-L1(3)];
                q1 = [0; L1(3)-Q1(3); Q1(2)-L1(2)];
                vartheta1 = acosd(abs(p1(1)*q1(1)+p1(2)*q1(2)+p1(3)*q1(3))/(sqrt(p1(1)^2+p1(2)^2+p1(3)^2)*sqrt(q1(1)^2+q1(2)^2+q1(3)^2)));
                p2 = [M2_120(1)-L2_120(1); M2_120(2)-L2_120(2); M2_120(3)-L2_120(3)];
                q2 = [0; L2_120(3)-Q2_120(3); Q2_120(2)-L2_120(2)];
                vartheta2 = acosd(abs(p2(1)*q2(1)+p2(2)*q2(2)+p2(3)*q2(3))/(sqrt(p2(1)^2+p2(2)^2+p2(3)^2)*sqrt(q2(1)^2+q2(2)^2+q2(3)^2)));
                p3 = [M3_240(1)-L3_240(1); M3_240(2)-L3_240(2); M3_240(3)-L3_240(3)];
                q3 = [0; L3_240(3)-Q3_240(3); Q3_240(2)-L3_240(2)];
                vartheta3 = acosd(abs(p3(1)*q3(1)+p3(2)*q3(2)+p3(3)*q3(3))/(sqrt(p3(1)^2+p3(2)^2+p3(3)^2)*sqrt(q3(1)^2+q3(2)^2+q3(3)^2)));
                if (vartheta1<varthetaMax) && (vartheta2<varthetaMax) && (vartheta3<varthetaMax)
                    %Если все условия выполнены - добавляем точку в
                    %вектор облака точек реальной РЗ
                    stt(k) = 1;
                    %Находим максимальный и минимальный углы поворота рычагов,
                    %которые могут быть с учёном ограничений
                    minTheta_wish = min([minTheta_wish, Theta1, Theta2, Theta3]);
                    maxTheta_wish = max([maxTheta_wish, Theta1, Theta2, Theta3]);
                end
            end
        end
    end
end
if min(stt) == 1
    status = 1;
    Dg = (OQ+R_l)*2;
    hg = R_l*sind(180-minTheta_wish);
    if hg < 0
        hg = 0;
    end
    Hg = hg - WZ_Z;
else
    status = 0;
end
end