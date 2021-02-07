% SK                    %%%%
% 27-October-2017       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Simulator_SK
%%INPUT DATA
% Clear all
close all;
clear all;
clc;

% Degrees to Radians unit conversion
RadToDeg = (pi/180);

% Input Data
a = 40;        % L2, [ft]
b = 100;       % L3, [ft]
c = 80;        % L4, [ft]
d = 100;       % L1, [ft]
Theta2 = 20 * RadToDeg;   % Radians
Omega2 = 25;   % [Rad/sec]
Alpha2 = 15;   % [Rad/sec^2]


Calc4BarPosition(a, b, c, d, Theta2);

    function [discr3, discr4, rc] = Calc4BarPosition(a, b, c, d, Theta2)
        
    rc = 1;
    FlagOpenCrossed = 1;
    
    %Calc for Flag OpenCrossed
    if (FlagOpenCrossed == 1)
       PlusMinus = -1;
    elseif (FlagOpenCrossed == 2)
        PlusMinus = +1;
    else
        disp('Input Error for Flag OpenCrossed')
        return;
    end
   
    
    % Calculate Defined Parameters
    % Calculate Link Ratios, from Eq. 4.8a
    K1 = d/a;
    K2 = d/c;
    K3 = ((a^2) - (b^2) + (c^2) + (d^2))/(2*a*c);

    % Find Intermediate Parameters: A, B, C from eq. 4.10a
    A = cos(Theta2) - K1 - K2*cos(Theta2) + K3;
    B = -2*sin(Theta2);
    C = K1 - (K2 + 1)*cos(Theta2) + K3;

    % Calculate Link Ratios, from Eq. 4.11b
    K4 = (d/b);
    K5 = ((c^2)-(d^2)-(a^2)-(b^2))/(2*a*b);

    % Find Intermediate Parameters: D, E, F from eq. 4.12
    D = cos(Theta2) - K1 + K4*cos(Theta2) + K5;
    E = -2*sin(Theta2);
    F = K1 + (K4-1)*cos(Theta2) + K5;
    
    NumPlotPoints = 1000;
    Theta2List = linspace(0, 2*pi, NumPlotPoints);
    % Generate Linkage Plot
    N = length(Theta2List);
    for kk = 1:N
        Theta2 = Theta2List(kk);

    end
    
    % Calc for Theta4
    if(rc == 0)
        disp('Error, change rc value.');
        return;
    end
    
    discr4 = B^2 - 4*A*C;
    if(discr4 < 0)
        disp('Error: Discriminant for Theta4 Solution is Negative.')
        return;
    end

    FactorSign = PlusMinus;
    arg4num = -B + FactorSign * sqrt(discr4);
    arg4den = 2 * A;
    Theta4 = arg4num / arg4den;

    tmpAngle = NortonAtan2(arg4num, arg4den);
    discr4 = 2 * tmpAngle;
    if(discr4 > (2*pi));
        discr4 = discr4 - 2 * pi;
    end

    % Calc for Theta3 eq. 4.13
    discr3 = E^2 - 4*D*F;
    if(discr3 < 0)
        disp('Error: Discriminant for Theta3 Solution is Negative.')
        return;
    end

    FactorSign = PlusMinus;
    arg3num = -E + FactorSign * sqrt(discr3);
    arg3den = 2 * D;
    Theta3 = arg3num / arg3den;
    
    tmpAngle = NortonAtan2(arg3num, arg3den);
    discr3 = 2 * tmpAngle;
    if(discr3 > (2*pi));
        discr3 = discr3 - 2 * pi;
    end
    
% Allocate Storage
val_Theta2 = zeros(NumPlotPoints, 1);
val_Theta3 = zeros(NumPlotPoints, 1);
val_Theta4 = zeros(NumPlotPoints, 1);
val_Omega3 = zeros(NumPlotPoints, 1);
val_Omega4 = zeros(NumPlotPoints, 1);
val_Alpha3 = zeros(NumPlotPoints, 1);
val_Alpha4 = zeros(NumPlotPoints, 1);
% val_CouplerPointX = zeros(NumPlotPoints, 1);
% val_CouplerPointY = zeros(NumPlotPoints, 1);
val_vA = zeros(NumPlotPoints, 1);
val_vB = zeros(NumPlotPoints, 1);
val_vBA = zeros(NumPlotPoints, 1);
val_aA = zeros(NumPlotPoints, 1);
val_aB = zeros(NumPlotPoints, 1);
val_aBA = zeros(NumPlotPoints, 1);
    
    val_Theta2(kk) = Theta2;
    val_Theta3(kk) = Theta3;
    val_Theta4(kk) = Theta4;
    
Calc4BarVelocity(a, b, c, Theta2, Theta3, Theta4, Omega2);

        function[Omega3, Omega4, vA, vAx, vAy, vB, vBx, vBy, vBA, vBAx, vBAy, ...
             rc] = Calc4BarVelocity(a, b, c, Theta2, Theta3, Theta4, Omega2)
             %% Velocity Calculations
    
              % Equations to find Omega3 and Omega4, eq. 6.18
                Omega3 = (a*Omega2 * sin(Theta4 - Theta2))/((b* sin(Theta3 - Theta4)));
                Omega4 = (a*Omega2 * sin(Theta2 - Theta3))/((c* sin(Theta4 - Theta3)));

               % Equations to find linear velocities of points A and B
                vA = a * Omega2*(-sin(Theta2) + cos(Theta2));
                vBA = b * Omega3 * (-sin(Theta3) + cos(Theta3));
                vB = c * Omega4 * (-sin(Theta4) + cos(Theta4));

    
                 val_vA(kk) = vA;
                 val_vB(kk) = vB;
                 val_vBA(kk) = vBA;
                 val_Omega3(kk) = Omega3;
                 val_Omega4(kk) = Omega4;
                 
                 Calc4BarAccel(a, b, c, Theta2, Theta3, Theta4, Omega2, ...
                    Omega3, Omega4, Alpha2)
    
                 function[Alpha3, Alpha4, aA, aAx, aAy, aB, aBx, aBy, aBA, aBAx, aBAy, ...
                      rc] = Calc4BarAccel(a, b, c, Theta2, Theta3, Theta4, Omega2, ...
                      Omega3, Omega4, Alpha2)
            %% Acceleration Calculations
            %rc = rcAcc
  
            %if(rcAcc == 0)
            %    disp('Error: Faulty RC for acceleration calc.');
            %    return;
            %end
    
            % eq. 7.12c
            A = c*sin(Theta4);
            B = b*sin(Theta3);
            C = a*Alpha2*sin(Theta2) + a * (Omega2^2)*cos(Theta2) + b * ...
                (Omega3^2)*cos(Theta3) - c*(Omega4^2)*cos(Theta4);
            D = c * cos(Theta4);
            E = b * cos(Theta3);
            F = a*Alpha2*cos(Theta2) - a * (Omega2^2)*sin(Theta2) - b * ...
                (Omega3^2)*sin(Theta3) + c*(Omega4^2)*sin(Theta4);
       
            %eq. 7.12a & b
            Alpha3 = ((C*D - A*F)/(A*E - B*D));
            Alpha4 = ((C*E - B*F)/(A*E - B*D));
    
            % Euler Identity Equations
            aA = a*Alpha2*(-sin(Theta2) + cos(Theta2)) - a*(Omega2^2)*...
                (cos(Theta2) + sin(Theta2));
            aAx = -a*Alpha2*sin(Theta2) - a*(Omega2^2)*cos(Theta2);
            aAy = a*Alpha2*cos(Theta2) - a*(Omega2^2)*sin(Theta2);
            aBA = b*Alpha3*(-sin(Theta3) + cos(Theta3)) - b*(Omega3^2)*...
                (cos(Theta3) + sin(Theta3));
            aBAx = -b*Alpha3*sin(Theta3) - b*(Omega3^2)*cos(Theta3);
            aBAy = b*Alpha3*cos(Theta3) - b*(Omega3^2)*sin(Theta3);
            aB = c*Alpha4*(-sin(Theta4) + cos(Theta4)) - c*(Omega4^2)*...
                (cos(Theta4) + sin(Theta4));
            aBx = -c*Alpha4*sin(Theta4) - c*(Omega4^2)*cos(Theta4);
            aBy = c*Alpha4*cos(Theta4) - c*(Omega4^2)*sin(Theta4);
    
            val_aA(kk) = aA;
            val_aAx(kk) = aAx;
            val_aAy(kk) = aAy;
            val_aB(kk) = aB;
            val_aBx(kk) = aBx;
            val_aBy(kk) = aBy;
            val_aBA(kk) = aBA;
            val_aBAx(kk) = aBAx;
            val_aBAy(kk) = aBAy;
            val_Alpha3(kk) = Alpha3;
            val_Alpha4(kk) = Alpha4;
            end
        end
    end

    function [AngleRad] = NortonAtan2(y,x)
    %% Norton Function

    %Calc AngleRad and Return 0 to 2*pi per Norton 4-bar Convention

    % Note: the order of the argument list: y, then x

    AngleRad = atan2(y,x);

    %Convert Negative angle to Positive angle clockwise
    if(AngleRad < 0)
    AngleRad = AngleRad + 2*pi;
    end
    end


%% Generate Plots

figure(1);
plot(val_Theta2 * RadToDeg, val_Theta3 * RadToDeg, 'r');
hold on;
plot(val_Theta2 * RadToDeg, val_Theta4 * RadToDeg, 'b');
grid on;

figure(2);
plot(val_Theta2 * RadToDeg, val_Omega3, 'r');
hold on
plot(val_Theta2 * RadToDeg, val_Omega4, 'b');
grid on;

figure(3);
plot(val_Theta2 * RadToDeg, val_vA, 'r');
hold on
plot(val_Theta2 * RadToDeg, val_vA, 'b');
plot(val_Theta2 * RadToDeg, val_vBA, 'k');
grid on;

figure(4);
plot(val_Theta2 * RadToDeg, val_aA, 'r');
hold on
plot(val_Theta2 * RadToDeg, val_aB, 'b');
plot(val_Theta2 * RadToDeg, val_aBA, 'k');
xlabel( 'Theta2 (deg)' );
ylabel( 'Acceleration (in/sec^2)');
title( 'Linkage Acceleration' );
grid on;

end

