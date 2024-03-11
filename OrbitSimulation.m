%%============
% This matlab program simulate the part of Journal,
% "New Approach to Satellite Formation-Keeping: Exact Solution to the Full Nonlinear Problem"
% by Hancheol Cho and Adam Yu, ASCE/OCTOBER 2009
% DOI: 10.1061/(ASCE) AS.1943-5525.0000013
%%============

clc; clear;

%% parameters
G = 6.67430e-11;  % J*m/kg^2 Gravitational constant
M = 5.9722e24; % kg, Earth Mass
rho = 50000;
w0 = 0;
i = deg2rad(80);
ohm = deg2rad(30);
r0 = 7e6; % radius of chief
r0vector = [r0;0;0];

%% time
n = sqrt(G*M/(r0^3)); % rad/s
p = 2*pi/n; % period
tspan = [0 p*3];  % Time span for the simulation

%% Initial conditions
% x0 y0 z0 xdot0 ydot0 zdot0 직접 적분해서 대입
alpha0 = 0.55;
x0 = rho*sin(alpha0)/2;                                       
y0 = rho*cos(alpha0);                                               
z0 = rho*sin(alpha0);  
xdot0 = rho*n*cos(alpha0)/2;                                          
ydot0 = -rho*n*sin(alpha0);  
zdot0 = rho*n*cos(alpha0); 

xyz0 = [x0;y0;z0];
xyzdot0 = [xdot0;ydot0;zdot0];

hillxyz0 = xyz0 + r0vector;

%% rotation matrix
syms nval time ival ohmval real
t0 = 0;
wval = w0 + nval*(time-t0);

R3w = [cos(wval) sin(wval) 0; -sin(wval) cos(wval) 0; 0 0 1];
R1i = [1 0 0; 0 cos(ival) sin(ival); 0 -sin(ival) cos(ival)];
R3ohm = [cos(ohmval) sin(ohmval) 0; -sin(ohmval) cos(ohmval) 0; 0 0 1];

R_cal = R3w * (R1i * R3ohm);

% hill x y z -> ECI X Y Z
R_0 = subs(R_cal, [nval, ival, ohmval, time], [n, i, ohm, 0]);
ECIXYZ0 = inv(R_0) * hillxyz0;

% hill x' y' z' -> ECI X' Y' Z'
Rdot = diff(R_cal, time);
Rdot0 = subs(Rdot, [nval, ival, ohmval, time], [n, i, ohm, 0]);
ECIXYZdot0 = inv(R_0) * (xyzdot0 - Rdot0 * ECIXYZ0);

initial_state = double(vpa([ECIXYZ0; ECIXYZdot0], 12)); %소수점 아래 12자리

%% Solve
opts = odeset('RelTol',1e-12,'AbsTol',1e-14);
[t, result] = ode45(@(t, yval) constrained_ECI(t, yval, G, M, i, ohm, n), tspan, initial_state, opts);

X = result(:, 1);
Y = result(:, 2);
Z = result(:, 3);
ECImatrix = [X Y Z];

N = length(t);

%% Force
Force = zeros(N,3);
mag_Force = zeros(N, 1);
for indexnum1 = 1:N
    time = t(indexnum1);
    [~, Force(indexnum1, 1:3)] = constrained_ECI(time, result(indexnum1,:)', G, M, i, ohm, n);
    mag_Force(indexnum1, 1) = norm(Force(indexnum1, 1:3));
end

%% eci to hill
hill = zeros(N:3);
for indexnum2 = 1:N
    time = t(indexnum2);
    R3w = [cos(n*time) sin(n*time) 0; 
        -sin(n*time) cos(n*time) 0; 
        0 0 1];
    R1i = [1 0 0; 
        0 cos(i) sin(i); 
        0 -sin(i) cos(i)];
    R3ohm = [cos(ohm) sin(ohm) 0; 
        -sin(ohm) cos(ohm) 0; 
        0 0 1];

    R = R3w * (R1i * R3ohm);

    element = ECImatrix(indexnum2, :);
    hilldata = (R * element') - r0vector;
    hill = [hill hilldata];
end

% hill x y z
x = hill(1, :);
y = hill(2, :);
z = hill(3, :);

%% plot Force
period = t./p; 
figure();
subplot(2, 2, 1);
plot(period, Force(:, 1));
xlim([0, 3]);
ylim([-1e-3, 1e-3]);
grid on;

subplot(2,2,2);
plot(period, Force(:,2));
xlim([0, 3]);
ylim([-5e-4, 15e-4]);
grid on;

subplot(2,2,3);
plot(period, Force(:,3));
xlim([0, 3]);
ylim([-1e-3, 1e-3]);
grid on;

subplot(2,2,4);
plot(period, mag_Force);
xlim([0, 3]);
ylim([0.5e-3, 1.5e-3]);
grid on;

%% Plot orbit
figure();
subplot(1, 2, 1);
plot(x./rho, z./rho);hold on;
plot(x(1)/rho, z(1)/rho, 'o')
title('Constrained Motion in the xz-plane');
xlabel('x/rho');
ylabel('z/rho');
grid on;

subplot(1, 2, 2)
plot(y./rho, z./rho);hold on;
plot(y(1)/rho, z(1)/rho, 'o')
xlim([-1 1])
title('Constrained Motion in the yz-plane');
xlabel('y/rho');
ylabel('z/rho');
grid on;

function [state, Force] = constrained_ECI(t, yval, G, M, ival, ohmval, nval)
    %% mass of deputy
    m = 1;

    vel = [yval(4);yval(5);yval(6)]; % X dot
    r = [yval(1);yval(2);yval(3)]; % X

    R_bar       = zeros(2,3);
    R_bar_dot   = zeros(2,3);
    R_bar_ddot  = zeros(2,3);
    R_tilde     = zeros(2,3);
    R_tilde_dot = zeros(2,3); 
    R_tilde_ddot= zeros(2,3);

    R_bar = [- sin(t*nval)*cos(ohmval) - cos(t*nval)*cos(ival)*sin(ohmval), cos(t*nval)*cos(ival)*cos(ohmval) - sin(t*nval)*sin(ohmval), cos(t*nval)*sin(ival); ...
        sin(ival)*sin(ohmval), -cos(ohmval)*sin(ival), cos(ival)];
    R_bar_dot = [-nval*cos(ohmval)*cos(t*nval)+nval*sin(ohmval)*cos(ival)*sin(t*nval), -nval*cos(ohmval)*cos(ival)*sin((t*nval))-nval*sin(ohmval)*cos((t*nval)), -nval*sin(ival)*sin(t*nval); ...
        0 0 0];
    R_bar_ddot = [(nval^2)*cos(ohmval)*sin(t*nval)+(nval^2)*sin(ohmval)*cos(ival)*cos(t*nval), -(nval^2)*cos(ohmval)*cos(ival)*cos(t*nval)+(nval^2)*sin(ohmval)*sin(t*nval), -(nval^2)*sin(ival)*cos(t*nval); ...
        0 0 0];
    
    R_tilde = [cos(t*nval)*cos(ohmval) - sin(t*nval)*cos(ival)*sin(ohmval), cos(t*nval)*sin(ohmval) + sin(t*nval)*cos(ival)*cos(ohmval), sin(t*nval)*sin(ival); ...
        sin(ival)*sin(ohmval), -cos(ohmval)*sin(ival), cos(ival)];
    R_tilde_dot = [-nval*cos(ohmval)*sin(t*nval)-nval*sin(ohmval)*cos(ival)*cos(t*nval), -nval*sin(ohmval)*sin(t*nval)+nval*cos(ohmval)*cos(ival)*cos(t*nval), nval*sin(ival)*cos(t*nval); ...
        0 0 0];
    R_tilde_ddot = [-(nval^2)*cos(ohmval)*cos(t*nval)+(nval^2)*sin(ohmval)*cos(ival)*sin(t*nval), -(nval^2)*sin(ohmval)*cos(t*nval)-(nval^2)*cos(ohmval)*cos(ival)*sin(t*nval), -(nval^2)*sin(ival)*sin(t*nval); ...
        0 0 0];

    %% [A11 A12 A13]
    Arow1 = (r')*((R_bar')*R_bar);
    A11 = Arow1(1);
    A12 = Arow1(2);
    A13 = Arow1(3);

    %% [A21 A22 A23]
    Arow2 = [2, -1]*R_tilde;
    A21 = Arow2(1);
    A22 = Arow2(2);
    A23 = Arow2(3);
    Amatrix = [A11 A12 A13;A21 A22 A23];

    %% [b1 b2]
    b1 = -(r')*(R_bar')*R_bar_ddot*(r) - 2*(r')*(R_bar')*R_bar_dot*(vel)...
        -(r')*(R_bar_dot')*R_bar_dot*(r) - 2*(r')*(R_bar_dot')*R_bar*(vel)...
        -(vel')*(R_bar')*R_bar*(vel);
    b2 = -[2 -1]*R_tilde_ddot*(r) - [4 -2]*R_tilde_dot*(vel);
    Bmatrix = [b1;b2];

    %% Moore-Penrose inverse
    A1MPinverse  = (1/(A11^2+A12^2+A13^2))*(Amatrix(1,:)');
    Bbeta = (A11*A21 + A12*A22 + A13*A23)/(A11^2+A12^2+A13^2);
    cMPinverse = (1/(A21^2+A22^2+A23^2-Bbeta^2*(A11^2 + A12^2+A13^2)))*((Amatrix(2,:)')-Bbeta*(Amatrix(1,:)'));

    AMPinverse = [(A1MPinverse - Bbeta*cMPinverse) cMPinverse];

    rnorm = sqrt(r(1)^2+r(2)^2+r(3)^2);

    unconstrained_accel = -(G*M/rnorm^3)*r;
    constrained_accel = (unconstrained_accel) ...
        + (AMPinverse*Bmatrix) ...
        + ((G*M/rnorm^3)*AMPinverse*Amatrix*r);

    Forcetest = m*(constrained_accel - unconstrained_accel);
    Force = Forcetest';

    state = [vel;constrained_accel];

end