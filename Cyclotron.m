%% Introduction to electrodynamics. Dynamic particle simulation.
close all
clear all
clc
%% General parameters
R = 0.05; % Radius (m)
Ds = 0.005; % Gap between D's (m)
B = [0,0,1]; % Magnetic field (T)
V = 5000; % Potential diff (V)
m = 1.6726219*10^-27; % Mass (kg)(proton)
q =  1.602176634*10^-19; % Charge (C)(proton)
dt = 10^-11; % Delta time (seconds)
tf = 4e-06; % Final time (s) 0.7854 microseconds
t = [0:dt:tf]; % Time (s)
w = q*norm(B)/m; % Frequency
% Oscillating electric field used asumming B in all space.
E = (-V/Ds)*sin(w.*t); 
% Oscillating electric field used assuming B just in the D's.
V = @(t) 5000/Ds*cos(t.*w.*(1+exp(-t./1e-7)));
c = 2.998e+8; % Light speed


%% Initial vectors
r = [0 -0.005 0]; % Position
rx = r(1); % Position x
ry = r(2); % Position y
rz = r(3); % Position z
v = [0 0 0]; % Velocity
KE = [0]; % Kinetic energy

%% Cyclotron
loops = 0;
sig = -1;
for i = 1:length(t)
    % Particle in the gap -> Adapting electric force.
    if r(1) < Ds/2 && r(1) > -Ds/2
        % Case 1: Asumming B in all space.
        f = [q*E(i),0,0] + q*cross(v,B); % Lorentz force
        % Case 2: Assuming B just into the D's.
%         DP = V(t(i));
%         f = [q*DP,0,0]; % Electric force
    else 
    % Particle into the D's -> Adapting magnetic force.
        f = q*cross(v,B); % Magnetic force
    end 
    v = v + f*dt/m; % New velocity
    r = r + dt*v; % New position
    rx(i+1) = r(1); % x component
    ry(i+1) = r(2); % y component
    rz(i+1) = r(3); % z component
    KE(i+1) = 0.5*m*(norm(v)^2); % New kinetic energy
    if rx(i+1) > 0 && sig < 0 % Lap counter
        loops = loops + 1;
        sig = 1;
    elseif rx(i+1) < 0
        sig = -1;
    end
    % Break statement. If the radius of the trajectory is larger than R of
    % the cyclotron.
    if norm(r)<R 
        continue
    else
        break
    end
end

%% Cyclotron plot
figure(1)
th = linspace(pi/2, -pi/2, 100);
xR = R*cos(th)+Ds/2;
yR = R*sin(th);
plot(xR,yR,'r','LineWidth',2);
axis equal;
hold on
plot([Ds/2 Ds/2], [-R R],'r','LineWidth',2)
xL = -(R*cos(th)+Ds/2) ;
yL = -(R*sin(th)) ;
hold on
plot(xL,yL,'r','LineWidth',2)
hold on
plot([-Ds/2 -Ds/2], [-R R],'r','LineWidth',2)
xlim([-R-2*Ds R+2*Ds])
ylim([-R-2*Ds R+2*Ds])
hold on

%% Magnetic field plot
x = -R:Ds:R;
y = -R:Ds:R;
[X,Y] = meshgrid(x,y);
plot(X,Y,'b.')
title("Bird's eye view of ciyclotron")
xlabel('Z')
ylabel('X')
hold on

% x2 = Ds:Ds:R;
% y2 = -R:Ds:R;
% [X2,Y2] = meshgrid(x2,y2);
% plot(X2,Y2,'b.')
%% Plot of the particle in certain t.
plot3(rx',ry',rz','b')
title("Particle trajectory")
xlabel('Z')
ylabel('X')
axis square
grid on

%% Kinetic energy analysis
figure(2)
eV = 1.60218*10^-19;
Ekv = KE/eV;
tEk = t(1,1:length(KE));
plot(tEk,Ekv,'r','LineWidth',1);
title('Kinetic energy of the proton')
xlabel('Time (seconds)')
ylabel('Kinetic Energy (eV)')
legend('Ek')
grid on

fprintf('The kinetic energy in Joules is: %2.10f \n',KE(end))
fprintf('The kinetic energy in eV is: %2.10f \n',Ekv(end))


%% Time
Time = length(rx)*dt; % Total time
fprintf('The elapsed time was: %d \n',Time)

%% Final velocity
finalv = norm(v);
perOfC = finalv/c*100;
fprintf('The velocity is %d m/s \n',finalv)
fprintf('The velocity is %f of c \n',perOfC)



