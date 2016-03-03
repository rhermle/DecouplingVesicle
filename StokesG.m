clear all
close all
clc

%2D Stokes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pressure Differential (Right side is 100)
p0 = 0;

%Viscosity
mu = 1;

%grav. constant
g = 9.8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
width = 20;
height = 20;
R = 5;
L = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute Stokes2DG with M=50 and use it as the "analytical" solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num = 50;
%[P U V X Y] = Stokes2DG(g, num, p0, mu, 1);

%Compute the L2E with 20 different grid spacing values
testPoints = num-20:2:num;

L2EP = zeros(length(testPoints),1);
L2EU = zeros(length(testPoints),1);
L2EV = zeros(length(testPoints),1);
d = zeros(length(testPoints),1);

% Compute the L2 error
%2D: || u(x,y) || = sqrt( 1/M^2 sum_{j=1}^M sum_{k=1]^M u_{j,k}^2 }
for i = 1:length(testPoints)
    
    [ p u v x y] = Stokes2DG(width, height, R, L, g, testPoints(i), p0, mu, 0);
    
    d(i) = y(2,1) - y(1,1);
    
    [P U V] = pTest(x,y,R,L);
        
    L2EP(i) = L2EP(i) + sum(sum((p - P).^2));
    L2EU(i) = L2EU(i) + sum(sum((u - U).^2));
    L2EV(i) = L2EV(i) + sum(sum((v - V).^2));  
    
    L2EP(i) = sqrt(L2EP(i) / (testPoints(i)^2));
    L2EU(i) = sqrt(L2EU(i) / (testPoints(i)^2));
    L2EV(i) = sqrt(L2EV(i) / (testPoints(i)^2));
end

figure()
subplot(2,2,1);
loglog(d,L2EP,'-',d,d.^2,'--');
title('L2 Error for P (Pressure)');
xlabel('dx');
ylabel('L2 Error');

subplot(2,2,2);
loglog(d,L2EU,'-',d,d.^2,'--');
title('L2 Error for U (Horizontal Velocity)');
xlabel('dx');
ylabel('L2 Error');

subplot(2,2,3);
loglog(d,L2EV,'-',d,d.^2,'--');
title('L2 Error for V (Vertical Velocity)');
xlabel('dx');
ylabel('L2 Error');

figure()
surf(x,y,u);
title('U');
xlabel('X');
ylabel('Y');

figure()
surf(x,y,v);
title('V');
xlabel('X');
ylabel('Y');

figure()
surf(x,y,p);
title('P');
xlabel('X');
ylabel('Y');

figure();
quiver(x(1:4:end,1:4:end),y(1:4:end,1:4:end),u(1:4:end,1:4:end),v(1:4:end,1:4:end));
streamline(x,y,u,v,1.0,0.7);

%%%%%%%%%%%%%benchmark%%%%%%%%%%%%%%%%
REPS = 10;
num = [25 50 75 100 200];
tic;
for i = 1:length(num)
    for j = 1:REPS
        j
        [ p u v x y] = Stokes2DG(width, height, R, L, g, num(i), p0, mu, 0);
    end
    averageTime(i) = toc / REPS;
end

averageTime

figure()
plot(num,averageTime);
title('Execution Time vs. Number of Discretization Points');
xlabel('M');
ylabel('Execution Time (s)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

