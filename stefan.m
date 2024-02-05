clear all
close all
% function [u,h,t] = stefan(u0,h0,St,tm,Tmax)
% Solving Stefan problem in dimensionless form:
% u_t = u_xx for u(x,t) on 0 < x <= h(t) 
% u(x,t) = tm for u(x,t) on h(t) < x < 1
% BC: u(0,t) = 1 and u(h(t),t) = tm
% IC: u(x,0) = u0 and h(0) = h0
% dh(t)/dt = -St u_x at x = h(t), where h(0) = h0
% St: Stefan number = cT/L
% c: Specific heat capacity
% T: The reference temperature (ambient temerature)
% L: Latent heat of melting 

% Number of grid point
n = 200;
% Initial layer thickness
h0 = 0.01;
% Initial temperature distribution
u0 = zeros(1,n); 
for i=1:n
    u0(i) = 0.0;
end
% Stefan number 
St = 0.5;
% Melting temperature
tm = 0.0;
% Maximum time
Tmax = 10;

if (h0 <= 0 || h0 > 1) 
    error('Moving boundary is outside domain'), 
end

cfl = 10;
nx = length(u0);
x = linspace (0 ,1 ,nx); 
dx = x(2)-x(1); % Spacial step size
dt = cfl*dx^2; % Time step size
t = 0:dt:Tmax; % Number of time steps
nt = length(t);
c = dt/dx^2;

% Initialize moving boundary position
h = zeros(1,nt); h(1) = h0;
ih = floor (h0/dx)+1;
% Initialize temperature
u = zeros(nx,nt); u(:,1) = u0;

for j=2:nt
% Update h and find lower bound on grid
h(j) = h(j-1) - dt * St/dx * (u(ih+1,j-1) - u(ih,j-1)); 
ih = floor(h(j)/dx)+1;
% If ih < 1, all the ice is melted; if ih >= nx, everything is frozen
 if (ih < 1 || ih > nx-1)
   disp('Stop the computation − the entire domain is melted/frozen') 
  break
 end
 
% Make A for implicit (backward Euler) method
Ad = ones(1,nx);
Ad(2:ih) = 2*c+1;
Asub = zeros(1,nx-1);
Asub(1:ih-1) = -c;
Asup = zeros(1,nx-1);
Asup(2:ih) = -c;
A = diag(Ad) + diag(Asub,-1) + diag(Asup,+1);
% Make b right hand side
b = ones(nx,1);
b(2:ih) = u(2:ih,j-1); % Edge of ice is at ih 
b(ih+1:end) = tm;
% Solve for Au = b at current time
u(:,j) = A\b; 
end

figure
subplot (2,1,1)
pcolor(t,x,u), shading interp
colorbar
xlabel('t') 
ylabel('x') 
set(gca,'Fontsize',14) 
title('u(x,t)')
axis([0 1 0 1])

subplot (2,1,2) 
plot(t,h) 
hold on
plot(t,sqrt(2*t*St)) % approximated solution for st << 1 
grid
axis([0 1 0 1]) 
legend({'Numerical, St=0.5','Analytical approximation, St<<1'},'Location','southeast')
legend('boxoff')
xlabel('t')
ylabel('h')
title ('h(t)')
set(gca,'FontSize',14)