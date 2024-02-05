%---
clear all
close all
hold on
%---

NLR=1; % number of layers
L=2.0;  % save length
Dt=0.01;  % time step
Nstep=128*2*128; % number of steps
th0=0.*pi;  % inclination angle
gac=1.0; % acceleration of gravity
NSG=2*16; % number of divisions
ICU=2;  % backward differences
ICU=1;  % central differences

mu(1)=1.0;
rho(1)=1.0;
gamma(1)=0.000;

%rho(2)=0.5;
%mu(2)=1.0;
%gamma(2)=0.00;

%---
% prepare
%---

rho(NLR+1)=0.0;

cs0 = cos(th0);
sn0 = sin(th0);
gx =  gac*sn0;
gy = -gac*cs0;

ROT = [cs0,sn0;-sn0,cs0];

%---
% wall and initial profiles
%---

Dx=L/NSG;

for k=1:NSG+1
 x(k)=(k-1)*Dx;
 arg=2*pi*x(k)/L;
 wall(k)=0.0;
 y(k,1)=0.2+0.1*cos(arg);
 if(NLR>=2)
 y(k,2)=0.4+0.1*cos(arg);
 end
 if(NLR>=3)
 y(k,3)=0.6+0.1*cos(arg);
 end
end

for step=1:Nstep

 dydt =films_pde ...
 ...
   (NLR,NSG,Dx,wall,y,rho,mu,gamma,gx,gy,ICU);

 ysave = y;
 dydtsave = dydt;
 y = y+dydt*Dt;
 dydt =films_pde ...
 ...
   (NLR,NSG,Dx,wall,y,rho,mu,gamma,gx,gy,ICU);
 y = ysave+0.5*(dydt+dydtsave)*Dt;

 for k=1:NSG+1
   xx = ROT*[x(k) wall(k)]';
   xplot0(k)=xx(1)/L;
   yplot0(k)=xx(2)/L;
   xx = ROT*[x(k) y(k,1)]';
   xplot1(k)=xx(1)/L;
   yplot1(k)=xx(2)/L;
   if(NLR>=2)
   xx = ROT*[x(k) y(k,2)]';
   xplot2(k)=xx(1)/L;
   yplot2(k)=xx(2)/L;
   end
   if(NLR>=3)
   xx = ROT*[x(k) y(k,3)]';
   xplot3(k)=xx(1)/L;
   yplot3(k)=xx(2)/L;
   end
 end

if(step==1 || mod(step,500)==0)
  Handle0 = plot(xplot0,yplot0,'k');
  hold on
  Handle1 = plot(xplot1,yplot1,'r.-');
  if(NLR>=2)
  Handle2 = plot(xplot2,yplot2,'r.-');
  end
  if(NLR>=3)
  Handle3 = plot(xplot3,yplot3,'r.-');
  end
%  set(Handle0,'EraseMode','xor')
 % set(Handle1,'EraseMode','xor')
 % set(Handle2,'EraseMode','xor')
 % set(Handle3,'EraseMode','xor')
  xlabel('x/L','fontsize',15)
  ylabel('y/L','fontsize',15)
  set(gca,'fontsize',15)
  axis equal
  %axis([0 1.2*cs0 -sn0 sn0])
  box
end
 set(Handle0,'XData',xplot0,'YData',yplot0)
 set(Handle1,'XData',xplot1,'YData',yplot1)
 if(NLR>=2)
 set(Handle2,'XData',xplot2,'YData',yplot2)
 end
 if(NLR>=3)
 set(Handle3,'XData',xplot3,'YData',yplot3)
 end
 drawnow


end
