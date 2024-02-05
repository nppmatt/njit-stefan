function dydt =films_pde ...
...
   (NLR,NSG,Dx,wall,y,rho,mu,gamma,gx,gy,ICU)

%---------------------------------------
% Compute the rate of change of y position
% of the interfaces:
% dy_i/dt, for i=1,...,NLR
%
% NLR: number of layers
% NSG: number of segments
% Dx: x finite-difference interval
% wall: wall profile
% y: position of interfaces
% ICU=1 for central differences
% ICU=2 for backward differences
%---------------------------------------

%---
% prepare
%---

Dx2 = 2.0*Dx;

%---
% compute derivatives
%---

for i=1:NLR  % run over interfaces

%---
% first derivative
%---

 y1(1,i)=(y(2,i)-y(NSG,i))/Dx2;
 for k=2:NSG
   y1(k,i) = (y(k+1,i)-y(k-1,i))/Dx2;
 end 
 y1(NSG+1,i) = y1(1,i);

%---
% second derivative
%---

 y2(1,i)=(y1(2,i)-y1(NSG,i))/Dx2;
 for k=2:NSG
   y2(k,i)=(y1(k+1,i)-y1(k-1,i))/Dx2;
 end
 y2(NSG+1,i) = y2(1,i);

%---
% third derivative
%---

 y3(1,i)=(y2(2,i)-y2(NSG,i))/Dx2;
 for k=2:NSG 
   y3(k,i) = (y2(k+1,i)-y2(k-1,i))/Dx2;
 end
 y3(NSG+1,i) = y3(1,i);

%---
end
%---


%--------------------------------------
% Compute the effective pressure drop G
% using equations (9.3.10) and (9.3.12)
%--------------------------------------

for i=1:NLR        % over interfaces
  for k=1:NSG       % over points
    dpdx = 0.0;
    for j=i:NLR
     dpdx= dpdx - gamma(j)*y3(k,j) ...
        + gy*(rho(j+1)-rho(j))*y1(k,j);
     end
      G(k,i) = 0.50*(-dpdx + rho(i)*gx)/mu(i);
  end
end

%----------------------------
% Compute the coefficients B 
%----------------------------

 for k=1:NSG
    B(k,NLR) = 2.0*G(k,NLR)*y(k,NLR);
 end

 for i=NLR-1:-1:1
  for k=1:NSG
    B(k,i) = 2.0*G(k,i)*y(k,i) ...
       +mu(i+1)/mu(i)*(B(k,i+1)-2.0*G(k,i+1)*y(k,i));
  end
 end

%----------------------------
% Compute the coefficients A 
%----------------------------

 for k=1:NSG
   A(k,1) = wall(k)*(wall(k)*G(k,1)-B(k,1));
 end

 for i=1:NLR-1
   for k=1:NSG
    A(k,i+1) = A(k,i) ...
         +(B(k,i)-B(k,i+1))*y(k,i) ...
         -(G(k,i)-G(k,i+1))*y(k,i)*y(k,i);
  end
 end

%--------------------------
% Compute the flow rates Q 
%--------------------------

for k=1:NSG
 Q(k,1)  = A(k,1)*(y(k,1)-wall(k)) ...
    +0.5*B(k,1)*(y(k,1)*y(k,1)-wall(k)*wall(k));
    -G(k,1)*(y(k,1)^3-wall(k)^3)/3.0;
end
Q(NSG+1,1) = Q(1,1);

for i=2:NLR
  for k=1:NSG
    Q(k,i) = A(k,i)*(y(k,i)-y(k,i-1)) ...
    +0.5*B(k,i)*(y(k,i)*y(k,i)-y(k,i-1)*y(k,i-1)) ...
    -G(k,i)*(y(k,i)^3-y(k,i-1)^3)/3.0;
  end
 Q(NSG+1,i) = Q(1,i);
end

%---
% compute dQ/dx
%
% by central differences (ICU =  1)
% or upwind  differences (ICU ne 1)
%---

%---
  if(ICU==1)   % central differences
%---

 for i=1:NLR
  dQdx(1,i)=(Q(2,i)-Q(NSG,i))/Dx2;
  for k=2:NSG
   dQdx(k,i)=(Q(k+1,i)-Q(k-1,i))/Dx2;
  end
 end

%---
 else   % backward differences
%---

 for i=1:NLR
  dQdx(1,i)=(Q(1,i)-Q(NSG,i))/Dx;
  for k=2:NSG
    dQdx(k,i)=(Q(k,i)-Q(k-1,i))/Dx;
  end
 end

%---
  end
%---

%---
% finally compute dydt
%---

 for i=1:NLR
   for k=1:NSG
   dydt(k,i) = 0.0;
   for j=1:i
     dydt(k,i) = dydt(k,i)-dQdx(k,j);
   end
   dydt(NSG+1,i) = dydt(1,i);
   end
 end


%---
% Done
%---

return
