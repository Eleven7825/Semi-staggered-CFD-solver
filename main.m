Nx = 100;
Ny = 100;
Ulid = 1;
Re = 3000;
Lx = 1;
Ly = 1;
dx = Lx/Nx;
dy = Ly/Ny;
u = zeros(Ny+2, Nx+2); v= u; 
p = u;Dp = p;
U = zeros(Ny, Nx+1);
V = zeros(Ny+1, Nx);
dt = 0.005;
ipx = 2:(Nx+1); % collocated; phyiscal
ipy = 2:(Ny+1); % collocated; phyiscal
us = u(ipy, ipx); % u star; physical
vs = v(ipy, ipx); % v star; physical
Us = U; Vs = V;   % Us and Vs for boundaries
u(Ny+2, :) = 2 * Ulid;
Nmax = 30000;

% setup for poisson
kx = (0 : Nx-1);
ky = (0 : Ny-1);
[kx,ky] = meshgrid(kx,ky);
mwx = 2 * (cos(pi * kx / Nx)-1)/dx^2;
mwy = 2 * (cos(pi * ky / Ny)-1)/dy^2;

for i = 1:Nmax
%% for the first time step:
% interplation from u->U
u_w = (u(ipy, ipx-1)+u(ipy, ipx))/2; u_e = (u(ipy, ipx+1)+u(ipy, ipx))/2; 
u_n = (u(ipy+1, ipx)+u(ipy, ipx))/2; u_s = (u(ipy-1, ipx)+u(ipy, ipx))/2; 

v_n = (v(ipy+1, ipx)+v(ipy, ipx))/2; v_s = (v(ipy-1, ipx)+v(ipy, ipx))/2; 
v_w = (v(ipy, ipx-1)+v(ipy, ipx))/2; v_e = (v(ipy, ipx+1)+v(ipy, ipx))/2; 

conu = (U(:, 2:Nx+1).*u_e - U(:, 1:Nx).*u_w)/dx + ...
       (V(2:Ny+1,: ).*u_n - V(1:Ny, :).*u_s)/dy;

conv = (U(:, 2:Nx+1).*v_e - U(:, 1:Nx).*v_w)/dx + ...
       (V(2:Ny+1,: ).*v_n - V(1:Ny, :).*v_s)/dy;

diffu = (u(ipy+1,ipx)+u(ipy-1,ipx)-2*u(ipy,ipx))/Re/dy^2+...
        (u(ipy,ipx+1)+u(ipy,ipx-1)-2*u(ipy,ipx))/Re/dx^2;

diffv = (v(ipy+1,ipx)+v(ipy-1,ipx)-2*v(ipy,ipx))/Re/dy^2+...
        (v(ipy,ipx+1)+v(ipy,ipx-1)-2*v(ipy,ipx))/Re/dx^2;

us = u(ipy,ipx) + dt * (-conu + diffu);
vs = v(ipy,ipx) + dt * (-conv + diffv);

% interpolation
Us(:, 2 : Nx) = (us(:,1:Nx-1)+us(:,2:Nx))/2;
Vs(2 : Ny, :) = (vs(1:Ny-1,:)+vs(2:Ny,:))/2;

F = ((Us(:,2:Nx+1) - Us(:,1:Nx))/dx + (Vs(2:Ny+1,:) - Vs(1:Ny,:))/dy)/dt;

%% Solve the poisson equation
B = dct2(F);
A =  B ./ ( mwx+mwy ); A(1,1) = 0;
p(ipy, ipx) = idct2(A);

% BC for p:
p(:, 1) = p(:,2);
p(:,Nx+2) = p(:, Nx+1);
p(1,:) = p(2,:);
p(Ny+2, :) = p(Ny+1,:);

%% The third step:
  % for cell centered values:
u(ipy, ipx) = us - dt * (p(ipy,ipx+1)-p(ipy,ipx-1))/(2*dx);
v(ipy, ipx) = vs - dt * (p(ipy+1,ipx)-p(ipy-1,ipx))/(2*dy);
  % for face centered values:
U(:, 2:Nx) =   Us(:, 2:Nx) -dt * (p(ipy, 3:Nx+1) - p(ipy, 2:Nx))/dx;
V(2:Ny, :) =   Vs(2:Ny, :) -dt * (p(3:Ny+1,ipx) -  p(2:Ny,ipx))/dy;

% Impose boundary for u:
u(1,:) = -u(2,:);
u(:, Nx+2) = -u(:, Nx+1);
u(:, 1) = -u(:, 2);
u(Ny+2, :) = 2 * Ulid-u(Ny+1,:);

% Impose boundary for v:
v(:, 1) = -v(:, 1);
v(:,Nx+2) = -v(:,Nx+1);
v(1,:) = -v(2,:);
v(Ny+2, :) = -v(Ny+1, :);

if mod(i,50)==0
    vortamp = log(abs((u(ipy+1, ipx)-u(ipy-1, ipx))/(2*dy)-...
       (v(ipy, ipx+1)-v(ipy, ipx-1))/(2*dx)));
    vortamp(vortamp<-9) = -9;
    contourf(vortamp,-19:.5:4, 'LineColor','none')
    clim([-19 4])
    drawnow; 
end
end

% 
% vort = (u(ipy+1, ipx)-u(ipy-1, ipx))/(2*dy)-...
%        (v(ipy, ipx+1)-v(ipy, ipx-1))/(2*dx);
% 
% div = (v(ipy+1, ipx)-v(ipy-1, ipx))/(2*dy)+...
%       (u(ipy, ipx+1)-u(ipy, ipx-1))/(2*dx);
% 
% Div = (V(2:Ny+1,:)-V(1:Ny,:))/dy + (U(:,2:Nx+1)-U(:,1:Nx))/dx;
% 
% uabs = sqrt(u(ipy,ipx).^2+v(ipy,ipx).^2);
figure
plot(u(ipy,Ny/2),1:Nx)
figure
plot(1:Ny,v(Nx/2,ipx))
% 
