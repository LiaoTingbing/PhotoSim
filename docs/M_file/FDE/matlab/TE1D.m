

clear;
L = 10e-6;
dx = 0.1e-6;
w = 2e-6;
index = 1.2;
nmodes = 12;
lambda = 1.55e-6;

%% 

eps0 = 8.854187817e-12;
mu0 = pi * 4e-7;
Z0 = sqrt(mu0/eps0);

x = (-L/2:dx:L/2)';
xm = x + 0.5 * dx ;
cx = round(L/dx);
nx = cx + 1;
k0 = 2*pi /lambda;
s1 = dsearchn(x , -w/2);
s2 = dsearchn(x , w/2 );
I = ones(nx,1);


eps_y = I;
eps_y(s1:s2) = index^2;
eps_y(s1) = 0.5*( 1+ index^2);
eps_y(s2) = 0.5*( 1+ index^2);
 

EPS_Y = spdiags(eps_y,0,nx,nx);
MU_X = speye(nx,nx);
MU_Z = speye(nx,nx);


%% TM

LXHZ = spdiags([-I,I]/dx , [-1 0] , nx, nx);
LXHZ(1,1)=0;
LXEY = spdiags([-I,I]/dx , [0 1] , nx,nx);
LXEY(end,:) = LXEY(end-1,:);

A = MU_X*LXHZ*inv(MU_Z)*LXEY + k0*k0*MU_X*EPS_Y;

[ey,d]=eigs(A,nmodes,index^2 * k0^2);

neff =  sqrt(diag(d))/k0;

%% 过滤
 
idx =  (abs(neff-1))<1e-9;
ey(:,idx)=[];
neff(idx)=[];
 

%% 其它场分量

for midx = 1:length(neff)
    beta = k0* neff(midx);
    hx(:,midx) = - 1j*beta/k0*inv(MU_X) *ey(:,midx);
end

hz =1/k0 * inv(MU_Z) * LXEY*ey;
hz = [hz(1,:) ; 0.5*(hz(1:nx-1,:)+hz(2:nx,:))];
 
hx = 1j*hx/Z0;
hz = 1j*hz/Z0;

%% 绘图

s = "images/TE/";
mkdir(s)

plot(x , eps_y);
xlabel("X (m)")
ylabel("介电常数")
title("EPS Y")
exportgraphics(gcf , s + "EPS_Y.png")



%% 

 
plot(x , abs(ey));
xlabel("X (m)")
ylabel("Amplititude")
title("Ey Amplititude")
exportgraphics(gcf , s + "Ey_Amplititude.png")
 

plot(x , abs(hx));
xlabel("X (m)")
ylabel("Amplititude")
title("Hx Amplititude")
exportgraphics(gcf , s + "Hx_Amplititude.png")


plot(x , abs(hz));
xlabel("X (m)")
ylabel("Amplititude")
title("Hz Amplititude")
exportgraphics(gcf , s + "Hz_Amplititude.png")