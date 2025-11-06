

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


eps_z = I;
eps_z(s1:s2) = index^2;
eps_z(s1) = 0.5*( 1+ index^2);
eps_z(s2) = 0.5*( 1+ index^2);
eps_x = I;
eps_x(s1:s2-1) = index^2;

EPS_X = spdiags(eps_x,0,nx,nx);
EPS_Z = spdiags(eps_z , 0 ,nx ,nx);
MU_Y = speye(nx,nx);


%% TM

LXHY = spdiags([-I,I]/dx , [-1 0] , nx, nx);
LXHY(1,1)=0;
LXEZ = spdiags([-I,I]/dx , [0 1] , nx,nx);
LXEZ(end,:) = LXEZ(end-1,:);

A = EPS_X*LXEZ*inv(EPS_Z)*LXHY + k0*k0*EPS_X*MU_Y;

[hy,d]=eigs(A,nmodes,index^2 * k0^2);

neff =  sqrt(diag(d))/k0;

%% 过滤

idx =  (abs(neff-1))<1e-9;
hy(:,idx)=[];
neff(idx)=[];

%% 其它场分量

for midx = 1:length(neff)
    beta = k0* neff(midx);
    ex(:,midx) = - 1j*beta/k0*inv(EPS_X) *hy(:,midx);
end

ez =1/k0 * inv(EPS_Z) * LXHY*hy;
hy = [hy(1,:) ; 0.5*(hy(1:nx-1,:)+hy(2:nx,:))];
ex = [ex(1,:) ; 0.5*(ex(1:nx-1,:)+ex(2:nx,:))];
hy = 1j*hy/Z0;

%% 绘图

s = "images/TM/";
mkdir(s)

plot(xm , eps_x);
xlabel("X (m)")
ylabel("介电常数")
title("EPS X")
exportgraphics(gcf , s + "EPS_X.png")

plot(x , eps_z);
xlabel("X (m)")
ylabel("介电常数")
title("EPS Z")
exportgraphics(gcf , s + "EPS_Z.png")

%% 

 
plot(x , abs(hy));
xlabel("X (m)")
ylabel("Amplititude")
title("Hy Amplititude")
exportgraphics(gcf , s + "Hy_Amplititude.png")
 

plot(x , abs(ex));
xlabel("X (m)")
ylabel("Amplititude")
title("Ex Amplititude")
exportgraphics(gcf , s + "Ex_Amplititude.png")


plot(x , abs(ez));
xlabel("X (m)")
ylabel("Amplititude")
title("Ez Amplititude")
exportgraphics(gcf , s + "Ez_Amplititude.png")