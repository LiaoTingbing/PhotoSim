
%% 

close all ;
clear ;
clc ; 


%% 物理参数
EPS0 = 8.85e-12;
MU0 = pi * 4e-7;
C0 = 1 / sqrt(MU0*EPS0);
Z0 = sqrt(MU0/EPS0);

%% 全局参数
lx = 8;
ly = 0;
dx = 0.1;
dy = 0.1;
lambda = 1.55;
w = 3;
h = 2;
index = 1.8;
backIndex = 1.1;
nmodes = 20;
nmodes_max = 100;
toler = 1e-6;
search_index = 1.8;
%% 参数计算

k0 = 2 * pi / lambda;

cx = round(lx/dx);
cy = round(ly/dy);
nx = cx + 1;
ny = cy + 1;
nt = nx * ny;
I = ones(nt, 1);

%% 建模

x = (-lx / 2:dx:lx / 2)';
y = (-ly / 2:dy:ly / 2)';
xm = x + dx / 2;
ym = y + dy / 2;

s1 = [dsearchn(x, -w/2), dsearchn(x, w/2), dsearchn(y, -h/2), dsearchn(y, h/2)]';

eps_x = backIndex^2 * ones(nx, ny);
eps_y = backIndex^2 * ones(nx, ny);
eps_z = backIndex^2 * ones(nx, ny);

eps_x(s1(1):s1(2)-1, s1(3):s1(4)) = index^2;
eps_y(s1(1):s1(2), s1(3):s1(4)) = index^2;
eps_z(s1(1):s1(2), s1(3):s1(4)) = index^2;

% eps_x(s1(1):s1(2)-1, s1(3)) = 0.5 * (backIndex^2 + index^2);
% eps_x(s1(1):s1(2)-1, s1(4)) = 0.5 * (backIndex^2 + index^2);
%
% eps_y(s1(1), s1(3):s1(4)-1) = 0.5 * (backIndex^2 + index^2);
% eps_y(s1(2), s1(3):s1(4)-1) = 0.5 * (backIndex^2 + index^2);
%
% eps_z(s1(1):s1(2), s1(3)) = 0.5 * (backIndex^2 + index^2);
% eps_z(s1(1):s1(2), s1(4)) = 0.5 * (backIndex^2 + index^2);
% eps_z(s1(1), s1(3):s1(4)) = 0.5 * (backIndex^2 + index^2);
% eps_z(s1(2), s1(3):s1(4)) = 0.5 * (backIndex^2 + index^2);
%
% eps_z(s1(1), s1(3)) = 0.25 * (3 * backIndex^2 + index^2);
% eps_z(s1(1), s1(4)) = 0.25 * (3 * backIndex^2 + index^2);
% eps_z(s1(2), s1(3)) = 0.25 * (3 * backIndex^2 + index^2);
% eps_z(s1(2), s1(4)) = 0.25 * (3 * backIndex^2 + index^2);
%% 构建算子

UX = spdiags([-I, I]/dx, [0, 1], nt, nt);
for i = nx:nx:nt
    UX(i, :) = UX(i-1, :);
end

UY = sparse(nt, nt);

VX = spdiags([-I, I]/dx, [-1, 0], nt, nt);
for i = 1:nx:nt
    VX(i, :) = 0;
end

VY = sparse(nt, nt);

%% 构建矩阵

EPS_X = spdiags(eps_x(:), 0, nt, nt);
EPS_Y = spdiags(eps_y(:), 0, nt, nt);
EPS_Z = spdiags(eps_z(:), 0, nt, nt);

MU_X = speye(nt);
MU_Y = speye(nt);
MU_Z = speye(nt);

Pxx = -VX / k0 * inv(MU_Z) * UY;
Pxy = VX / k0 * inv(MU_Z) * UX + k0 * EPS_Y;
Pyx = -VY / k0 * inv(MU_Z) * UY - k0 * EPS_X;
Pyy = VY / k0 * inv(MU_Z) * UX;

P = [Pxx, Pxy; Pyx, Pyy];

Qxx = -UX / k0 * inv(EPS_Z) * VY;
Qxy = UX / k0 * inv(EPS_Z) * VX + k0 * MU_Y;
Qyx = -UY / k0 * inv(EPS_Z) * VY - k0 * MU_X;
Qyy = UY / k0 * inv(EPS_Z) * VX;

Q = [Qxx, Qxy; Qyx, Qyy];

[et_, d] = eigs(-Q*P, nmodes_max, search_index*search_index*k0*k0);
neff_ = sqrt(diag(d)) / k0;
%% 删除特征值

flag = ones(nmodes_max,1);
for i = 1:nmodes_max-1
    if neff_(i) - neff_(i+1) <toler
        flag(i) = 0;
        flag(i+1) = 0;
    end
end
neff = neff_(flag > 0);
et = et_(:, flag > 0);

neff = neff(1:nmodes);
et = et(:, 1:nmodes);
%% 计算场分量

ex = et(1:nt, :);
ey = et(nt+1:2*nt, :);

for midx = 1:nmodes
    beta = neff(midx) * k0;
    hx(:, midx) = -1i / beta * (Pxx * ex(:, midx) + Pxy * ey(:, midx));
    hy(:, midx) = -1i / beta * (Pyx * ex(:, midx) + Pyy * ey(:, midx));
end

ez = 1 / k0 * inv(EPS_Z) * (VX * hy - VY * hx);
hz = 1 / k0 * inv(MU_Z) * (UX * ey - UY * ex);

F.Ex = reshape(ex, nx, ny, nmodes);
F.Ey = reshape(ey, nx, ny, nmodes);
F.Ez = reshape(ez, nx, ny, nmodes);
F.Hx = reshape(hx, nx, ny, nmodes) * 1i / Z0;
F.Hy = reshape(hy, nx, ny, nmodes) * 1i / Z0;
F.Hz = reshape(hz, nx, ny, nmodes) * 1i / Z0;

%% 绘图
fp = "images/fde1d/";
mkdir(fp);

for midx = 1:nmodes

    subplot(231)
    plot(x,real (F.Ex(:, :, midx))');
    xlabel("X (m)")
    ylabel("Y (m)")
    subtitle("Ex real")
    colorbar

    subplot(232)
    plot(x,real (F.Ey(:, :, midx))');
    xlabel("X (m)")
    ylabel("Y (m)")
    subtitle("Ey real")
    colorbar

    subplot(233)
    plot(x, imag(F.Ez(:, :, midx))');
    xlabel("X (m)")
    ylabel("Y (m)")
    subtitle("Ez imag")
    colorbar

    subplot(234)
    plot(x,  real(F.Hx(:, :, midx))');
    xlabel("X (m)")
    ylabel("Y (m)")
    subtitle("Hx real")
    colorbar

    subplot(235)
    plot(x,  real(F.Hy(:, :, midx))');
    xlabel("X (m)")
    ylabel("Y (m)")
    subtitle("Hy real")
    colorbar

    subplot(236)
    plot(x, imag(F.Hz(:, :, midx))');
    xlabel("X (m)")
    ylabel("Y (m)")
    subtitle("Hz imag")
    colorbar

    sgtitle("Mode #"+num2str(midx)+"="+num2str(neff(midx)))

    set(gcf ,"OuterPosition",get(0, "ScreenSize"))

    exportgraphics(gcf , fp + "Mode" + num2str(midx) + ".png")
 

end

%% 结构
close

subplot(221)
plot(x, eps_x');
xlabel("X (m)")
ylabel("Y (m)")
subtitle("eps x")
colorbar

subplot(222)
plot(x, eps_y');
xlabel("X (m)")
ylabel("Y (m)")
subtitle("eps y")
colorbar

subplot(223)
plot(x, eps_z');
xlabel("X (m)")
ylabel("Y (m)")
subtitle("eps z")
colorbar

sgtitle("相对介电常数分布")

set(gcf , "outerposition" , get(0, "ScreenSize"))

exportgraphics(gcf , fp + "structure.png")

save(fp+"neff.txt", "neff", "-ascii", "-double");
