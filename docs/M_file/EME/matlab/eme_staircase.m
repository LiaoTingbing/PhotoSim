%%
clc ;
close all;
clear;

%%
lx = 10;
ly = 8;
lz = 4;
ncells = 3;
nmodes = 6;
port_mode_index = [1, 1];
dy = 0.02;
dz = 0.02;
lambda = 1.55;
%%
dx = lx / 3;
x = linspace(-lx/2, lx/2, ncells+1)';
cell_pos = 0.5 * (x(1:ncells) + x(2:ncells+1));
port_pos = [-lx / 2, lx / 2]';
y = (-ly / 2:dy:ly / 2)';
z = (-lz / 2:dz:lz / 2)';
ny = length(y);
nz = length(z);
%%

for cidx = 1:ncells
    cell{cidx}.mode = load("../lumerical/cell_"+num2str(cidx)+"_mode.mat");
    cell{cidx}.neff = load("../lumerical/cell_"+num2str(cidx)+"_neff.mat");
end

for pidx = 1:length(port_mode_index)
    port{pidx}.mode = load("../lumerical/port_"+num2str(pidx)+"_mode.mat");
    port{pidx}.neff = load("../lumerical/port_"+num2str(pidx)+"_neff.mat");
end
%% 加载矩阵

for cidx = 1:ncells
    cell_neff(:, cidx) = squeeze(cell{cidx}.neff.neff.neff);
    for midx = 1:nmodes
        s = sprintf("[cell{%d}.mode.mode.E%d cell{%d}.mode.mode.H%d]", cidx, midx, cidx, midx);
        cell_modes(:, :, midx, cidx) = eval(s);

    end
end

% port
for pidx = 1:length(port_mode_index)
    port_neff(:, pidx) = squeeze(port{pidx}.neff.neff.neff);
    for midx = port_mode_index(pidx)
        s = sprintf("[port{%d}.mode.mode.E%d port{%d}.mode.mode.H%d]", pidx, midx, pidx, midx);
        port_modes(:, :, midx, pidx) = eval(s);

    end
end
%% 归一化 已经归一化
ds = dz * dy * 1e-12;
for cidx = 1:ncells
    for midx = 1:nmodes
        e = cell_modes(:, 1:3, midx, cidx);
        h = cell_modes(:, 4:6, midx, cidx);
        S = real(cross(e, h));
        Px = trapz(S(:, 1)) * ds;
        zx(midx) = Px;
        cell_modes(:, :, midx, cidx) = cell_modes(:, :, midx, cidx) / sqrt(Px);
    end
end
%
for pidx = 1:length(port_mode_index)
    for midx = port_mode_index(pidx)
        e = port_modes(:, 1:3, midx, pidx);
        h = port_modes(:, 4:6, midx, pidx);
        S = real(cross(e, h));
        Px = trapz(S(:, 1)) * ds;
        port_modes(:, :, midx, pidx) = port_modes(:, :, midx, pidx) / sqrt(Px);
    end
end
%% port1 积分
% port1 cell
e_p1 = port_modes(:, 1:3, port_mode_index(1), 1);
h_p1 = port_modes(:, 4:6, port_mode_index(1), 1);
for midx = 1:nmodes
    e_c = cell_modes(:, 1:3, midx, 1);
    h_c = cell_modes(:, 4:6, midx, 1);
    S = cross(e_p1, (h_c));
    port_1.overlap_farward(1, midx) = trapz(S(:, 1)) * ds;
    S = cross(e_c, (h_p1));
    port_1.overlap_backward(midx, 1) = trapz(S(:, 1)) * ds;
end
%% 第一个cell

for cidx = 1
    for i = 1:nmodes
        for j = 1:nmodes
            ei = cell_modes(:, 1:3, i, cidx);
            hj = cell_modes(:, 4:6, j, cidx+1);
            S = cross(ei, hj);
            cell_{cidx}.overlap_farward(i, j) = trapz(S(:, 1)) * ds;
        end
    end
end
%% 中间的Cell

for cidx = 1 + 1:ncells - 1
    for i = 1:nmodes
        for j = 1:nmodes
            ei = cell_modes(:, 1:3, i, cidx);
            hj = cell_modes(:, 4:6, j, cidx+1);
            S = cross(ei, hj);
            cell_{cidx}.overlap_farward(i, j) = trapz(S(:, 1)) * ds;

            ei = cell_modes(:, 1:3, i, cidx);
            hj = cell_modes(:, 4:6, j, cidx-1);
            S = cross(ei, hj);
            cell_{cidx}.overlap_backward(i, j) = trapz(S(:, 1)) * ds;
        end
    end
end

%% 最后一个cell

for cidx = ncells
    for i = 1:nmodes
        for j = 1:nmodes
            ei = cell_modes(:, 1:3, i, cidx);
            hj = cell_modes(:, 4:6, j, cidx-1);
            S = cross(ei, hj);
            cell_{cidx}.overlap_backward(i, j) = trapz(S(:, 1)) * ds;
        end
    end
end

%% port 2 积分
% cell_end port2
e_p2 = port_modes(:, 1:3, port_mode_index(2), 2);
h_p2 = port_modes(:, 4:6, port_mode_index(2), 2);
for midx = 1:nmodes
    e_c = cell_modes(:, 1:3, midx, ncells);
    h_c = cell_modes(:, 4:6, midx, ncells);
    S = cross(e_p2, (h_c));
    port_2.overlap_backward(1, midx) = trapz(S(:, 1)) * ds;
    S = cross(e_c, (h_p2));
    port_2.overlap_farward(midx, 1) = trapz(S(:, 1)) * ds;
end
%% 界面散射矩阵

for cidx = 1:ncells - 1

    O12 = cell_{cidx}.overlap_farward;
    O21 = cell_{cidx+1}.overlap_backward;
    T12 = (0.5 * (O12 + O21'));
    R12 = 0.5 * (O21' - O12) / T12;
    
    % T =O12 \( eye(nmodes) - R12)

    % T12 = 2 * pinv((O12 + O21'));
    % R12 = 0.5 * (O21' - O12) * T12;

    % M = O21'/O12 ; 
    % R12 = (eye(nmodes)+M)\(M-eye(nmodes));

        T21 = 0.5 *  ((O21 + O12'));
    R21 = 0.5 * (O12' - O21) / T21;

    % T21 = 2 * pinv((O21 + O12'));
    % R21 = 0.5 * (O12' - O21) * T21;

    cell_{cidx}.S_farward.S = [R12, T21; T12, R21];
    cell_{cidx}.S_farward.R12 = R12;
    cell_{cidx}.S_farward.T21 = T21;
    cell_{cidx}.S_farward.T12 = T12;
    cell_{cidx}.S_farward.R21 = R21;

end

%% 传播矩阵
k0 = 2 * pi / lambda;
cell_span = dx;
for cidx = 1:ncells

    neff_cidx = squeeze(cell{cidx}.neff.neff.neff);
    T12 = diag(exp(1i*k0*neff_cidx*cell_span));
    % 1->2
    R12 = zeros(nmodes);

    % 2-1>
    T21 = diag(exp(-1i*k0*neff_cidx*cell_span));
    R21 = zeros(nmodes);

    cell_{cidx}.prop.S = [R12, T21; T12, R21];
    cell_{cidx}.prop.R12 = R12;
    cell_{cidx}.prop.T21 = T21;
    cell_{cidx}.prop.T12 = T12;
    cell_{cidx}.prop.R21 = R21;
end
%% port 1 参数矩阵
O12 = zeros(nmodes);
O12(port_mode_index(1), :) = port_1.overlap_farward;

O21 = zeros(nmodes);
O21(:, port_mode_index(1)) = port_1.overlap_backward;

T12 = 0.5*((O12 + O21'))';
R12 = 0.5 * (O21' - O12)*pinv(T12);

 

1- O12(1,:)*(T12(:,1))

T21 = 2*pinv((O21 + O12'));
R21 = 0.5 * (O12' - O21)*( T21);

port_1.S = [R12, T21; T12, R21];
port_1.R12 = R12;
port_1.T21 = T21;
port_1.T12 = T12;
port_1.R21 = R21;
%% port2

O12 = zeros(nmodes);
O12(:, port_mode_index(2)) = port_2.overlap_farward;
O21 = zeros(nmodes);
O21(port_mode_index(2), :) = port_2.overlap_backward;

T12 = 2*pinv((O12 + O21'));
R12 = 0.5 * (O21' - O12)* (T12);

T21 =  2*pinv((O21 + O12'));
R21 = 0.5 * (O12' - O21) *( T21);


port_2.S = [R12, T21; T12, R21];
port_2.R12 = R12;
port_2.T21 = T21;
port_2.T12 = T12;
port_2.R21 = R21;
%% 内部S

% is.R12 = zeros(nmodes);
% is.T21 = eye(nmodes);
% is.T12 = eye(nmodes);
% is.R21 = zeros(nmodes);

is = cell_{1}.prop.S;

for cidx = 1:ncells - 1
    is = star_product(is, cell_{cidx}.S_farward.S);
    is = star_product(is, cell_{cidx+1}.prop.S);


end


%% 

is = star_product( port_1.S , is);
is = star_product(is , port_2.S);


%% 

us = abs(is)
