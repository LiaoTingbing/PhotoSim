function S_AB = star_product(S_A, S_B)
    % 计算两个散射矩阵S_A和S_B的Redheffer星积S_AB = S_A ⊗ S_B
    % 输入：S_A, S_B - 两个待级联的散射矩阵（尺寸需满足分块运算要求）
    % 输出：S_AB - 星积后的散射矩阵
    
    % 提取分块子矩阵
    S11_A = S_A(1:end/2, 1:end/2);  % S_A的11分块
    S12_A = S_A(1:end/2, end/2+1:end);  % S_A的12分块
    S21_A = S_A(end/2+1:end, 1:end/2);  % S_A的21分块
    S22_A = S_A(end/2+1:end, end/2+1:end);  % S_A的22分块
    
    S11_B = S_B(1:end/2, 1:end/2);  % S_B的11分块
    S12_B = S_B(1:end/2, end/2+1:end);  % S_B的12分块
    S21_B = S_B(end/2+1:end, 1:end/2);  % S_B的21分块
    S22_B = S_B(end/2+1:end, end/2+1:end);  % S_B的22分块
    
    % 计算重复项D和F
    D = S12_A * inv(eye(size(S11_B)) - S11_B * S22_A);
    F = S21_B * inv(eye(size(S22_A)) - S22_A * S11_B);
    
    % 计算星积后的各分块
    S11_AB = S11_A + D * S11_B * S21_A;
    S12_AB = D * S12_B;
    S21_AB = F * S21_A;
    S22_AB = S22_B + F * S22_A * S12_B;
    
    % 组合分块得到总矩阵
    S_AB = [S11_AB, S12_AB;
            S21_AB, S22_AB];
end