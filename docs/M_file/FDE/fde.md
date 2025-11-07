---
sidebar_position: 4
---

# FDE 2D 案例

## 矩形波导

计算的分量：Ex Ey Ez Hx Hy Hz

对对角介电常数和对角磁导率

#### 矩形波导参数
```matlab
%% 全局参数
lx = 8;
ly = 6;
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
```

### 结构

![alt text](matlab/images/fde/structure.png)

### 有效折射率

```matlab
1.7556826063901874e+00
1.7510039383666316e+00
1.7078231749348238e+00
1.7051830771757968e+00
1.6702098530713665e+00
1.6436053137711506e+00
1.6346904714717176e+00
1.6191671395326253e+00
1.6165329314886545e+00
1.5972754141960221e+00
1.5326418063207210e+00
1.5221029861551907e+00
1.5207756370334444e+00
1.5205899123715299e+00
1.4925200271322971e+00
1.4627175584545149e+00
1.4577585555870196e+00
1.4090636582240483e+00
1.4058159208264256e+00
1.3996551948712739e+00
```

### 模场分布

![alt text](matlab/images/fde/Mode1.png) ![alt text](matlab/images/fde/Mode2.png) ![alt text](matlab/images/fde/Mode3.png) ![alt text](matlab/images/fde/Mode4.png) ![alt text](matlab/images/fde/Mode5.png) ![alt text](matlab/images/fde/Mode6.png) ![alt text](matlab/images/fde/Mode7.png) ![alt text](matlab/images/fde/Mode8.png) ![alt text](matlab/images/fde/Mode9.png) ![alt text](matlab/images/fde/Mode10.png) ![alt text](matlab/images/fde/Mode11.png) ![alt text](matlab/images/fde/Mode12.png) ![alt text](matlab/images/fde/Mode13.png) ![alt text](matlab/images/fde/Mode14.png) ![alt text](matlab/images/fde/Mode15.png) ![alt text](matlab/images/fde/Mode16.png) ![alt text](matlab/images/fde/Mode17.png) ![alt text](matlab/images/fde/Mode18.png) ![alt text](matlab/images/fde/Mode19.png) ![alt text](matlab/images/fde/Mode20.png) [alt text](fde.md) 