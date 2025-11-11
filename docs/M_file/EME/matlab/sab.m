function AB = sab(A, B)


SA = {A.R12, A.T21; A.T12, A.R21};
SB = {B.R12, B.T21; B.T12, B.R21};
D = SA{1, 2} * pinv(eye(size(SA{1,2}))-SB{1,1}*SA{2,2});
F = SB{2,1} * pinv(eye(size(SA{2,1}))-SA{2,2}*SB{1,1});

AB.R12 = SA{1,1} + D * SB{1,1} * SA{2,1};
AB.T21 = D * SB{1,2};
AB.T12 = F * SA{2,1};
AB.R21 = SB{2,2} + F * SA{2,2} * SB{1,2};


end
