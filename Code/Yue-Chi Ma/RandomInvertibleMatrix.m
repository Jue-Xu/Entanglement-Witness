function Matrix = RandomInvertibleMatrix( N )
%RANDOMINVERTIBLEMATRIX 返回随机生成的矩阵（默认其一定可逆）
%   此处显示详细说明
Matrix = randn(N)+1i*randn(N);
% disp(['determinant = ',num2str(det(Matrix))]);
end

