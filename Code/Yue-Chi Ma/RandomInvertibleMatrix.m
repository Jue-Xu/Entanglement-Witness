function Matrix = RandomInvertibleMatrix( N )
%RANDOMINVERTIBLEMATRIX ����������ɵľ���Ĭ����һ�����棩
%   �˴���ʾ��ϸ˵��
Matrix = randn(N)+1i*randn(N);
% disp(['determinant = ',num2str(det(Matrix))]);
end

