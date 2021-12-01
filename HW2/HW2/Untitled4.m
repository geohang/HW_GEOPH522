clear


% 先限定三维图中的x,y轴坐标范围
X=1:700;
Y=1:500;
% 标准差
sigma = 500;


% Z = zeros( 51, 51 );
sig=[500:-0.5:0];
sig1=[250:1:400];
for col = 1 : 1 : 500
    Z(:, col ) =mynormpdf(1:700,350,sig(col)).';
end
for col = 1 : 1 : 500
    
if (50-col)>0
    
        Z(1:50-col+1, 500-col ) =mynormpdf([1:50-col+1],25,sig1(col)).';
        Z(650:700-col+1, 451+col ) =mynormpdf([650:700-col+1],675,sig1(col)).';
    end
    
end
z=Z./max(max(Z));
% z(:,end)=1

% 显示高斯函数的三维曲面
[x,y]=meshgrid(Y,X);
contourf(x,y,z)
colorbar



% [X,Y]=meshgrid(x,y);
%
% data=rand(500,700);
%
% image(x,y,data)
