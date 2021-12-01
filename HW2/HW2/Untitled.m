

%% Title: A_Star_4_Path
% Authors: ruogu7
% Email: 380545156@qq.com
% Start time: 8:30 am,Dec 11th,2019
% Latest update: 12th,Dec,2019
%% Calculation steps
% 1. Read .xls file and visualize points data
% 2. Path routing
% 3. Smoothing the path
% 4. Segment the smoothed path into minimal steps, and calculate the slope
% rate of each step.
% 5. Visualize the result generated by steps above
% 6. Make a vedio or GIF file to record the dynamic process
%% 1. 初始化参数
close all;clc;clear;
% 地图大小
m = 30;
n = 30;
% 起始位置
Spoint = [3 3];	  %起始点坐标
Epoint = [29 22];	%目标点坐标

%% 2. 构建地图
% -inf表示不可到达的障碍物点
%%构建地图
for i = 1:m+2
    if i == 1
        for j = 1:n+2
            Matrix(i,j) = -inf;
        end
    elseif i == m+2
        for j = 1:n+2
            Matrix(i,j) = -inf;
        end
    else
        for j = 1:n+2
            if ((j == 1)|(j == n+2))
                Matrix(i,j) = -inf;
            else
                Matrix(i,j) = inf;
            end
        end
    end
end
% %%障碍
% for j=2:10
%     Matrix(5,j)=-inf;
%     for j=2:15
%         Matrix(24,j)=-inf;
%         for j=9:24
%             %for j=6:24
%             Matrix(10,j)=-inf;
%             for j=20:31
%                 Matrix(15,j)=-inf;
%                 for j=5:20
%                     Matrix(20,j)=-inf;
%                     for j=18:27
%                         Matrix(28,j)=-inf;
%                         for i=2:6
%                             Matrix(i,18)=-inf;
%                             for i=17:20
%                                 Matrix(i,5)=-inf;
%                                 for i=23:25
%                                     Matrix(i,20)=-inf;
%                                     for i=13:17
%                                         Matrix(i,13)=-inf;
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% 
% Matrix=zeros(32,32);
% Matrix(10,10)=-inf;


% 显示地图
%subplot(2,2,1);
h1 = plot(Spoint(1),Spoint(2),'gO');
hold on;
h2 = plot(Epoint(1),Epoint(2),'rO');
title('Route planing with A* algorithms');
% axes('pos',[.1 .6 .5 .3])
% [bottomleftcornerXposition bottomleftcornerYposition width height]
% imshow('coins.jpg');
annotation('arrow', [.3 .5], [.6 .5]);

%% 3. A*算法搜索路径
%%寻路
Matrix(Spoint(1),Spoint(2))=0;
Matrix(Epoint(1),Epoint(2))=inf;
G=Matrix;
F=Matrix;
openlist=Matrix;
closelist=Matrix;
parentx=Matrix;
parenty=Matrix;
openlist(Spoint(1),Spoint(2)) =0;
%closelist(Epoint(1),Epoint(2))=inf;

for i = 1:n+2
    for j = 1:m+2
        k = Matrix(i,j);
        if(k == -inf)
            %subplot(2,2,1);
            h3 = plot(i,j,'k.');
            %         elseif(k == inf)  % show green feasible point
            %             %subplot(2,2,1);
            %             plot(i,j,'gh');
            %         else
            %             %subplot(2,2,1);
            %             plot(i,j,'gh');
        end
        hold on
    end
end
axis([0 m+3 0 n+3]);
%subplot(2,2,1);
plot(Epoint(1),Epoint(2),'b+');
%subplot(2,2,1);
plot(Spoint(1),Spoint(2),'b+');
while(1)
    num=inf;
    for p=1:m+2
        for q=1:n+2
            if(openlist(p,q)==0&&closelist(p,q)~=1)
                Outpoint=[p,q];
                if(F(p,q)>=0&&num>F(p,q))
                    num=F(p,q);
                    Nextpoint=[p,q];
                end
            end
        end
    end
    closelist(Nextpoint(1),Nextpoint(2))=1;
    for i = 1:3
        for j = 1:3
            k = G(Nextpoint(1)-2+i,Nextpoint(2)-2+j);
            if(i==2&&j==2||closelist(Nextpoint(1)-2+i,Nextpoint(2)-2+j)==1)
                continue;
            elseif (k == -inf)
                G(Nextpoint(1)-2+i,Nextpoint(2)-2+j) = G(Nextpoint(1)-2+i,Nextpoint(2)-2+j);
                closelist(Nextpoint(1)-2+i,Nextpoint(2)-2+j)=1;
            elseif (k == inf)
                distance=((i-2)^2+(j-2)^2)^0.5;
                G(Nextpoint(1)-2+i,Nextpoint(2)-2+j)=G(Nextpoint(1),Nextpoint(2))+distance;
                openlist(Nextpoint(1)-2+i,Nextpoint(2)-2+j)=0;
                % H=((Nextpoint(1)-2+i-Epoint(1))^2+(Nextpoint(2)-2+j-Epoint(2))^2)^0.5;%欧几里德距离启发函数
                H_diagonal=min(abs(Nextpoint(1)-2+i-Epoint(1)),abs(Nextpoint(2)-2+j-Epoint(2)));%比较复杂的对角线启发函数
                H_straight=abs(Nextpoint(1)-2+i-Epoint(1))+abs(Nextpoint(2)-2+j-Epoint(2));
                H=sqrt(2)*H_diagonal+(H_straight-2*H_diagonal);
                % H=max(abs(Nextpoint(1)-2+i-Epoint(1)),abs(Nextpoint(2)-2+j-Epoint(2)));%比较简单的对角线函数
                F(Nextpoint(1)-2+i,Nextpoint(2)-2+j)=G(Nextpoint(1)-2+i,Nextpoint(2)-2+j)+H;
                parentx(Nextpoint(1)-2+i,Nextpoint(2)-2+j)=Nextpoint(1);
                parenty(Nextpoint(1)-2+i,Nextpoint(2)-2+j)=Nextpoint(2);
            else distance=((i-2)^2+(j-2)^2)^0.5;
                if(k>(distance+G(Nextpoint(1),Nextpoint(2))))
                    k=distance+G(Nextpoint(1),Nextpoint(2));
                    % H=((Nextpoint(1)-2+i-Epoint(1))^2+(Nextpoint(2)-2+j-Epoint(2))^2)^0.5;  %欧几里德距离启发函数
                    H_diagonal=min(abs(Nextpoint(1)-2+i-Epoint(1)),abs(Nextpoint(2)-2+j-Epoint(2)));%比较复杂的对角线启发函数
                    H_straight=abs(Nextpoint(1)-2+i-Epoint(1))+abs(Nextpoint(2)-2+j-Epoint(2));
                    H=sqrt(2)*10*H_diagonal+10*(H_straight-2*H_diagonal);
                    % H=max(abs(Nextpoint(1)-2+i-Epoint(1)),abs(Nextpoint(2)-2+j-Epoint(2)));%比较简单的对角线函数
                    F(Nextpoint(1)-2+i,Nextpoint(2)-2+j)=k+H;
                    parentx(Nextpoint(1)-2+i,Nextpoint(2)-2+j)=Nextpoint(1);
                    parenty(Nextpoint(1)-2+i,Nextpoint(2)-2+j)=Nextpoint(2);
                end
            end
            if(((Nextpoint(1)-2+i)==Epoint(1)&&(Nextpoint(2)-2+j)==Epoint(2))|num==inf)
                parentx(Epoint(1),Epoint(2))=Nextpoint(1);
                parenty(Epoint(1),Epoint(2))=Nextpoint(2);
                break;
            end
        end
        if(((Nextpoint(1)-2+i)==Epoint(1)&&(Nextpoint(2)-2+j)==Epoint(2))|num==inf)
            parentx(Epoint(1),Epoint(2))=Nextpoint(1);
            parenty(Epoint(1),Epoint(2))=Nextpoint(2);
            break;
        end
    end
    if(((Nextpoint(1)-2+i)==Epoint(1)&&(Nextpoint(2)-2+j)==Epoint(2))|num==inf)
        parentx(Epoint(1),Epoint(2))=Nextpoint(1);
        parenty(Epoint(1),Epoint(2))=Nextpoint(2);
        break;
    end
end
P=[];
s=1;
while(1)
    if(num==inf)
        break;
    end
    %subplot(2,2,1);
    h4 = plot(Epoint(1),Epoint(2),'b+');
    P(s,:)=Epoint;
    s=s+1;
    % pause(1);
    xx=Epoint(1);
    Epoint(1)=parentx(Epoint(1),Epoint(2));
    Epoint(2)=parenty(xx,Epoint(2));
    if(parentx(Epoint(1),Epoint(2))==Spoint(1)&&parenty(Epoint(1),Epoint(2))==Spoint(2))
        %subplot(2,2,1);
        plot(Epoint(1),Epoint(2),'b+');
        P(s,:)=Epoint;
        break;
    end
end
P(s+1,:)=Spoint;
legend([h1,h2,h3,h4],'Start Point','End Point','Obstacle','Routing');

count=0;
for i=2:12
    for j=2:12
        if(G(i,j)~=inf&&G(i,j)~=-inf)
            count=count+1;
        end
    end
end
% count

%% 4. 路径优化
%将得到的折现曲线拟合成光滑的曲线
P=P';
a=[];
b=[];
a=P(1,:);
figure
%subplot(2,2,3);
plot(a,b);
axis([0,n+3,0,n+3]);
title('Route planing');
values = spcrv([[a(1) a a(end)];[b(1) b b(end)]],3);
%% 5. 结果可视化
figure
%subplot(2,2,4);
plot(values(1,:),values(2,:),'r');
axis([0,m+3,0,m+3]);
title('Route planing after smoothing');