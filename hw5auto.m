%hw5auto.m 自动计算收敛阶

L = 4; % 总共的计算次数
Table = zeros(L,4); % 表格

N0 = 5; % 初始的剖分段数

%% 填表
for m = 1:L
    
    N = N0*2^(m - 1); % 加密网格
    
    Nx = N; Ny = N;
    
    if N == 80
        draw_solution = 1;
    else
        draw_solution = 0;
    end
    
    tic
    main % 计算结果
    toc
    
    Table(m,1) = L2; % 把L2误差填入
    Table(m,3) = H1; % 把H1误差填入
    
    if m > 1 % 确保有旧数据
        Table(m,2) = log(Table(m - 1,1)./Table(m,1))./log(2); % 计算L2 Error的阶
        Table(m,4) = log(Table(m - 1,3)./Table(m,3))./log(2); % 计算H1 Error的阶
    end
    
end