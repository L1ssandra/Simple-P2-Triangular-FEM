xa = 0;
xb = 1;
ya = 0;
yb = 1;
%Nx = 10;
%Ny = 10;
hx = (xb - xa)/Nx;
hy = (yb - ya)/Ny;
[x,y] = meshgrid(linspace(xa,xb,Nx + 1),linspace(ya,yb,Ny + 1));
x = reshape(x,(Nx + 1)*(Ny + 1),1);
y = reshape(y,(Nx + 1)*(Ny + 1),1);
DT = delaunayTriangulation(x,y);

draw_mesh = 0; % 绘制网格
draw_solution = 1; % 绘制数值解图像

if draw_mesh == 1
    triplot(DT)
    hold on
    vxlabels = arrayfun(@(n) {sprintf('P%d', n)}, (1:(Nx + 1)*(Ny + 1))');
    Hpl = text(x, y, vxlabels, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'BackgroundColor', 'none');
    ic = incenter(DT);
    numtri = size(DT,1);
    trilabels = arrayfun(@(x) {sprintf('T%d', x)}, (1:numtri)');
    Htl = text(ic(:,1), ic(:,2), trilabels, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', 'blue');
    axis([xa - 0.2*hx,xb + 0.2*hx,ya - 0.2*hy,yb + 0.2*hy])
    hold off
end

f = @(x,y) ((x - 1).*sin(x) - 2.*cos(x)).*(y - 1).*sin(y) + ((y - 1).*sin(y) - 2.*cos(y)).*(x - 1).*sin(x);
u = @(x,y) (x - 1).*(y - 1).*sin(x).*sin(y);
ux = @(x,y) (sin(x) + (x - 1).*cos(x)).*(y - 1).*sin(y);
uy = @(x,y) (sin(y) + (y - 1).*cos(y)).*(x - 1).*sin(x);

phi1 = @(x,y) x.*(2.*x - 1);
phi2 = @(x,y) y.*(2.*y - 1);
phi3 = @(x,y) (1 - x - y).*(2.*(1 - x - y) - 1);
phi4 = @(x,y) 4.*x.*y;
phi5 = @(x,y) 4.*y.*(1 - x - y);
phi6 = @(x,y) 4.*(1 - x - y).*x;

phi1x = @(x,y) 2.*x - 1 + 2.*x;
phi2x = @(x,y) 0.*x;
phi3x = @(x,y) -(2.*(1 - x - y) - 1) - 2.*(1 - x - y);
phi4x = @(x,y) 4.*y;
phi5x = @(x,y) -4.*y;
phi6x = @(x,y) 4.*( -1.*x + (1 - x - y) );

phi1y = @(x,y) 0.*x;
phi2y = @(x,y) 2.*y - 1 + 2.*y;
phi3y = @(x,y) -(2.*(1 - x - y) - 1) - 2.*(1 - x - y);
phi4y = @(x,y) 4.*x;
phi5y = @(x,y) 4.*( 1 - x - y - y );
phi6y = @(x,y) -4.*x;

Phix = { phi1x,phi2x,phi3x,phi4x,phi5x,phi6x };
Phiy = { phi1y,phi2y,phi3y,phi4y,phi5y,phi6y };

Iijxx = zeros(6,6);
Iijxy = zeros(6,6);
Iijyx = zeros(6,6);
Iijyy = zeros(6,6);

for i = 1:6
    for j = 1:6
        
        phiix = Phix{i};
        phiiy = Phiy{i};
        phijx = Phix{j};
        phijy = Phiy{j};
        
        phiijxx = @(x,y) phiix(x,y).*phijx(x,y);
        phiijxy = @(x,y) phiix(x,y).*phijy(x,y);
        phiijyx = @(x,y) phiiy(x,y).*phijx(x,y);
        phiijyy = @(x,y) phiiy(x,y).*phijy(x,y);
        
        Iijxx(i,j) = quad2d(phiijxx,0,1,0,@(x) 1 - x);
        Iijxy(i,j) = quad2d(phiijxy,0,1,0,@(x) 1 - x);
        Iijyx(i,j) = quad2d(phiijyx,0,1,0,@(x) 1 - x);
        Iijyy(i,j) = quad2d(phiijyy,0,1,0,@(x) 1 - x);
        
    end
end

% 构造局部刚度矩阵和右端项
[N,~] = size(DT.Points); % 节点总个数
[Nk,~] = size(DT.ConnectivityList); % 单元总个数
C2 = [DT.ConnectivityList,zeros(Nk,3)]; % 节点和中点的关系
P1 = DT.Points; % 节点
P2 = zeros(0,2); % 中点

for count = 1:Nk % 遍历所有单元，记录中点信息
    i = DT.ConnectivityList(count,1);
    j = DT.ConnectivityList(count,2);
    k = DT.ConnectivityList(count,3);
    
    xi = DT.Points(i,1); yi = DT.Points(i,2);
    xj = DT.Points(j,1); yj = DT.Points(j,2);
    xk = DT.Points(k,1); yk = DT.Points(k,2);
    
    Pij = 0.5*[xi + xj,yi + yj];
    Pjk = 0.5*[xj + xk,yj + yk];
    Pki = 0.5*[xk + xi,yk + yi];
    
    [N1,~] = size(P2);
    ij = 1;
    jk = 1;
    ki = 1;
    if N1 > 0
        for count2 = 1:N1
            if P2(count2,:) == Pij
                ij = 0;
                C2(count,4) = count2;
            end
            if P2(count2,:) == Pjk
                jk = 0;
                C2(count,5) = count2;
            end
            if P2(count2,:) == Pki
                ki = 0;
                C2(count,6) = count2;
            end
        end
    end
    if ij == 1
        P2 = [P2;Pij];
        [N1,~] = size(P2);
        C2(count,4) = N1;
    end
    if jk == 1
        P2 = [P2;Pjk];
        [N1,~] = size(P2);
        C2(count,5) = N1;
    end
    if ki == 1
        P2 = [P2;Pki];
        [N1,~] = size(P2);
        C2(count,6) = N1;
    end
end

C2(:,4:6) = C2(:,4:6) + N;

P2 = [DT.Points;P2];

[N2,~] = size(P2);

A = zeros(N2,N2);
F = zeros(N2,1);

for count = 1:Nk % 按单元循环
    
    i = C2(count,1);
    j = C2(count,2);
    k = C2(count,3);
    ij = C2(count,4);
    jk = C2(count,5);
    ki = C2(count,6);
    
    xi = P2(i,1); yi = P2(i,2);
    xj = P2(j,1); yj = P2(j,2);
    xk = P2(k,1); yk = P2(k,2);
    xij = P2(ij,1); yij = P2(ij,2);
    xjk = P2(jk,1); yjk = P2(jk,2);
    xki = P2(ki,1); yki = P2(ki,2);
    
    Delta = det([xi,yi,1;xj,yj,1;xk,yk,1])/2;
    
    Xx = (yj - yk)/(2*Delta); Xy = (xk - xj)/(2*Delta);
    Yx = (yk - yi)/(2*Delta); Yy = (xi - xk)/(2*Delta);
    
    A(i,i) = A(i,i) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(1,1) + (Xx*Yx + Xy*Yy)*(Iijxy(1,1) + Iijyx(1,1)) + (Yx^2 + Yy^2)*Iijyy(1,1) );
    A(i,j) = A(i,j) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(1,2) + (Xx*Yx + Xy*Yy)*(Iijxy(1,2) + Iijyx(1,2)) + (Yx^2 + Yy^2)*Iijyy(1,2) );
    A(i,k) = A(i,k) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(1,3) + (Xx*Yx + Xy*Yy)*(Iijxy(1,3) + Iijyx(1,3)) + (Yx^2 + Yy^2)*Iijyy(1,3) );
    A(i,ij) = A(i,ij) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(1,4) + (Xx*Yx + Xy*Yy)*(Iijxy(1,4) + Iijyx(1,4)) + (Yx^2 + Yy^2)*Iijyy(1,4) );
    A(i,jk) = A(i,jk) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(1,5) + (Xx*Yx + Xy*Yy)*(Iijxy(1,5) + Iijyx(1,5)) + (Yx^2 + Yy^2)*Iijyy(1,5) );
    A(i,ki) = A(i,ki) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(1,6) + (Xx*Yx + Xy*Yy)*(Iijxy(1,6) + Iijyx(1,6)) + (Yx^2 + Yy^2)*Iijyy(1,6) );
    
    A(j,i) = A(j,i) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(2,1) + (Xx*Yx + Xy*Yy)*(Iijxy(2,1) + Iijyx(2,1)) + (Yx^2 + Yy^2)*Iijyy(2,1) );
    A(j,j) = A(j,j) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(2,2) + (Xx*Yx + Xy*Yy)*(Iijxy(2,2) + Iijyx(2,2)) + (Yx^2 + Yy^2)*Iijyy(2,2) );
    A(j,k) = A(j,k) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(2,3) + (Xx*Yx + Xy*Yy)*(Iijxy(2,3) + Iijyx(2,3)) + (Yx^2 + Yy^2)*Iijyy(2,3) );
    A(j,ij) = A(j,ij) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(2,4) + (Xx*Yx + Xy*Yy)*(Iijxy(2,4) + Iijyx(2,4)) + (Yx^2 + Yy^2)*Iijyy(2,4) );
    A(j,jk) = A(j,jk) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(2,5) + (Xx*Yx + Xy*Yy)*(Iijxy(2,5) + Iijyx(2,5)) + (Yx^2 + Yy^2)*Iijyy(2,5) );
    A(j,ki) = A(j,ki) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(2,6) + (Xx*Yx + Xy*Yy)*(Iijxy(2,6) + Iijyx(2,6)) + (Yx^2 + Yy^2)*Iijyy(2,6) );
    
    A(k,i) = A(k,i) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(3,1) + (Xx*Yx + Xy*Yy)*(Iijxy(3,1) + Iijyx(3,1)) + (Yx^2 + Yy^2)*Iijyy(3,1) );
    A(k,j) = A(k,j) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(3,2) + (Xx*Yx + Xy*Yy)*(Iijxy(3,2) + Iijyx(3,2)) + (Yx^2 + Yy^2)*Iijyy(3,2) );
    A(k,k) = A(k,k) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(3,3) + (Xx*Yx + Xy*Yy)*(Iijxy(3,3) + Iijyx(3,3)) + (Yx^2 + Yy^2)*Iijyy(3,3) );
    A(k,ij) = A(k,ij) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(3,4) + (Xx*Yx + Xy*Yy)*(Iijxy(3,4) + Iijyx(3,4)) + (Yx^2 + Yy^2)*Iijyy(3,4) );
    A(k,jk) = A(k,jk) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(3,5) + (Xx*Yx + Xy*Yy)*(Iijxy(3,5) + Iijyx(3,5)) + (Yx^2 + Yy^2)*Iijyy(3,5) );
    A(k,ki) = A(k,ki) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(3,6) + (Xx*Yx + Xy*Yy)*(Iijxy(3,6) + Iijyx(3,6)) + (Yx^2 + Yy^2)*Iijyy(3,6) );
    
    A(ij,i) = A(ij,i) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(4,1) + (Xx*Yx + Xy*Yy)*(Iijxy(4,1) + Iijyx(4,1)) + (Yx^2 + Yy^2)*Iijyy(4,1) );
    A(ij,j) = A(ij,j) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(4,2) + (Xx*Yx + Xy*Yy)*(Iijxy(4,2) + Iijyx(4,2)) + (Yx^2 + Yy^2)*Iijyy(4,2) );
    A(ij,k) = A(ij,k) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(4,3) + (Xx*Yx + Xy*Yy)*(Iijxy(4,3) + Iijyx(4,3)) + (Yx^2 + Yy^2)*Iijyy(4,3) );
    A(ij,ij) = A(ij,ij) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(4,4) + (Xx*Yx + Xy*Yy)*(Iijxy(4,4) + Iijyx(4,4)) + (Yx^2 + Yy^2)*Iijyy(4,4) );
    A(ij,jk) = A(ij,jk) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(4,5) + (Xx*Yx + Xy*Yy)*(Iijxy(4,5) + Iijyx(4,5)) + (Yx^2 + Yy^2)*Iijyy(4,5) );
    A(ij,ki) = A(ij,ki) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(4,6) + (Xx*Yx + Xy*Yy)*(Iijxy(4,6) + Iijyx(4,6)) + (Yx^2 + Yy^2)*Iijyy(4,6) );
    
    A(jk,i) = A(jk,i) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(5,1) + (Xx*Yx + Xy*Yy)*(Iijxy(5,1) + Iijyx(5,1)) + (Yx^2 + Yy^2)*Iijyy(5,1) );
    A(jk,j) = A(jk,j) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(5,2) + (Xx*Yx + Xy*Yy)*(Iijxy(5,2) + Iijyx(5,2)) + (Yx^2 + Yy^2)*Iijyy(5,2) );
    A(jk,k) = A(jk,k) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(5,3) + (Xx*Yx + Xy*Yy)*(Iijxy(5,3) + Iijyx(5,3)) + (Yx^2 + Yy^2)*Iijyy(5,3) );
    A(jk,ij) = A(jk,ij) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(5,4) + (Xx*Yx + Xy*Yy)*(Iijxy(5,4) + Iijyx(5,4)) + (Yx^2 + Yy^2)*Iijyy(5,4) );
    A(jk,jk) = A(jk,jk) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(5,5) + (Xx*Yx + Xy*Yy)*(Iijxy(5,5) + Iijyx(5,5)) + (Yx^2 + Yy^2)*Iijyy(5,5) );
    A(jk,ki) = A(jk,ki) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(5,6) + (Xx*Yx + Xy*Yy)*(Iijxy(5,6) + Iijyx(5,6)) + (Yx^2 + Yy^2)*Iijyy(5,6) );
    
    A(ki,i) = A(ki,i) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(6,1) + (Xx*Yx + Xy*Yy)*(Iijxy(6,1) + Iijyx(6,1)) + (Yx^2 + Yy^2)*Iijyy(6,1) );
    A(ki,j) = A(ki,j) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(6,2) + (Xx*Yx + Xy*Yy)*(Iijxy(6,2) + Iijyx(6,2)) + (Yx^2 + Yy^2)*Iijyy(6,2) );
    A(ki,k) = A(ki,k) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(6,3) + (Xx*Yx + Xy*Yy)*(Iijxy(6,3) + Iijyx(6,3)) + (Yx^2 + Yy^2)*Iijyy(6,3) );
    A(ki,ij) = A(ki,ij) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(6,4) + (Xx*Yx + Xy*Yy)*(Iijxy(6,4) + Iijyx(6,4)) + (Yx^2 + Yy^2)*Iijyy(6,4) );
    A(ki,jk) = A(ki,jk) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(6,5) + (Xx*Yx + Xy*Yy)*(Iijxy(6,5) + Iijyx(6,5)) + (Yx^2 + Yy^2)*Iijyy(6,5) );
    A(ki,ki) = A(ki,ki) + 2*Delta*( (Xx^2 + Xy^2)*Iijxx(6,6) + (Xx*Yx + Xy*Yy)*(Iijxy(6,6) + Iijyx(6,6)) + (Yx^2 + Yy^2)*Iijyy(6,6) );
    
    f1 = @(x,y) f( (xi - xk).*x + (xj - xk).*y + xk , (yi - yk).*x + (yj - yk).*y + yk );
    fphi1 = @(x,y) f1(x,y).*phi1(x,y);
    fphi2 = @(x,y) f1(x,y).*phi2(x,y);
    fphi3 = @(x,y) f1(x,y).*phi3(x,y);
    fphi4 = @(x,y) f1(x,y).*phi4(x,y);
    fphi5 = @(x,y) f1(x,y).*phi5(x,y);
    fphi6 = @(x,y) f1(x,y).*phi6(x,y);
    
    F(i) = F(i) + quad2d(fphi1,0,1,0,@(x) 1 - x)*(2*Delta);
    F(j) = F(j) + quad2d(fphi2,0,1,0,@(x) 1 - x)*(2*Delta);
    F(k) = F(k) + quad2d(fphi3,0,1,0,@(x) 1 - x)*(2*Delta);
    F(ij) = F(ij) + quad2d(fphi4,0,1,0,@(x) 1 - x)*(2*Delta);
    F(jk) = F(jk) + quad2d(fphi5,0,1,0,@(x) 1 - x)*(2*Delta);
    F(ki) = F(ki) + quad2d(fphi6,0,1,0,@(x) 1 - x)*(2*Delta);

end

% 赋边界条件，按节点循环
[N2,~] = size(P2);
for count = 1:N2
    xi = P2(count,1); yi = P2(count,2);
    if (xi == 0) || (xi == 1) || (yi == 0) || (yi == 1)
        A(count,:) = 0;
        A(:,count) = 0;
        A(count,count) = 1;
        F(count) = 0;
    end
end

U = A\F;
U1 = U(1:N);

for i = 1:N2
    for j = 1:N2
        if abs(A(i,j)) < 1e-7
            A(i,j) = 0;
        end
    end
end

if draw_solution == 1
    trimesh(DT.ConnectivityList,x,y,U1);
    colormap(cool)
end


% 计算误差的L2和H1模，按单元循环
L2 = 0;
H1 = 0;
for count = 1:Nk
    
    i = C2(count,1);
    j = C2(count,2);
    k = C2(count,3);
    ij = C2(count,4);
    jk = C2(count,5);
    ki = C2(count,6);
    
    xi = P2(i,1); yi = P2(i,2);
    xj = P2(j,1); yj = P2(j,2);
    xk = P2(k,1); yk = P2(k,2);
    xij = P2(ij,1); yij = P2(ij,2);
    xjk = P2(jk,1); yjk = P2(jk,2);
    xki = P2(ki,1); yki = P2(ki,2);
    
    % 逆变换
    Lx = @(x,y) (xi - xk).*x + (xj - xk).*y + xk;
    Ly = @(x,y) (yi - yk).*x + (yj - yk).*y + yk;
    
    Delta = det([xi,yi,1;xj,yj,1;xk,yk,1])/2;
    
    Xx = (yj - yk)/(2*Delta); Xy = (xk - xj)/(2*Delta);
    Yx = (yk - yi)/(2*Delta); Yy = (xi - xk)/(2*Delta);
    
    uh = @(x,y) U(i).*phi1(x,y) + U(j).*phi2(x,y) + U(k).*phi3(x,y) + U(ij).*phi4(x,y) + U(jk).*phi5(x,y) + U(ki).*phi6(x,y);
    e2 = @(x,y) ( uh(x,y) - u(Lx(x,y),Ly(x,y)) ).^2;
    
    uhx = @(x,y) (U(i).*phi1x(x,y) + U(j).*phi2x(x,y) + U(k).*phi3x(x,y) + U(ij).*phi4x(x,y) + U(jk).*phi5x(x,y) + U(ki).*phi6x(x,y)).*Xx...
        + (U(i).*phi1y(x,y) + U(j).*phi2y(x,y) + U(k).*phi3y(x,y) + U(ij).*phi4y(x,y) + U(jk).*phi5y(x,y) + U(ki).*phi6y(x,y)).*Yx;
    ex2 = @(x,y) ( uhx(x,y) - ux(Lx(x,y),Ly(x,y)) ).^2;
    
    uhy = @(x,y) (U(i).*phi1x(x,y) + U(j).*phi2x(x,y) + U(k).*phi3x(x,y) + U(ij).*phi4x(x,y) + U(jk).*phi5x(x,y) + U(ki).*phi6x(x,y)).*Xy...
        + (U(i).*phi1y(x,y) + U(j).*phi2y(x,y) + U(k).*phi3y(x,y) + U(ij).*phi4y(x,y) + U(jk).*phi5y(x,y) + U(ki).*phi6y(x,y)).*Yy;
    ey2 = @(x,y) ( uhy(x,y) - uy(Lx(x,y),Ly(x,y)) ).^2;
    
    L2 = L2 + quad2d(e2,0,1,0,@(x) 1 - x)*(2*Delta);
    
    H1 = H1 + (quad2d(ex2,0,1,0,@(x) 1 - x) + quad2d(ey2,0,1,0,@(x) 1 - x))*(2*Delta);
end
L2 = sqrt(L2);
H1 = sqrt(L2.^2 + H1);