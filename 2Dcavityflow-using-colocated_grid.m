% 初期化
clear all;

% グローバル変数宣言
global x y dt nx ny ite_max dLx dLy Lx Ly rho nu

% パラメーター
Lx = 1;% 計算領域
Ly = 1;% 計算領域
nx = 40;% 要素数
ny = 40;% 要素数
dLx = Lx / nx;% 要素サイズ
dLy = Ly / ny;% 要素サイズ
nu = 0.001;% 動粘性係数
rho = 1;% 密度
t_max = 50;% SIM計算時間
dt = 0.001;% タイムステップ
ite_max = t_max / dt;% 反復数

% 座標の生成
x = 0 : Lx / nx : Lx;
y = Ly : - Ly / ny : 0;
[X, Y] = meshgrid(x, y);

% 配列確保
u = zeros(nx, ny);% セル中心での速度
v = zeros(nx, ny);% セル中心での速度
p = zeros(nx, ny);% セル中心での圧力
dp = zeros(nx, ny);% セル中心での補正圧力
diffx = zeros(nx, ny);
diffy = zeros(nx, ny);
advx = zeros(nx, ny);
advy = zeros(nx, ny);
usb = zeros(nx - 1, ny - 2);% セル界面での中間速度
vsb = zeros(nx - 2, ny - 1);% セル界面での中間速度

% 境界条件
bc.ue = 0;% 東側
bc.us = 0;% 南側
bc.uw = 0;% 西側
bc.un = 1;% 北側
bc.ve = 0;% 東側
bc.vs = 0;% 南側
bc.vw = 0;% 西側
bc.vn = 0;% 北側

% 境界条件の設定
[u, v] = setBC(u, v, bc);

% 時間発展
for ite = 2 : ite_max
    
    % 拡散項の計算
    [diffx, diffy] = diffusion(u, v, diffx, diffy);
    
    % 移流項の計算
    [advx, advy] = advection(u, v, advx, advy);
    
    % 移流拡散による時間進行を計算
    u = u + (nu * diffx + advx) * dt;
    v = v + (nu * diffy + advy) * dt;
    
    % セル界面での中間速度を計算
    [usb, vsb] = midvel(usb, vsb, u, v, p);
    
    % 補正圧力の計算
    [dp, ds] = poisson(usb, vsb, dp);
    
    % 次ステップ圧力の計算
    p = p + dp;
    
    % 次ステップ流速の計算
    [u, v] = velupdate(u, v, p);
    
    % 境界条件の設定
    [u, v] = setBC(u, v, bc);
    
    % コマンドウィンドウへの出力
    txt = ['ite = ', num2str(ite), ' / ', num2str(ite_max)];
    disp(txt);
    
    % リアルタイム可視化
%     vis_contour('u', ite, u, -0.8, 0.8, 1);
%     vis_contour('v', ite, v, -0.5, 0.5, 2);
%     vis_vector('vec', ite, u, v, 3)
    %vis_contour('p', ite, p, -3, 3, 4);
    %vis_contour('ds', ite, ds, -0.01, 0.01, 5);
    
end

%% 以下関数
function[u, v] = setBC(u, v, bc)

% グローバル変数呼び出し
global nx ny

u(nx, :) =  bc.ue;
u(:, 1) =  bc.us;
u(1, :) =  bc.uw;
u(:, ny) =  bc.un;

v(nx, :) =  bc.ve;
v(:, 1) =  bc.vs;
v(1, :) =  bc.vw;
v(:, ny) =  bc.vn;

end

function[diffx, diffy] = diffusion(u, v, diffx, diffy)

% グローバル変数呼び出し
global nx ny dLx dLy

for i = 2 : nx - 1
    for j = 2 : ny - 1
        
        diffx(i, j) =  (u(i + 1, j) - 2 * u(i, j) + u(i - 1, j))/(dLx)^2 + (u(i , j + 1) - 2 * u(i, j) + u(i, j - 1))/(dLy)^2 ;
        diffy(i, j) =  (v(i + 1, j) - 2 * v(i, j) + v(i - 1, j))/(dLx)^2 + (v(i , j + 1) - 2 * v(i, j) + v(i, j - 1))/(dLy)^2 ;
        
    end
end

end

function[advx, advy] = advection(u, v, advx, advy)

% グローバル変数呼び出し
global nx ny dLx dLy

for i = 2 : nx - 1
    for j = 2 : ny - 1
        
        advx(i, j) =  -(u(i, j) * (u(i + 1, j) - u(i - 1, j)) / (2 * dLx) + abs(u(i, j)) * (2 * u(i, j) -u(i - 1, j) - u(i + 1, j)) / (2 * dLx) ...
            + v(i, j) * (u(i, j + 1) - u(i, j - 1)) / (2 * dLy) + abs(v(i, j)) * (2 * u(i, j) -u(i , j - 1) - u(i, j + 1)) / (2 * dLy));
        
        advy(i, j) =  -(u(i, j) * (v(i + 1, j) - v(i - 1, j)) / (2 * dLx) + abs(u(i, j)) * (2 * v(i, j) -v(i - 1, j) - v(i + 1, j)) / (2 * dLx) ...
            + v(i, j) * (v(i, j + 1) - v(i, j - 1)) / (2 * dLy) + abs(v(i, j)) * (2 * v(i, j) -v(i , j - 1) -v(i, j + 1)) / (2 * dLy));
        
    end
end

end

function[usb, vsb] = midvel(usb, vsb, u, v, p)

% グローバル変数呼び出し
global nx ny dLx dLy dt rho

for i = 1 : nx - 1
    for j = 1 : ny - 2
        
        usb(i, j) =  (u(i + 1, j + 1) + u(i, j + 1)) / 2 - (dt / rho) * (p(i + 1, j + 1) - p(i, j + 1)) / dLx;
        
    end
end

for i = 1 : nx - 2
    for j = 1 : ny -1
        
        vsb(i, j) =  (v(i + 1, j + 1) + v(i + 1, j )) / 2 - (dt / rho) * (p(i + 1, j + 1) - p(i + 1, j )) / dLy;
        
    end
end

end

function[dp, ds] = poisson(usb, vsb, dp)

% グローバル変数呼び出し
global nx ny dLx dLy dt rho

% ポアソン方程式の係数定義
ac = zeros(nx, ny);
an = zeros(nx, ny);
ae = zeros(nx, ny);
as = zeros(nx, ny);
aw = zeros(nx, ny);

ac(:, :) = 2 * ((1/(dLx)^2) +(1/(dLy)^2));
an(:, :) = 1 / (dLy)^2;
ae(:, :) = 1 / (dLx)^2;
as(:, :) = 1 / (dLy)^2;
aw(:, :) = 1 / (dLx)^2;

% SOR法
ds = zeros(nx, ny);
eps = 10^(- 5);
ite_pmax = nx * ny * 20;% 反復回数
alpha = 1.7;% 緩和係数

for ite_p = 1 : ite_pmax% SOR法により圧力補正値を求める。
    
    error = 0;
    
    for i = 1 : nx
        for j = 1 : ny
            
            if i == 2 && j == 2% 圧力の基準点を固定しておく。
                
                dp(i, j) = 0;
                
            elseif  i == 1% 西側境界条件
                
                dp(i, j) =  dp(i + 1, j);
                
            elseif  i == nx% 東側境界条件
                
                dp(i, j) =  dp(nx - 1, j);
                
            elseif  j == 1% 南側境界条件
                
                dp(i, j) =  dp(i, j + 1);
                
            elseif  j == ny% 北側境界条件
                
                dp(i, j) =  dp(i, ny - 1);
                
            else% 内部領域
                
                ds(i, j) = (usb(i, j - 1) - usb(i - 1, j - 1))/dLx + (vsb(i - 1, j) - vsb(i - 1, j - 1))/dLy;
                dp_new= ( 1 / ac(i, j) ) * ...
                    ( ae(i, j) * dp(i + 1, j) + aw(i, j) * dp(i - 1, j) + an(i, j)* dp(i, j + 1) + as(i, j)* dp(i, j - 1) - (rho / dt) * ds(i, j) );
                error = max(abs(dp_new - dp(i, j)), error);
                dp(i, j) = dp(i, j) + alpha * (dp_new - dp(i, j));
                
            end
        end
    end
    
    if error < eps % 収束条件が満たされたらループを抜ける。
        break
    end
    
    if ite_p >= ite_pmax
        disp('最大反復回数に達しました。収束条件を満たしていません。');
    end
    
end

end

function[u, v]= velupdate(u, v, p)

% グローバル変数呼び出し
global nx ny dLx dLy dt rho

for i = 2 : nx - 1
    for j = 2 : ny - 1
        
        % セル中心速度には圧力は加えてなかったので、更新した圧力を加算する。（界面補間速度には更新前の圧力を加えていた。）
        u(i, j) = u(i, j) -(dt / rho) * ((p(i + 1, j) - p(i - 1, j))/(2 * dLx));
        v(i, j) = v(i, j) -(dt / rho) * ((p(i, j + 1) - p(i, j - 1))/(2 * dLy));
        
    end
end

end

function[] = vis_contour(filename, timestep, u, maxrange, minrange, fignum)

% グローバル変数呼び出し
global nx ny dt x y Lx Ly nu rho 

figure(fignum);
h = imagesc(x, y, flip(u));
xlabel('y')
ylabel('x')
view(270, 90); % 視点の設定
xticks(0 : round(Lx / 5, 1) : Lx);
yticks(0 : round(Ly / 5, 1) : Ly);

h.AlphaData = isfinite(u); % NaNやInfを透明にする
title(filename);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
axis equal;
axis tight;
axis on;
colorbar;
caxis([maxrange minrange])

% コメント
str_1 =['time = ', num2str(timestep * dt, '%.3f')];
str_2 =['Re = ', num2str(1 / nu, '%.f')];
str_3 =['rho = ', num2str(rho, '%.1f'),', nu = ', num2str(nu, '%.3f')];
str_4 =['nx = ', num2str(nx, '%.0f'),', ny = ', num2str(ny, '%.0f')];
str = {str_1, str_2, str_3, str_4};
text(0.1, 0.01, str)

% pngで保存
folder = [filename, '_re', num2str(1/nu),'_dt', num2str(dt),'_nx', num2str(nx),'_ny', num2str(ny)];

if timestep == 2
    mkdir(folder);
    filename2 = [folder, '/', filename, '_', num2str(timestep), '.png'];
    saveas(gcf, filename2)
elseif rem(timestep, 10) == 0
    filename2 = [folder, '/', filename, '_', num2str(timestep), '.png'];
    saveas(gcf, filename2)
end


end

function[] = vis_vector(filename, timestep, u, v, fignum)

% グローバル変数呼び出し
global dt nu nx ny

figure(fignum);
h = quiver(flipud(rot90(u)),flipud(rot90(v))) ;
h.Color = 'r';
h.AutoScaleFactor= 6;

title(['time = ', num2str(timestep * dt, '%.3f'),', Re = ', num2str(1 / nu, '%.f')]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
axis equal; axis tight; axis on;
xlim([1 nx]);
ylim([1 ny]);

% pngで保存
folder = [filename, '_re', num2str(1/nu),'_dt', num2str(dt),'_nx', num2str(nx),'_ny', num2str(ny)];

if timestep == 2
    mkdir(folder);
    filename2 = [folder, '/', filename, '_', num2str(timestep), '.png'];
    saveas(gcf, filename2)
elseif rem(timestep, 10) == 0
    filename2 = [folder, '/', filename, '_', num2str(timestep), '.png'];
    saveas(gcf, filename2)
end

end
