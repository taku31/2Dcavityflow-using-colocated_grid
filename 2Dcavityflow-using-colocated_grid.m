% ������
clear all;

% �O���[�o���ϐ��錾
global x y dt nx ny ite_max dLx dLy Lx Ly rho nu

% �p�����[�^�[
Lx = 1;% �v�Z�̈�
Ly = 1;% �v�Z�̈�
nx = 40;% �v�f��
ny = 40;% �v�f��
dLx = Lx / nx;% �v�f�T�C�Y
dLy = Ly / ny;% �v�f�T�C�Y
nu = 0.001;% ���S���W��
rho = 1;% ���x
t_max = 50;% SIM�v�Z����
dt = 0.001;% �^�C���X�e�b�v
ite_max = t_max / dt;% ������

% ���W�̐���
x = 0 : Lx / nx : Lx;
y = Ly : - Ly / ny : 0;
[X, Y] = meshgrid(x, y);

% �z��m��
u = zeros(nx, ny);% �Z�����S�ł̑��x
v = zeros(nx, ny);% �Z�����S�ł̑��x
p = zeros(nx, ny);% �Z�����S�ł̈���
dp = zeros(nx, ny);% �Z�����S�ł̕␳����
diffx = zeros(nx, ny);
diffy = zeros(nx, ny);
advx = zeros(nx, ny);
advy = zeros(nx, ny);
usb = zeros(nx - 1, ny - 2);% �Z���E�ʂł̒��ԑ��x
vsb = zeros(nx - 2, ny - 1);% �Z���E�ʂł̒��ԑ��x

% ���E����
bc.ue = 0;% ����
bc.us = 0;% �쑤
bc.uw = 0;% ����
bc.un = 1;% �k��
bc.ve = 0;% ����
bc.vs = 0;% �쑤
bc.vw = 0;% ����
bc.vn = 0;% �k��

% ���E�����̐ݒ�
[u, v] = setBC(u, v, bc);

% ���Ԕ��W
for ite = 2 : ite_max
    
    % �g�U���̌v�Z
    [diffx, diffy] = diffusion(u, v, diffx, diffy);
    
    % �ڗ����̌v�Z
    [advx, advy] = advection(u, v, advx, advy);
    
    % �ڗ��g�U�ɂ�鎞�Ԑi�s���v�Z
    u = u + (nu * diffx + advx) * dt;
    v = v + (nu * diffy + advy) * dt;
    
    % �Z���E�ʂł̒��ԑ��x���v�Z
    [usb, vsb] = midvel(usb, vsb, u, v, p);
    
    % �␳���͂̌v�Z
    [dp, ds] = poisson(usb, vsb, dp);
    
    % ���X�e�b�v���͂̌v�Z
    p = p + dp;
    
    % ���X�e�b�v�����̌v�Z
    [u, v] = velupdate(u, v, p);
    
    % ���E�����̐ݒ�
    [u, v] = setBC(u, v, bc);
    
    % �R�}���h�E�B���h�E�ւ̏o��
    txt = ['ite = ', num2str(ite), ' / ', num2str(ite_max)];
    disp(txt);
    
    % ���A���^�C������
%     vis_contour('u', ite, u, -0.8, 0.8, 1);
%     vis_contour('v', ite, v, -0.5, 0.5, 2);
%     vis_vector('vec', ite, u, v, 3)
    %vis_contour('p', ite, p, -3, 3, 4);
    %vis_contour('ds', ite, ds, -0.01, 0.01, 5);
    
end

%% �ȉ��֐�
function[u, v] = setBC(u, v, bc)

% �O���[�o���ϐ��Ăяo��
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

% �O���[�o���ϐ��Ăяo��
global nx ny dLx dLy

for i = 2 : nx - 1
    for j = 2 : ny - 1
        
        diffx(i, j) =  (u(i + 1, j) - 2 * u(i, j) + u(i - 1, j))/(dLx)^2 + (u(i , j + 1) - 2 * u(i, j) + u(i, j - 1))/(dLy)^2 ;
        diffy(i, j) =  (v(i + 1, j) - 2 * v(i, j) + v(i - 1, j))/(dLx)^2 + (v(i , j + 1) - 2 * v(i, j) + v(i, j - 1))/(dLy)^2 ;
        
    end
end

end

function[advx, advy] = advection(u, v, advx, advy)

% �O���[�o���ϐ��Ăяo��
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

% �O���[�o���ϐ��Ăяo��
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

% �O���[�o���ϐ��Ăяo��
global nx ny dLx dLy dt rho

% �|�A�\���������̌W����`
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

% SOR�@
ds = zeros(nx, ny);
eps = 10^(- 5);
ite_pmax = nx * ny * 20;% ������
alpha = 1.7;% �ɘa�W��

for ite_p = 1 : ite_pmax% SOR�@�ɂ�舳�͕␳�l�����߂�B
    
    error = 0;
    
    for i = 1 : nx
        for j = 1 : ny
            
            if i == 2 && j == 2% ���͂̊�_���Œ肵�Ă����B
                
                dp(i, j) = 0;
                
            elseif  i == 1% �������E����
                
                dp(i, j) =  dp(i + 1, j);
                
            elseif  i == nx% �������E����
                
                dp(i, j) =  dp(nx - 1, j);
                
            elseif  j == 1% �쑤���E����
                
                dp(i, j) =  dp(i, j + 1);
                
            elseif  j == ny% �k�����E����
                
                dp(i, j) =  dp(i, ny - 1);
                
            else% �����̈�
                
                ds(i, j) = (usb(i, j - 1) - usb(i - 1, j - 1))/dLx + (vsb(i - 1, j) - vsb(i - 1, j - 1))/dLy;
                dp_new= ( 1 / ac(i, j) ) * ...
                    ( ae(i, j) * dp(i + 1, j) + aw(i, j) * dp(i - 1, j) + an(i, j)* dp(i, j + 1) + as(i, j)* dp(i, j - 1) - (rho / dt) * ds(i, j) );
                error = max(abs(dp_new - dp(i, j)), error);
                dp(i, j) = dp(i, j) + alpha * (dp_new - dp(i, j));
                
            end
        end
    end
    
    if error < eps % �����������������ꂽ�烋�[�v�𔲂���B
        break
    end
    
    if ite_p >= ite_pmax
        disp('�ő唽���񐔂ɒB���܂����B���������𖞂����Ă��܂���B');
    end
    
end

end

function[u, v]= velupdate(u, v, p)

% �O���[�o���ϐ��Ăяo��
global nx ny dLx dLy dt rho

for i = 2 : nx - 1
    for j = 2 : ny - 1
        
        % �Z�����S���x�ɂ͈��͉͂����ĂȂ������̂ŁA�X�V�������͂����Z����B�i�E�ʕ�ԑ��x�ɂ͍X�V�O�̈��͂������Ă����B�j
        u(i, j) = u(i, j) -(dt / rho) * ((p(i + 1, j) - p(i - 1, j))/(2 * dLx));
        v(i, j) = v(i, j) -(dt / rho) * ((p(i, j + 1) - p(i, j - 1))/(2 * dLy));
        
    end
end

end

function[] = vis_contour(filename, timestep, u, maxrange, minrange, fignum)

% �O���[�o���ϐ��Ăяo��
global nx ny dt x y Lx Ly nu rho 

figure(fignum);
h = imagesc(x, y, flip(u));
xlabel('y')
ylabel('x')
view(270, 90); % ���_�̐ݒ�
xticks(0 : round(Lx / 5, 1) : Lx);
yticks(0 : round(Ly / 5, 1) : Ly);

h.AlphaData = isfinite(u); % NaN��Inf�𓧖��ɂ���
title(filename);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
axis equal;
axis tight;
axis on;
colorbar;
caxis([maxrange minrange])

% �R�����g
str_1 =['time = ', num2str(timestep * dt, '%.3f')];
str_2 =['Re = ', num2str(1 / nu, '%.f')];
str_3 =['rho = ', num2str(rho, '%.1f'),', nu = ', num2str(nu, '%.3f')];
str_4 =['nx = ', num2str(nx, '%.0f'),', ny = ', num2str(ny, '%.0f')];
str = {str_1, str_2, str_3, str_4};
text(0.1, 0.01, str)

% png�ŕۑ�
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

% �O���[�o���ϐ��Ăяo��
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

% png�ŕۑ�
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
