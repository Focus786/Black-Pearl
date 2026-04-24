function Simulation

%% =====================================================
% 0. 参数设定
%% =====================================================

clear; clc; close all;

params.stl_file   = 'ship.stl';
params.ship_length = 58;               %%   船长，单位cm
params.mass        = 2500;             %%   质量，单位g
params.rho         = 1.0;              %    水密度，单位g/cm³
params.CG_z        = 2.5;              %%   质心高度
params.dtheta      = 1;                %%   偏转步长
params.theta_min   = 0;                %%   最小偏转角
params.theta_max   = 70;               %%   最大偏转角

%% =====================================================
fprintf("\n====================================================\n");
fprintf("     SHIP STATIC STABILITY ANALYSIS START           \n");
fprintf("====================================================\n");

%% =====================================================
% 1. 加载
%% =====================================================
fprintf("\n========== STEP 1 ==========\n");

TR = stlread(params.stl_file);
F  = double(TR.ConnectivityList);
V  = double(TR.Points);

V = V(all(isfinite(V),2),:);
[V,~,idx] = unique(V,'rows');
F = reshape(idx(F(:)), size(F));

fprintf("Mesh loaded: V=%d F=%d\n", size(V,1), size(F,1));

%% =====================================================
% 2. SCALE
%% =====================================================
fprintf("\n========== STEP 2 ==========\n");

Lx    = max(V(:,1)) - min(V(:,1));
scale = params.ship_length / Lx;
V     = V * scale;
V(:,3)= V(:,3) - min(V(:,3));

fprintf("Scaled length = %.2f cm\n", max(V(:,1))-min(V(:,1)));

%% =====================================================
% 3. 统一法向朝外（修复负体积）
%% =====================================================
fprintf("\n========== STEP 3 - Normal Fix ==========\n");

p1n = V(F(:,1),:);
p2n = V(F(:,2),:);
p3n = V(F(:,3),:);
Cn  = (p1n+p2n+p3n)/3;
Nn  = cross(p2n-p1n, p3n-p1n, 2);
V_center    = mean(V,1);
flip_mask   = dot(Cn - V_center, Nn, 2) < 0;
F(flip_mask,2:3) = F(flip_mask,[3,2]);
fprintf("Normals fixed: %d faces flipped\n", sum(flip_mask));

%% =====================================================
% 4. CG
%% =====================================================
fprintf("\n========== STEP 4 ==========\n");

p1 = V(F(:,1),:);
p2 = V(F(:,2),:);
p3 = V(F(:,3),:);

C  = (p1+p2+p3)/3;
A  = 0.5 * vecnorm(cross(p2-p1,p3-p1,2),2,2);

CG_xy = sum(C(:,1:2).*A,1) / sum(A);
CG    = [CG_xy, params.CG_z];
CG(2) = 0;  % 横向对称

fprintf("CG = [%.2f %.2f %.2f]\n", CG);

%% =====================================================
% 5. 坐标系移到CG
%% =====================================================
fprintf("\n========== STEP 5 ==========\n");

V = V - CG;
C = C - CG;

fprintf("Coordinate shifted to CG origin\n");

%% =====================================================
% 6. 初始吃水
%% =====================================================
fprintf("\n========== STEP 6 ==========\n");

subVol0 = @(h) submerged_volume_exact(V,F,h) - params.mass/params.rho;
h0 = fzero(subVol0, 0);
fprintf("Equilibrium draft h0 = %.3f cm\n", h0);

%% =====================================================
% 7. 角度测试
%% =====================================================
fprintf("\n========== STEP 7 ==========\n");

theta_range = params.theta_min:params.dtheta:params.theta_max;
n = length(theta_range);

GZ   = nan(1,n);
RM   = nan(1,n);
Vsub_arr = nan(1,n);
h_arr    = nan(1,n);
CB_arr   = nan(n,3);

fprintf("\n----------------------------------------------------\n");
fprintf(" θ(°) |   h(cm) |   Vsub   |   GZ(cm)  |  RM(g·cm)\n");
fprintf("----------------------------------------------------\n");

for i = 1:n

    theta = theta_range(i);
    R     = rotX(theta);

    Vrot = (R * V')';

    % 重新求解平衡吃水
    subVol = @(h) submerged_volume_exact(Vrot,F,h) - params.mass/params.rho;
    try
        h_theta = fzero(subVol, h0);
    catch
        fprintf(" %3.0f |   FAIL  |    -    |    -    |    -\n", theta);
        continue;
    end

    % 进水判断（水面超过旋转后船体最高点）
    z_top = max(Vrot(:,3));
    if h_theta >= z_top * 0.99
        fprintf(" %3.0f | FLOODED (h=%.2f >= top=%.2f), stopping\n", ...
            theta, h_theta, z_top);
        break;
    end

    [Vsub, CB] = submerged_volume_exact(Vrot,F,h_theta);

    % CG在新坐标系下是原点，旋转后
    CGrot = (R * [0;0;0])';  % 始终是原点

    dy = CB(2) - CGrot(2);
    dz = CB(3) - CGrot(3);

    gz = -(dy*cosd(theta) - dz*sind(theta));

    if theta == 0
        gz = 0;
    end

    GZ(i)       = gz;
    RM(i)       = params.rho * Vsub * gz;
    Vsub_arr(i) = Vsub;
    h_arr(i)    = h_theta;
    CB_arr(i,:) = CB;

    fprintf(" %3.0f | %6.3f | %7.1f | %7.3f | %9.2f\n", ...
        theta, h_theta, Vsub, gz, RM(i));
end

fprintf("----------------------------------------------------\n");

%% =====================================================
% 8. 绘图
%% =====================================================

scr = get(0,'ScreenSize');
fig = uifigure('Name','Ship Stability Results', ...
    'Position', [0 0 scr(3) scr(4)]);
tg  = uitabgroup(fig, 'Units','normalized', 'Position',[0 0 1 1]);

%% 每个角度一个3D Tab
for i = 1:n

    if isnan(GZ(i)); continue; end

    theta = theta_range(i);
    R     = rotX(theta);
    Vrot  = (R * V')';
    h_theta = h_arr(i);
    CB    = CB_arr(i,:);

    tab = uitab(tg,'Title',sprintf('%.0f°',theta));
    ax = uiaxes(tab, 'Units','normalized', 'Position',[0.05 0.05 0.9 0.9]);
    hold(ax,'on');

    % 船体
    trisurf(triangulation(F,Vrot), ...
        'Parent',ax, ...
        'FaceColor',[0.6 0.8 1], ...
        'FaceAlpha',0.35, ...
        'EdgeColor','none');

    % CG（原点）
    plot3(ax,0,0,0,'ro','MarkerFaceColor','r','MarkerSize',8);

    % CB
    plot3(ax,CB(1),CB(2),CB(3),'bo','MarkerFaceColor','b','MarkerSize',8);

    % 水面（足够大覆盖船体）
    margin = max(abs(V(:,2))) * 1.5;
    [Xw,Yw] = meshgrid(...
        linspace(min(Vrot(:,1)),max(Vrot(:,1)),15), ...
        linspace(-margin, margin, 15));
    surf(ax,Xw,Yw,h_theta*ones(size(Xw)), ...
        'FaceColor',[0 0.8 1],'FaceAlpha',0.25,'EdgeColor','none');

    axis(ax,'equal'); grid(ax,'on');
    view(ax,35,20);
    camlight(ax); lighting(ax,'phong');
    rotate3d(ax,'on');
    title(ax,sprintf('Heel = %.0f°  |  h=%.2f  |  GZ=%.3f', ...
        theta, h_theta, GZ(i)));
    legend(ax,{'Hull','CG','CB','Water'},'Location','best');
end

%% GZ Tab
tab_gz = uitab(tg,'Title','GZ Curve');
ax_gz  = uiaxes(tab_gz,  'Units','normalized', 'Position',[0.05 0.05 0.9 0.9]);
valid  = ~isnan(GZ);
plot(ax_gz, theta_range(valid), GZ(valid), 'k-o','LineWidth',1.5);
yline(ax_gz,0,'r--');
grid(ax_gz,'on');
title(ax_gz,'GZ Curve (Righting Arm)');
xlabel(ax_gz,'Heel Angle (°)');
ylabel(ax_gz,'GZ (cm)');

%% RM Tab
tab_rm = uitab(tg,'Title','RM Curve');
ax_rm  = uiaxes(tab_rm,  'Units','normalized', 'Position',[0.05 0.05 0.9 0.9]);
plot(ax_rm, theta_range(valid), RM(valid), 'm-o','LineWidth',1.5);
yline(ax_rm,0,'r--');
grid(ax_rm,'on');
title(ax_rm,'Righting Moment');
xlabel(ax_rm,'Heel Angle (°)');
ylabel(ax_rm,'RM (g·cm)');

%% Summary Tab
[RMmax, i_max] = max(RM(valid));
theta_stable   = theta_range(find(valid,1) + i_max - 1);

idx_limit = find(RM(2:end)<=0,1,'first')+1;
if isempty(idx_limit)
    theta_limit = NaN;
else
    theta_limit = theta_range(idx_limit);
end

tab_sum = uitab(tg,'Title','Summary');
ax_sum = uiaxes(tab_sum, 'Units','normalized', 'Position',[0.05 0.05 0.9 0.9]);
axis(ax_sum,'off');
text(ax_sum,0.1,0.85,sprintf('Draft h0       = %.3f cm', h0),      'FontSize',13);
text(ax_sum,0.1,0.70,sprintf('Max RM angle   = %.2f°',  theta_stable),'FontSize',13);
text(ax_sum,0.1,0.55,sprintf('Max RM         = %.2f g·cm', RMmax),  'FontSize',13);
if isnan(theta_limit)
    text(ax_sum,0.1,0.40,'Limit angle    = NOT REACHED','FontSize',13);
else
    text(ax_sum,0.1,0.40,sprintf('Limit angle    = %.2f°', theta_limit),'FontSize',13);
end

%% =====================================================
% 9. 总结
%% =====================================================
fprintf("\n====================================================\n");
fprintf("                 STABILITY SUMMARY                 \n");
fprintf("====================================================\n");
fprintf("Draft h0       : %.3f cm\n",   h0);
fprintf("Max RM angle   : %.2f°\n",     theta_stable);
fprintf("Max RM         : %.2f g·cm\n", RMmax);
if isnan(theta_limit)
    fprintf("Limit angle    : NOT REACHED\n");
else
    fprintf("Limit angle    : %.2f°\n", theta_limit);
end
fprintf("====================================================\n");

end

%% =====================================================
% 函数
%% =====================================================

function R = rotX(t)
t = deg2rad(t);
R = [1 0 0;
     0 cos(t) -sin(t);
     0 sin(t)  cos(t)];
end

function [Vsub, CB] = submerged_volume_exact(V,F,h)

Vsub = 0;
M    = [0 0 0];

for i = 1:size(F,1)
    tri = V(F(i,:),:);
    [poly,ok] = clip_triangle_z(tri,h);
    if ~ok; continue; end
    [v,c] = poly_tetra_volume(poly);
    Vsub = Vsub + v;
    M    = M + c*v;
end

Vsub = abs(Vsub);

if Vsub > 1e-9
    CB = M / Vsub;
else
    CB = [0 0 0];
end

end

function [poly,ok] = clip_triangle_z(tri,h)

inside = tri(:,3) <= h;

if all(inside);  poly=tri; ok=true;  return; end
if all(~inside); poly=[];  ok=false; return; end

pts = [];
for i = 1:3
    j  = mod(i,3)+1;
    p1 = tri(i,:);
    p2 = tri(j,:);
    if p1(3) <= h
        pts = [pts; p1];
    end
    if (p1(3)-h)*(p2(3)-h) < 0
        t   = (h-p1(3))/(p2(3)-p1(3));
        pts = [pts; p1+t*(p2-p1)];
    end
end

poly = pts;
ok   = size(poly,1) >= 3;

end

function [vol,cent] = poly_tetra_volume(poly)

if size(poly,1) < 3
    vol=0; cent=[0 0 0]; return;
end

p0   = poly(1,:);
vol  = 0;
cent = [0 0 0];

for i = 2:size(poly,1)-1
    a = poly(i,:);
    b = poly(i+1,:);
    v = dot(p0,cross(a,b))/6;
    vol  = vol  + v;
    cent = cent + v*(p0+a+b)/4;
end

if abs(vol) > 1e-12
    cent = cent/vol;
else
    cent = [0 0 0];
end

end