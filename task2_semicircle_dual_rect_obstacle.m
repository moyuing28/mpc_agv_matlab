% 任务2：半圆弧全局轨迹 + 双矩形风险通道避障 + 时变误差状态MPC
%
% 结构：
% 半圆弧全局参考
%      ↓
% NORMAL / AVOID / RECOVER 状态机
%      ↓
% AVOID：fmincon 优化“相对名义曲率角速度”的增量 δω
% RECOVER：Frenet 坐标下三次 Hermite 横向回归
%      ↓
% 多速率插值 T -> dt
%      ↓
% 下层 TV error-state MPC（含 d_k 和 -Delta(u_ref)）
%
% 任务2相对任务1最关键的变化：
% 1) 参考线由直线变成半圆弧，omega_ref = V/R 不再为0；
% 2) 所有“横向/纵向”判断改到局部 Frenet 坐标；
% 3) corridor 不再是固定 y=c，而是沿曲线法向偏移 e_corridor；
% 4) 局部规划优化变量改为 δω = omega - omega_ref，
%    控制时域以后 δω=0，即恢复名义圆弧曲率，而不是 omega=0；
% 5) RECOVER 在曲线 Frenet 坐标中生成平滑回归轨迹。

% ！！！这个代码完全AI改的，延续的是task1_dual_rect_obstacle_v3_fast_recover的思路，
% AI似乎加入了一些改进但是整体思路不变
clear; close all; clc;

%% ======================== 1. 基本参数 ========================
L = 0.5;                 % 轮距 (m)
dt = 0.01;               % 下层MPC采样时间
T  = 0.05;               % 上层规划采样时间

V_target = 0.8;          % 巡航速度 (m/s)

Np = 60;                 % 下层MPC预测步长：0.6 s
Nc = 30;                 % 下层MPC控制步长

plan_Np = 20;            % 上层规划预测步长：1.0 s
plan_Nc = 10;            % 上层规划控制步长：0.5 s

Nx = 3;                  % x,y,theta
Nu = 2;                  % v,omega

O = zeros(Nu,Nx);
E = eye(Nu);

%% ======================== 2. 半圆弧全局参考 ========================
% 选择一个非常简单的右半圆：
%
% 起点 (0,0), theta=0
%      →
%       )
%      )
%     )
% 终点 (0,6), theta=pi
%
R_circle = 3.0;
circle_center = [0, R_circle];

phi0   = -pi/2;
phi_end = pi/2;

arc_length = R_circle*(phi_end-phi0);
arc_time   = arc_length/V_target;
T_max = arc_time + 1.5;

% 初始状态直接位于半圆起点
x = 0;
y = 0;
theta = 0;

v = V_target;
omega = V_target/R_circle;

%% ======================== 3. 双矩形障碍物 ========================
% 把障碍物放在半圆中点 phi=0 附近。
%
% 半圆中点：p=(3,3)，切向竖直向上。
% 局部 Frenet 法向指向圆心（左侧）。
%
% 高风险障碍物：位于全局曲线外侧，离参考轨迹近；
% 低风险障碍物：位于曲线内侧，离参考轨迹更远。
%
phi_obs_anchor = 0;
p_obs_anchor = circle_center + R_circle*[cos(phi_obs_anchor), sin(phi_obs_anchor)];
theta_obs_anchor = phi_obs_anchor + pi/2;
tangent_anchor = [cos(theta_obs_anchor), sin(theta_obs_anchor)];
normal_anchor  = [-sin(theta_obs_anchor), cos(theta_obs_anchor)];

% Frenet横向中心偏移：
% e<0：曲线外侧（高风险）
% e>0：曲线内侧（低风险）
e_obs_high = -0.38;
e_obs_low  =  0.85;

p_high = p_obs_anchor + e_obs_high*normal_anchor;
p_low  = p_obs_anchor + e_obs_low *normal_anchor;

% 因为中点切向为竖直方向，所以矩形：
% x方向宽度 = 0.36（法向）
% y方向长度 = 1.00（切向）
obs_rect = [ ...
    p_high(1), p_high(2), 0.36, 1.00; ...
    p_low(1),  p_low(2),  0.36, 1.00];

Nobs = size(obs_rect,1);

vehicle_radius = 0.16;
safe_margin    = 0.10;
required_clearance = vehicle_radius + safe_margin;

d_influence = 0.70;

% 对每个矩形计算它相对半圆参考轨迹的 Frenet 信息
obs_s = zeros(Nobs,1);
obs_e = zeros(Nobs,1);
obs_half_tangent = zeros(Nobs,1);
obs_half_normal  = zeros(Nobs,1);
obs_path_clearance = zeros(Nobs,1);

for j = 1:Nobs
    [obs_s(j), obs_e(j), theta_path_j] = ...
        obstacle_frenet_info(obs_rect(j,:), circle_center, R_circle, phi0, phi_end);

    tj = [cos(theta_path_j), sin(theta_path_j)];
    nj = [-sin(theta_path_j), cos(theta_path_j)];

    hx = obs_rect(j,3)/2;
    hy = obs_rect(j,4)/2;

    obs_half_tangent(j) = abs(tj(1))*hx + abs(tj(2))*hy;
    obs_half_normal(j)  = abs(nj(1))*hx + abs(nj(2))*hy;

    obs_path_clearance(j) = ...
        abs(obs_e(j)) - obs_half_normal(j) - vehicle_radius;
end

% 风险权重：越靠近全局曲线，权重越大
obs_risk_weight = 1 + 4*exp(-max(obs_path_clearance,0)/0.25);

% 计算安全通道在 Frenet 法向上的中心偏移
neg_idx = find(obs_e < 0,1,'first');
pos_idx = find(obs_e > 0,1,'first');

if isempty(neg_idx) || isempty(pos_idx)
    error('任务2需要两块障碍物位于参考曲线法向两侧。');
end

% 障碍物加上车辆半径和安全裕量之后，面向通道的边界
neg_inner_edge = obs_e(neg_idx) + obs_half_normal(neg_idx) + required_clearance;
pos_inner_edge = obs_e(pos_idx) - obs_half_normal(pos_idx) - required_clearance;

if pos_inner_edge <= neg_inner_edge
    error('两个矩形之间的安全通道过窄，请增大障碍物间距或减小安全裕量。');
end

corridor_e = 0.5*(neg_inner_edge + pos_inner_edge);

% 沿弧长方向的通道作用范围
corridor_s_min = min(obs_s - obs_half_tangent) - 0.65;
corridor_s_max = max(obs_s + obs_half_tangent) + 0.65;

%% ======================== 4. 状态机与权重 ========================
normal_mode  = 1;
avoid_mode   = 2;
recover_mode = 3;

mode = normal_mode;
recover_side = 0;

% 先保持任务1比较稳妥的触发范围；后面只调这个就能改变“提前多久避障”
d_trigger_long = 0.85;       % 沿弧长距障碍物前缘的触发距离
d_trigger_lat  = 0.15;       % 只有真正靠近全局轨迹的障碍物负责触发
pass_margin    = 0.20;

omega_seed = zeros(plan_Nc,1);   % 这里表示 δω 的热启动

% 下层MPC权重
Q_single = diag([15,10,5]);
Q = kron(eye(Np),Q_single);

R_single = diag([0.5,0.05]);
R = kron(eye(Nc),R_single);

Row = 10;

% 上层局部规划权重
Q_plan = diag(linspace(40,120,plan_Np));
R_plan = 12*eye(plan_Nc);

S_avoid = 70;
Q_corridor = 500;
Q_terminal_avoid = 350;
Q_terminal_theta_avoid = 450;

% 绘图
draw = true;
plan_draw = true;

%% ======================== 5. 历史数据 ========================
history.x = [];
history.y = [];
history.theta = [];
history.v = [];
history.omega = [];
history.mode = [];
history.t = [];
history.clearance = [];
history.plan_x = {};
history.plan_y = {};

history_ref.x = [];
history_ref.y = [];

plan_num = 1;
num = 1;

kesi = zeros(Nx+Nu,1);

% 下层输入增量累加矩阵（固定，可提前构造）
A_t = tril(ones(Nc));
A_I = kron(A_t,eye(Nu));

%% ======================== 6. 主仿真循环 ========================
for t = 0:dt:T_max

    %% ---------- 6.1 生成全局半圆参考 ----------
    if mode == normal_mode
        pre_t = Np*dt;
    else
        pre_t = plan_Np*T;
    end

    t_global = t + (0:dt:pre_t)';

    [global_x,global_y,global_theta,global_v,global_omega,global_s] = ...
        semicircle_reference(t_global,circle_center,R_circle,phi0,phi_end,V_target);

    x_ref = global_x(1);
    y_ref = global_y(1);
    theta_ref = global_theta(1);

    %% ---------- 6.2 根据模式生成下层最终参考 ----------
    if mode == normal_mode

        ctrl_ref_x     = global_x(1:Np+1);
        ctrl_ref_y     = global_y(1:Np+1);
        ctrl_ref_theta = unwrap(global_theta(1:Np+1));
        ctrl_ref_v     = global_v(1:Np+1);
        ctrl_ref_omega = global_omega(1:Np+1);

    else
        % 把全局参考从 dt 重采样到规划器 T
        t_old = (0:length(global_x)-1)'*dt;
        t_new = (1:plan_Np)'*T;

        Xref = interp1(t_old,global_x,t_new,'linear');
        Yref = interp1(t_old,global_y,t_new,'linear');
        THETAref = interp1(t_old,unwrap(global_theta),t_new,'linear');
        Vref_plan = interp1(t_old,global_v,t_new,'linear');
        OMEGAref_plan = interp1(t_old,global_omega,t_new,'linear');
        Sref = interp1(t_old,global_s,t_new,'linear');

        if mode == avoid_mode
            %% -------- AVOID：优化 δω，而不是直接优化 omega --------
            omega_max = 4;
            omega_min = -4;

            % δω 的上下界，由实际 omega 限制换算得到
            lb_plan = omega_min - OMEGAref_plan(1:plan_Nc);
            ub_plan = omega_max - OMEGAref_plan(1:plan_Nc);

            options_plan = optimoptions('fmincon', ...
                'Algorithm','sqp', ...
                'Display','none', ...
                'MaxIterations',45, ...
                'MaxFunctionEvaluations',3000);

            % 实际角速度变化率约束：
            % actual_omega = omega_ref + delta_omega
            delta_omega_plan = 0.35;

            Domega = eye(plan_Nc);
            Domega(2:end,1:end-1) = ...
                Domega(2:end,1:end-1) - eye(plan_Nc-1);

            omega_ref_ctrl = OMEGAref_plan(1:plan_Nc);
            first_prev = [omega; zeros(plan_Nc-1,1)];

            Aineq_plan = [ Domega; -Domega];
            bineq_plan = [ ...
                delta_omega_plan*ones(plan_Nc,1) + first_prev - Domega*omega_ref_ctrl; ...
                delta_omega_plan*ones(plan_Nc,1) - first_prev + Domega*omega_ref_ctrl];

            % 终端约束1：最后的 δω 回到0，即重新恢复名义圆弧曲率
            %
            % 终端约束2：控制时域末端航向与全局圆弧参考航向一致
            Aeq_plan = [zeros(1,plan_Nc-1),1; ...
                        T*ones(1,plan_Nc)];

            theta_target = THETAref(plan_Nc);
            nominal_heading_change = T*sum(OMEGAref_plan(1:plan_Nc));

            beq_plan = [ ...
                0; ...
                theta_target - theta - nominal_heading_change];

            theta_limit = deg2rad(35);

            nonlcon = @(delta_target)rect_clearance_constraints_curve( ...
                delta_target,x,y,theta,v,plan_Np,plan_Nc,T, ...
                OMEGAref_plan,THETAref,obs_rect,vehicle_radius,safe_margin,theta_limit);

            [delta_omega_opt,~,exitflag_plan] = fmincon( ...
                @(delta_target)function_avoid_curve( ...
                    delta_target,x,y,theta,v,plan_Np,plan_Nc,T, ...
                    Xref,Yref,THETAref,OMEGAref_plan,Sref, ...
                    Q_plan,R_plan,S_avoid,Q_terminal_avoid,Q_terminal_theta_avoid, ...
                    obs_rect,vehicle_radius,d_influence,obs_risk_weight, ...
                    corridor_e,Q_corridor,corridor_s_min,corridor_s_max), ...
                omega_seed,Aineq_plan,bineq_plan,Aeq_plan,beq_plan, ...
                lb_plan,ub_plan,nonlcon,options_plan);

            if exitflag_plan <= 0 || any(~isfinite(delta_omega_opt))
                warning('曲线局部避障规划器在 t=%.2f s 未正常收敛，沿用上一拍 δω 序列。',t);
                delta_omega_opt = omega_seed;
            end

            % warm start
            omega_seed = [delta_omega_opt(2:end);0];

            [x_pre,y_pre,theta_pre,omega_pre] = rollout_plan_curve( ...
                delta_omega_opt,x,y,theta,v,plan_Np,plan_Nc,T,OMEGAref_plan);

        else
            %% -------- RECOVER：曲线 Frenet + Hermite 回归 --------
            [x_pre,y_pre,theta_pre,omega_pre] = generate_recover_plan_curve( ...
                x,y,theta,v,plan_Np,T, ...
                Xref,Yref,THETAref,x_ref,y_ref,theta_ref);

            omega_seed(:) = 0;
        end

        history.plan_x{plan_num,1} = x_pre;
        history.plan_y{plan_num,1} = y_pre;
        plan_num = plan_num + 1;

        %% -------- 多速率插值：规划器 T -> 下层 dt --------
        t_plan = [0;(1:plan_Np)'*T];

        x_plan = [x;x_pre];
        y_plan = [y;y_pre];
        theta_plan = unwrap([theta;theta_pre]);

        t_ctrl = (0:Np)'*dt;

        ctrl_ref_x = interp1(t_plan,x_plan,t_ctrl,'linear');
        ctrl_ref_y = interp1(t_plan,y_plan,t_ctrl,'linear');
        ctrl_ref_theta = interp1(t_plan,theta_plan,t_ctrl,'linear');

        ctrl_ref_v = V_target*ones(Np+1,1);

        ctrl_ref_omega = zeros(Np+1,1);
        ctrl_ref_omega(1:Np) = diff(ctrl_ref_theta)/dt;
        ctrl_ref_omega(Np+1) = ctrl_ref_omega(Np);
    end

    %% ---------- 6.3 d_k 与参考输入变化 ----------
    d_ref = zeros(Nx,Np);
    delta_u_ref = zeros(Nu,Np);

    for i = 1:Np
        xr_i = ctrl_ref_x(i);
        yr_i = ctrl_ref_y(i);
        theta_r_i = ctrl_ref_theta(i);

        vr_i = ctrl_ref_v(i);
        wr_i = ctrl_ref_omega(i);

        xr_model_next = xr_i + dt*vr_i*cos(theta_r_i);
        yr_model_next = yr_i + dt*vr_i*sin(theta_r_i);
        theta_model_next = theta_r_i + dt*wr_i;

        d_ref(:,i) = [ ...
            xr_model_next - ctrl_ref_x(i+1); ...
            yr_model_next - ctrl_ref_y(i+1); ...
            wrap_pi(theta_model_next - ctrl_ref_theta(i+1))];

        delta_u_ref(:,i) = [ ...
            ctrl_ref_v(i+1)-ctrl_ref_v(i); ...
            ctrl_ref_omega(i+1)-ctrl_ref_omega(i)];
    end

    %% ---------- 6.4 当前增广误差 ----------
    kesi(1) = x - ctrl_ref_x(1);
    kesi(2) = y - ctrl_ref_y(1);
    kesi(3) = wrap_pi(theta - ctrl_ref_theta(1));
    kesi(4) = v - ctrl_ref_v(1);
    kesi(5) = omega - ctrl_ref_omega(1);

    %% ---------- 6.5 时变线性模型 ----------
    A_cell = cell(Np,1);
    B_cell = cell(Np,1);

    C_tilde = [ ...
        1 0 0 0 0; ...
        0 1 0 0 0; ...
        0 0 1 0 0];

    for j = 1:3
        A_cell{j} = [ ...
            1,0,-dt*V_target*sin(theta); ...
            0,1, dt*V_target*cos(theta); ...
            0,0,1];

        B_cell{j} = [ ...
            dt*cos(theta),0; ...
            dt*sin(theta),0; ...
            0,dt];
    end

    for i = 4:Np
        theta_r_i = ctrl_ref_theta(i);

        A_cell{i} = [ ...
            1,0,-dt*V_target*sin(theta_r_i); ...
            0,1, dt*V_target*cos(theta_r_i); ...
            0,0,1];

        B_cell{i} = [ ...
            dt*cos(theta_r_i),0; ...
            dt*sin(theta_r_i),0; ...
            0,dt];
    end

    %% ---------- 6.6 构造 THETA ----------
    THETA_cell = cell(Np,Nc);

    for i = 1:Np
        for j = 1:Nc
            if j <= i
                A_prod = eye(Nx+Nu);

                for k = (j+1):i
                    A_tilde_k = [A_cell{k},B_cell{k};O,E];
                    A_prod = A_tilde_k*A_prod;
                end

                B_tilde_j = [B_cell{j};E];
                THETA_cell{i,j} = C_tilde*A_prod*B_tilde_j;
            else
                THETA_cell{i,j} = zeros(Nx,Nu);
            end
        end
    end

    THETA = cell2mat(THETA_cell);

    %% ---------- 6.7 含 d_k 的自由误差传播 ----------
    Xi_free = kesi;
    error_free = zeros(Nx*Np,1);

    for i = 1:Np
        A_tilde_i = [A_cell{i},B_cell{i};O,E];

        g_i = [ ...
            d_ref(:,i); ...
            -delta_u_ref(:,i)];

        Xi_free = A_tilde_i*Xi_free + g_i;

        error_free((i-1)*Nx+1:i*Nx) = C_tilde*Xi_free;
    end

    %% ---------- 6.8 QP目标函数 ----------
    H_cell = cell(2,2);
    H_cell{1,1} = THETA'*Q*THETA + R;
    H_cell{1,2} = zeros(Nu*Nc,1);
    H_cell{2,1} = zeros(1,Nu*Nc);
    H_cell{2,2} = Row;

    H = cell2mat(H_cell);
    H = (H+H')/2;

    f_cell = cell(1,2);
    f_cell{1,1} = error_free'*Q*THETA;
    f_cell{1,2} = 0;

    f = cell2mat(f_cell)';

    %% ---------- 6.9 下层控制约束 ----------
    U = [v;omega];

    umin = [0;-4.0];
    umax = [1.2;4.0];

    delta_umin = [-0.2;-0.7];
    delta_umax = [ 0.2; 0.7];

    Umin = kron(ones(Nc,1),umin);
    Umax = kron(ones(Nc,1),umax);
    Ut = kron(ones(Nc,1),U);

    A_cons = [ ...
        A_I, zeros(Nu*Nc,1); ...
       -A_I, zeros(Nu*Nc,1)];

    b_cons = [ ...
        Umax-Ut; ...
       -Umin+Ut];

    delta_Umin = kron(ones(Nc,1),delta_umin);
    delta_Umax = kron(ones(Nc,1),delta_umax);

    lb = [delta_Umin;0];
    ub = [delta_Umax;10];

    options_qp = optimset('Algorithm','interior-point-convex','Display','off');

    [X,~,exitflag_qp] = quadprog(H,f,A_cons,b_cons,[],[],lb,ub,[],options_qp);

    if exitflag_qp <= 0 || any(~isfinite(X))
        warning('下层QP在 t=%.2f s 未正常求解，本拍保持控制量。',t);
        d_v = 0;
        d_omega = 0;
    else
        d_v = X(1);
        d_omega = X(2);
    end

    v = v + d_v;
    omega = omega + d_omega;

    %% ---------- 6.10 AGV实际运动学 ----------
    vr = v + omega*L/2;
    vl = v - omega*L/2;

    v_act = (vr+vl)/2;
    w_act = (vr-vl)/L;

    x = x + v_act*cos(theta)*dt;
    y = y + v_act*sin(theta)*dt;
    theta = wrap_pi(theta + w_act*dt);

    %% ---------- 6.11 Frenet状态与碰撞距离 ----------
    [s_vehicle,e_lat,theta_path_vehicle] = ...
        circle_frenet_projection([x,y],circle_center,R_circle,phi0,phi_end);

    e_theta_path = wrap_pi(theta-theta_path_vehicle);

    rect_clearance = zeros(Nobs,1);
    for j = 1:Nobs
        rect_clearance(j) = ...
            point_rect_signed_distance([x,y],obs_rect(j,:)) - vehicle_radius;
    end

    %% ---------- 6.12 状态机 ----------
    switch mode
        case normal_mode
            should_avoid = false;

            for j = 1:Nobs
                % 车辆到障碍物“前缘”的沿弧长距离
                d_long = obs_s(j) - s_vehicle - obs_half_tangent(j);

                obstacle_ahead = ...
                    d_long > -pass_margin;

                % 只有确实靠近全局曲线的高风险障碍物负责触发
                lateral_risk = ...
                    obs_path_clearance(j) <= d_trigger_lat;

                if obstacle_ahead && ...
                   lateral_risk && ...
                   d_long <= d_trigger_long

                    should_avoid = true;
                end
            end

            if should_avoid
                mode = avoid_mode;

                % 当前车辆相对 corridor 的方向
                turn_sign = sign(corridor_e-e_lat);
                if turn_sign == 0
                    turn_sign = 1;
                end

                % S形的 δω 初值；正负成对，末端回到0
                omega_seed(:) = 0;

                seed_shape = [ ...
                    0.20;0.40;0.55;0.55;0.35; ...
                   -0.20;-0.45;-0.45;-0.25;0];

                nseed = min(length(seed_shape),plan_Nc);
                omega_seed(1:nseed) = turn_sign*seed_shape(1:nseed);
                omega_seed(end) = 0;
            end

        case avoid_mode
            passed_all = true;

            for j = 1:Nobs
                if s_vehicle < ...
                        obs_s(j)+obs_half_tangent(j)+pass_margin
                    passed_all = false;
                end
            end

            if passed_all
                mode = recover_mode;

                recover_side = sign(e_lat);
                if recover_side == 0
                    recover_side = sign(corridor_e);
                end
                if recover_side == 0
                    recover_side = 1;
                end

                omega_seed(:) = 0;
            end

        case recover_mode
            recover_converged = ...
                abs(e_lat) <= 0.07 && ...
                abs(e_theta_path) <= deg2rad(10);

            recover_crossed = ...
                recover_side ~= 0 && ...
                e_lat*recover_side <= 0 && ...
                abs(e_theta_path) <= deg2rad(15);

            if recover_converged || recover_crossed
                mode = normal_mode;
                recover_side = 0;
                omega_seed(:) = 0;
            end
    end

    %% ---------- 6.13 保存 ----------
    history_ref.x(num) = x_ref;
    history_ref.y(num) = y_ref;

    history.x(num) = x;
    history.y(num) = y;
    history.theta(num) = theta;
    history.v(num) = v;
    history.omega(num) = omega;
    history.mode(num) = mode;
    history.t(num) = t;
    history.clearance(num) = min(rect_clearance);

    num = num + 1;
end

%% ======================== 7. 绘图 ========================
if draw
    % 完整半圆参考
    t_plot = linspace(0,arc_time,800)';
    [x_plot,y_plot,~,~,~,~] = ...
        semicircle_reference(t_plot,circle_center,R_circle,phi0,phi_end,V_target);

    figure;
    plot(x_plot,y_plot,'LineWidth',1.4,'DisplayName','半圆弧全局参考');
    hold on;
    plot(history.x,history.y,'LineWidth',1.3,'DisplayName','MPC实际轨迹');

    if plan_draw
        for k = 1:length(history.plan_x)
            plot(history.plan_x{k},history.plan_y{k}, ...
                'LineWidth',0.5,'HandleVisibility','off');
        end
    end

    for j = 1:Nobs
        cx = obs_rect(j,1);
        cy = obs_rect(j,2);
        w = obs_rect(j,3);
        h = obs_rect(j,4);

        xr = [cx-w/2,cx+w/2,cx+w/2,cx-w/2];
        yr = [cy-h/2,cy-h/2,cy+h/2,cy+h/2];

        if j == 1
            patch(xr,yr,[0.25 0.25 0.25], ...
                'EdgeColor','k','LineWidth',1.2, ...
                'DisplayName','矩形障碍物');
        else
            patch(xr,yr,[0.25 0.25 0.25], ...
                'EdgeColor','k','LineWidth',1.2, ...
                'HandleVisibility','off');
        end
    end

    axis equal;
    grid on;
    xlabel('X / m');
    ylabel('Y / m');
    title('任务2：半圆弧双矩形风险通道避障');
    legend('Location','best');
    hold off;
    
    figure;
    plot(x_plot,y_plot,'LineWidth',1.4,'DisplayName','半圆弧全局参考');
    hold on;
    plot(history.x,history.y,'LineWidth',1.3,'DisplayName','MPC实际轨迹');

    for j = 1:Nobs
        cx = obs_rect(j,1);
        cy = obs_rect(j,2);
        w = obs_rect(j,3);
        h = obs_rect(j,4);

        xr = [cx-w/2,cx+w/2,cx+w/2,cx-w/2];
        yr = [cy-h/2,cy-h/2,cy+h/2,cy+h/2];

        if j == 1
            patch(xr,yr,[0.25 0.25 0.25], ...
                'EdgeColor','k','LineWidth',1.2, ...
                'DisplayName','矩形障碍物');
        else
            patch(xr,yr,[0.25 0.25 0.25], ...
                'EdgeColor','k','LineWidth',1.2, ...
                'HandleVisibility','off');
        end
    end

    axis equal;
    grid on;
    xlabel('X / m');
    ylabel('Y / m');
    title('任务2：半圆弧双矩形风险通道避障');
    legend('Location','best');
    hold off;

    figure;
    plot(history.t,history.mode,'LineWidth',1.2);
    grid on;
    xlabel('Time / s');
    ylabel('Mode');
    yticks([1 2 3]);
    yticklabels({'NORMAL','AVOID','RECOVER'});
    title('状态机');

    figure;
    plot(history.t,history.omega,'LineWidth',1.2);
    grid on;
    xlabel('Time / s');
    ylabel('\omega / (rad/s)');
    title('实际角速度');

    figure;
    plot(history.t,history.clearance,'LineWidth',1.2);
    hold on;
    yline(safe_margin,'--','安全裕量');
    yline(0,'--','碰撞边界');
    grid on;
    xlabel('Time / s');
    ylabel('矩形净距离 / m');
    title('车辆到矩形障碍物的最小净距离');
end

%% ========================================================================
%%                            局部函数
%% ========================================================================

function [xr,yr,thetar,vr,wr,s] = semicircle_reference( ...
    tvec,center,R,phi0,phi_end,V)

    tvec = tvec(:);

    s_max = R*(phi_end-phi0);
    s_raw = V*tvec;
    s = min(max(s_raw,0),s_max);

    phi = phi0 + s/R;

    xr = center(1) + R*cos(phi);
    yr = center(2) + R*sin(phi);

    thetar = phi + pi/2;

    vr = V*ones(size(tvec));
    wr = (V/R)*ones(size(tvec));

    finished = s_raw >= s_max;
    vr(finished) = 0;
    wr(finished) = 0;
end

function cost = function_avoid_curve( ...
    delta_target,x0,y0,theta0,v0,plan_Np,plan_Nc,T, ...
    Xref,Yref,THETAref,OMEGAref,Sref, ...
    Q,R,S_obst,Q_terminal,Q_terminal_theta, ...
    obs_rect,vehicle_radius,d_influence,obs_risk_weight, ...
    corridor_e,Q_corridor,corridor_s_min,corridor_s_max)

    [x_pre,y_pre,theta_pre,~] = rollout_plan_curve( ...
        delta_target,x0,y0,theta0,v0,plan_Np,plan_Nc,T,OMEGAref);

    e_long = zeros(plan_Np,1);
    e_lat  = zeros(plan_Np,1);

    J_obst = zeros(plan_Np,1);
    J_corridor = zeros(plan_Np,1);

    k_gate = 8;

    for i = 1:plan_Np
        tangent = [cos(THETAref(i)),sin(THETAref(i))];
        normal  = [-sin(THETAref(i)),cos(THETAref(i))];

        e_vec = [x_pre(i)-Xref(i),y_pre(i)-Yref(i)];

        e_long(i) = dot(e_vec,tangent);
        e_lat(i)  = dot(e_vec,normal);

        p = [x_pre(i),y_pre(i)];

        for j = 1:size(obs_rect,1)
            d_surface = point_rect_signed_distance(p,obs_rect(j,:));
            clearance = d_surface-vehicle_radius;

            if clearance < d_influence
                z = max((d_influence-clearance)/d_influence,0);
                J_obst(i) = J_obst(i) + ...
                    obs_risk_weight(j)*(z^2 + 0.3*z^4);
            end

            if clearance < 0
                J_obst(i) = J_obst(i) + ...
                    250*(abs(clearance)+0.02)^2;
            end
        end

        % 沿“弧长 s”而不是世界坐标 x 打开/关闭 corridor 代价
        w_in = 1/(1+exp(-k_gate*(Sref(i)-corridor_s_min)));
        w_out = 1/(1+exp( k_gate*(Sref(i)-corridor_s_max)));
        w_corridor = w_in*w_out;

        J_corridor(i) = ...
            w_corridor*(e_lat(i)-corridor_e)^2;
    end

    % AVOID期间仍然保持沿曲线向前的进度，
    % 但降低对原始横向 e_lat=0 的惩罚，让 corridor_e 主导避障横移。
    J_track_long = e_long'*Q*e_long;
    J_track_lat  = 0.12*(e_lat'*Q*e_lat);

    theta_terminal_error = ...
        wrap_pi(theta_pre(end)-THETAref(end));

    terminal_error = ...
        e_long(end)^2 + ...
        2.0*(e_lat(end)-corridor_e)^2;

    % delta_target 是相对名义圆弧角速度的 δω，
    % 因此这里只惩罚“额外转向”，不会惩罚正常圆弧本身需要的 V/R。
    cost = ...
        J_track_long + ...
        J_track_lat + ...
        delta_target'*R*delta_target + ...
        S_obst*sum(J_obst) + ...
        Q_corridor*sum(J_corridor) + ...
        Q_terminal*terminal_error + ...
        Q_terminal_theta*theta_terminal_error^2;
end

function [c,ceq] = rect_clearance_constraints_curve( ...
    delta_target,x0,y0,theta0,v0,plan_Np,plan_Nc,T, ...
    OMEGAref,THETAref,obs_rect,vehicle_radius,safe_margin,theta_limit)

    [x_pre,y_pre,theta_pre,~] = rollout_plan_curve( ...
        delta_target,x0,y0,theta0,v0,plan_Np,plan_Nc,T,OMEGAref);

    Nobs = size(obs_rect,1);
    required_clearance = vehicle_radius+safe_margin;

    c = zeros(plan_Np*Nobs + plan_Np,1);
    idx = 1;

    for i = 1:plan_Np
        for j = 1:Nobs
            d_surface = point_rect_signed_distance( ...
                [x_pre(i),y_pre(i)],obs_rect(j,:));

            c(idx) = required_clearance-d_surface;
            idx = idx+1;
        end
    end

    % 航向偏差相对“曲线上每一个参考航向”判断，
    % 而不是像任务1那样相对固定 theta=0 判断。
    for i = 1:plan_Np
        heading_error = wrap_pi(theta_pre(i)-THETAref(i));
        c(idx) = abs(heading_error)-theta_limit;
        idx = idx+1;
    end

    ceq = [];
end

function [x_pre,y_pre,theta_pre,omega_pre] = rollout_plan_curve( ...
    delta_target,x0,y0,theta0,v0,plan_Np,plan_Nc,T,OMEGAref)

    x_pre = zeros(plan_Np,1);
    y_pre = zeros(plan_Np,1);
    theta_pre = zeros(plan_Np,1);
    omega_pre = zeros(plan_Np,1);

    for i = 1:plan_Np
        if i <= plan_Nc
            delta_omega = delta_target(i);
        else
            % 控制时域后，不是 omega=0，
            % 而是 δω=0，恢复全局半圆的名义曲率。
            delta_omega = 0;
        end

        omega_pre(i) = OMEGAref(i)+delta_omega;

        if i == 1
            theta_pre(i) = theta0 + omega_pre(i)*T;
            x_pre(i) = x0 + v0*cos(theta_pre(i))*T;
            y_pre(i) = y0 + v0*sin(theta_pre(i))*T;
        else
            theta_pre(i) = theta_pre(i-1) + omega_pre(i)*T;
            x_pre(i) = x_pre(i-1) + v0*cos(theta_pre(i))*T;
            y_pre(i) = y_pre(i-1) + v0*sin(theta_pre(i))*T;
        end
    end
end

function [x_pre,y_pre,theta_pre,omega_pre] = generate_recover_plan_curve( ...
    x0,y0,theta0,v0,plan_Np,T, ...
    Xref,Yref,THETAref,x_ref_now,y_ref_now,theta_ref_now)

    x_pre = zeros(plan_Np,1);
    y_pre = zeros(plan_Np,1);
    theta_pre = zeros(plan_Np,1);
    omega_pre = zeros(plan_Np,1);

    normal_now = [-sin(theta_ref_now),cos(theta_ref_now)];

    e0 = dot([x0-x_ref_now,y0-y_ref_now],normal_now);

    e_theta0 = wrap_pi(theta0-theta_ref_now);
    e_theta_for_slope = ...
        max(min(e_theta0,deg2rad(25)),-deg2rad(25));

    m0 = tan(e_theta_for_slope);

    v_prog = max(abs(v0),0.55);
    L_recover = max(0.65,0.85*v_prog*plan_Np*T);

    % 先在每一个未来圆弧参考点的法向上施加一个逐渐归零的 Frenet 偏移
    for i = 1:plan_Np
        s_long = v_prog*i*T;
        tau = min(max(s_long/L_recover,0),1);

        h00 = 2*tau^3 - 3*tau^2 + 1;
        h10 = tau^3 - 2*tau^2 + tau;

        e_lat = h00*e0 + h10*L_recover*m0;

        if tau >= 1
            e_lat = 0;
        end

        normal_i = [-sin(THETAref(i)),cos(THETAref(i))];

        p = [Xref(i),Yref(i)] + e_lat*normal_i;

        x_pre(i) = p(1);
        y_pre(i) = p(2);
    end

    % 再根据几何轨迹本身求航向，而不是直接假设 theta=theta_ref+atan(de/ds)
    prev_x = x0;
    prev_y = y0;
    prev_theta = theta0;

    for i = 1:plan_Np
        raw_theta = atan2(y_pre(i)-prev_y,x_pre(i)-prev_x);

        theta_pre(i) = ...
            prev_theta + wrap_pi(raw_theta-prev_theta);

        omega_pre(i) = ...
            wrap_pi(theta_pre(i)-prev_theta)/T;

        prev_x = x_pre(i);
        prev_y = y_pre(i);
        prev_theta = theta_pre(i);
    end
end

function [s,e,theta_path] = obstacle_frenet_info( ...
    rect,center,R,phi0,phi_end)

    p = rect(1:2);

    phi = atan2(p(2)-center(2),p(1)-center(1));
    phi = min(max(phi,phi0),phi_end);

    p_ref = center + R*[cos(phi),sin(phi)];

    theta_path = phi+pi/2;
    normal = [-sin(theta_path),cos(theta_path)];

    e = dot(p-p_ref,normal);
    s = R*(phi-phi0);
end

function [s,e,theta_path] = circle_frenet_projection( ...
    p,center,R,phi0,phi_end)

    phi = atan2(p(2)-center(2),p(1)-center(1));
    phi = min(max(phi,phi0),phi_end);

    p_ref = center + R*[cos(phi),sin(phi)];

    theta_path = phi+pi/2;
    normal = [-sin(theta_path),cos(theta_path)];

    e = dot(p-p_ref,normal);
    s = R*(phi-phi0);
end

function d = point_rect_signed_distance(p,rect)
    cx = rect(1);
    cy = rect(2);

    hx = rect(3)/2;
    hy = rect(4)/2;

    q = abs(p-[cx,cy])-[hx,hy];

    outside = norm(max(q,0));
    inside = min(max(q(1),q(2)),0);

    d = outside+inside;
end

function a = wrap_pi(a)
    a = mod(a+pi,2*pi)-pi;
end
