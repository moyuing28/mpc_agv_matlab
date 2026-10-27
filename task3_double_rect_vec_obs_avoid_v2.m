% 双车联合避障基础版：保持原单车代码框架，显式复制1号车/2号车
% =========================================================================
% 当前版本的目的：
% 1. 1号车、2号车分别跟踪两条互相包含的正方形轨迹；
% 2. 两辆车都保留原来的 NORMAL -> AVOID -> RECOVER 矩形避障结构；
% 3. 轨迹跟踪MPC、fmincon规划、状态机在主循环中显式分成两套，
%    所有车辆相关变量都严格使用 _1 / _2 后缀，避免相互覆盖；
% 4. “车-车动态避碰”勉强实现，后续可以继续改进。加入更泛化的避障软约束
% 5. 矩形，车车的判断算法已经一般化，水平道路，竖直道路都可以实现模式切换。
%
% 车辆分配：
%   1号车：内层正方形 Waypoints_upper，第一条边 y=d_traj，靠近上方高风险障碍物
%   2号车：外层正方形 Waypoints_lowwer，第一条边 y=0
%   
%   目前依旧存在的问题：
%   1.回归模式依旧还没有泛化，目前只针对水平道路
%   2.对障碍物的检测过于刻意，代码上默许了控制器会盯着矩形障碍物以及临近车辆，
%   后面的改进思路是在仿真里面增加一个“传感器”用来区分动静态障碍物和风险障碍物，
%   获得风险障碍物的坐标进行存储便于车辆的规避，这个在实机上可能要部署一套针对实际情况的算法
%   3.为了算法的实现目前的车辆是简化为了圆形的，实际中可能会发生碰撞，
%   后续要么将模型也改成成长方形，要么增大安全裕度防止碰撞
% 需要 Optimization Toolbox：fmincon, quadprog

clear; close all; clc;

%% 1. 共有参数设置
L = 0.5;                  % 轮距 (m)
dt = 0.01;                % 下层控制采样时间 (s)
T = 0.05;                 % 上层规划采样时间 (s)
T_max = 10;               % 最大模拟时长 (s)
V_target = 0.8;           % 期望巡航速度 (m/s)

Np = 60;                  % MPC预测步长
Nc = 30;                  % MPC控制步长

plan_Np = 20;             % 规划器预测步长
plan_Nc = 10;             % 规划器控制步长

Nx = 3;
Nu = 2;

O = zeros(Nu,Nx);
E = eye(Nu);

%% 2. 两条正方形轨迹
squa_len = 7;             % 外轨边长
d_traj = 0.6;             % 两轨道间隔

Waypoints_lowwer = [ ...
    0,0; ...
    squa_len,0; ...
    squa_len,squa_len; ...
    0,squa_len; ...
    0,0];

Waypoints_upper = [ ...
    d_traj,d_traj; ...
    squa_len-d_traj,d_traj; ...
    squa_len-d_traj,squa_len-d_traj; ...
    d_traj,squa_len-d_traj; ...
    d_traj,d_traj];

% 1号车走内轨，2号车走外轨
Waypoints_1 = Waypoints_upper;
Waypoints_2 = Waypoints_lowwer;

% 两辆车必须使用各自真实的正方形边长
squa_len_1 = Waypoints_1(2,1)-Waypoints_1(1,1);   % 5.8 m
squa_len_2 = Waypoints_2(2,1)-Waypoints_2(1,1);   % 7.0 m

theta_ref_value = [0;pi/2;pi;-pi/2];
omega_ref = 0;
v_ref = V_target;

%% 3. 障碍物
% 每行：[中心x, 中心y, x方向长度, y方向宽度]
obs_rect = [ ...
    2.50,  0.98, 1.00, 0.36; ...   % 上方高风险障碍物
    2.50, -0.85, 1.00, 0.36];      % 下方低风险障碍物

Nobs = size(obs_rect,1);

vehicle_radius = 0.16;
% vehicle_safe_margin = 0.20;

vehicle_lateral_margin = 0.2;
vehicle_distance_margin = 0.2;
% vehicle_return_safe_margin  = 0.30;
vehicle_return_safe_margin  = 0.23;
safe_margin = 0.10;
d_influence = 0.70;

% 两矩形之间的安全通道中心
obs_upper_idx = find(obs_rect(:,2)>0,1,'first');
obs_lower_idx = find(obs_rect(:,2)<0,1,'first');

upper_lower_edge = ...
    obs_rect(obs_upper_idx,2)-obs_rect(obs_upper_idx,4)/2;

lower_upper_edge = ...
    obs_rect(obs_lower_idx,2)+obs_rect(obs_lower_idx,4)/2;

corridor_y = 0.5*(upper_lower_edge+lower_upper_edge);

Q_corridor = 500;


% 两辆车分别根据“自己的第一条参考线”计算静态障碍物风险
path_y_1 = Waypoints_1(1,2);
path_y_2 = Waypoints_2(1,2);

obs_risk_weight_1 = zeros(Nobs,1);
obs_risk_weight_2 = zeros(Nobs,1);

for j = 1:Nobs
    hy = obs_rect(j,4)/2;

    d_path_1 = max( ...
        abs(obs_rect(j,2)-path_y_1)-hy-vehicle_radius,0);

    d_path_2 = max( ...
        abs(obs_rect(j,2)-path_y_2)-hy-vehicle_radius,0);

    obs_risk_weight_1(j) = 1+4*exp(-d_path_1/0.25);
    obs_risk_weight_2(j) = 1+4*exp(-d_path_2/0.25);
end

% d_trigger = 1.35;
d_trigger = 0.8;
d_lateral_trigger = 0.4;
pass_margin = 0.20;

draw = true;
plan_draw = false;

%% 4. 两辆车初始状态
% 1号车必须从内轨起点开始，避免误差状态MPC一开始脱离参考线
x_1 = Waypoints_1(1,1);
y_1 = Waypoints_1(1,2);
theta_1 = 0;

% 2号车从外轨起点开始
x_2 = d_traj;
y_2 = Waypoints_2(1,2);
theta_2 = 0;

v_1 = 0;
omega_1 = 0;

v_2 = 0;
omega_2 = 0;

% 两辆车各自的参考轨迹进度
test_traj = false;%双轨迹测试标志量

arri_1 = true;
arri_2 = true;

t_arri_1 = 0;
t_arri_2 = 0;

the_end_1 = false;
the_end_2 = false;

current_idx_1 = 2;
current_idx_2 = 2;

% 界定两车中哪辆是避障车，哪辆是避开避障车的车
% 设立true为避障车
avoid_vec_1 = false;
avoid_vec_2 = false;
%% 5. 两辆车独立状态机
nornal_1 = 1;
avoid_1 = 2;
recover_1 = 3;
mode_1 = nornal_1;
recover_side_1 = 0;

nornal_2 = 1;
avoid_2 = 2;
recover_2 = 3;
mode_2 = nornal_2;
recover_side_2 = 0;

%% 6. 两辆车独立热启动
omega_seed_1 = zeros(plan_Nc,1);
omega_seed_2 = zeros(plan_Nc,1);

%% 7. 历史量
history_1.x = [];
history_1.y = [];
history_1.theta = [];
history_1.v = [];
history_1.omega = [];
history_1.d_obs_dot = [];
history_1.plan_x_test = {};
history_1.plan_y_test = {};
history_1.mode = [];
history_1.t = [];

history_2.x = [];
history_2.y = [];
history_2.theta = [];
history_2.v = [];
history_2.omega = [];
history_2.d_obs_dot = [];
history_2.plan_x_test = {};
history_2.plan_y_test = {};
history_2.mode = [];
history_2.t = [];

history_1_ref.x = [];
history_1_ref.y = [];

history_2_ref.x = [];
history_2_ref.y = [];

history_1_pre_ref.x = [];
history_1_pre_ref.y = [];
history_1_pre_ref.theta = [];
history_1_pre_ref.current_idx = [];
history_1_pre_ref.the_end = [];

history_2_pre_ref.x = [];
history_2_pre_ref.y = [];
history_2_pre_ref.theta = [];
history_2_pre_ref.current_idx = [];
history_2_pre_ref.the_end = [];

%状态机的一些关键标志位的存储，用于验证两车是否按照预定思路在进行模式转换
history_1.should_avoid = [];
history_1.should_avoid_vec = [];
history_1.avoid_vec = [];
history_1.lat_vec = [];
history_1.distance_vec = [];

history_2.should_avoid = [];
history_2.should_avoid_vec = [];
history_2.avoid_vec = [];
history_2.lat_vec = [];
history_2.distance_vec = [];


kesi_1 = zeros(Nx+Nu,1);
kesi_2 = zeros(Nx+Nu,1);

num_1 = 1;
num_2 = 1;

plan_num_1 = 1;
plan_num_2 = 1;

%% 8. MPC / 规划权重
Q_single = diag([15,10,5]);
Q = kron(eye(Np),Q_single);

R_single = diag([0.5,0.05]);
R = kron(eye(Nc),R_single);

Q_plan = diag(linspace(40,120,plan_Np));
R_plan = 12*eye(plan_Nc);

Q_terminal_avoid = 350;
Q_terminal_theta_avoid = 450;

Q_lat = 500*eye(plan_Np); 
Q_theta = 80*eye(plan_Np); 
Q_terminal_recover = 800; 
Q_terminal_theta = 120; 

S_avoid = 70;
S_recover = 0; 

Row = 10;

%% 9. 双车主循环
for t = 0:dt:T_max
    v_plan = V_target;
    %% =====================================================================
    %% 1号车：全局参考轨迹预测序列
    %% =====================================================================
    if mode_1 == nornal_1
        pre_t_1 = (Np+1)*dt;
    else
        pre_t_1 = plan_Np*T;
    end
    t_end_1 = t + pre_t_1;

    % 每一个控制时刻重新生成本车未来参考序列
    history_1_pre_ref.x = [];
    history_1_pre_ref.y = [];
    history_1_pre_ref.theta = [];
    history_1_pre_ref.current_idx = [];
    history_1_pre_ref.the_end = [];

    pre_num_1 = 1;
    current_idx_temp_1 = current_idx_1;
    the_end_temp_1 = the_end_1;

    for i_ref_1 = t:dt:t_end_1

        switch current_idx_temp_1

            case 2
                x_pre_ref_1 = Waypoints_1(current_idx_temp_1-1,1) ...
                    + V_target*(i_ref_1-t_arri_1) ...
                    - squa_len_1*(current_idx_temp_1-2);
                y_pre_ref_1 = Waypoints_1(current_idx_temp_1,2);

                ref_err_1 = [ ...
                    x_pre_ref_1-Waypoints_1(current_idx_temp_1,1), ...
                    y_pre_ref_1-Waypoints_1(current_idx_temp_1,2)];
                ref_distance_1 = norm(ref_err_1);

                if ref_distance_1 <= 0.02
                    current_idx_temp_1 = current_idx_temp_1+1;
                end

            case 3
                x_pre_ref_1 = Waypoints_1(current_idx_temp_1,1);
                y_pre_ref_1 = Waypoints_1(current_idx_temp_1-1,2) ...
                    + V_target*(i_ref_1-t_arri_1) ...
                    - squa_len_1*(current_idx_temp_1-2);

                ref_err_1 = [ ...
                    x_pre_ref_1-Waypoints_1(current_idx_temp_1,1), ...
                    y_pre_ref_1-Waypoints_1(current_idx_temp_1,2)];
                ref_distance_1 = norm(ref_err_1);

                if ref_distance_1 <= 0.02
                    current_idx_temp_1 = current_idx_temp_1+1;
                end

            case 4
                x_pre_ref_1 = Waypoints_1(current_idx_temp_1-1,1) ...
                    - (V_target*(i_ref_1-t_arri_1) ...
                    - squa_len_1*(current_idx_temp_1-2));
                y_pre_ref_1 = Waypoints_1(current_idx_temp_1,2);

                ref_err_1 = [ ...
                    x_pre_ref_1-Waypoints_1(current_idx_temp_1,1), ...
                    y_pre_ref_1-Waypoints_1(current_idx_temp_1,2)];
                ref_distance_1 = norm(ref_err_1);

                if ref_distance_1 <= 0.02
                    current_idx_temp_1 = current_idx_temp_1+1;
                end

            otherwise
                if ~the_end_temp_1
                    x_pre_ref_1 = Waypoints_1(current_idx_temp_1,1);
                    y_pre_ref_1 = Waypoints_1(current_idx_temp_1-1,2) ...
                        - abs(V_target*(i_ref_1-t_arri_1) ...
                        - squa_len_1*(current_idx_temp_1-2));

                    ref_err_1 = [ ...
                        x_pre_ref_1-Waypoints_1(current_idx_temp_1,1), ...
                        y_pre_ref_1-Waypoints_1(current_idx_temp_1,2)];
                    ref_distance_1 = norm(ref_err_1);

                    if ref_distance_1 <= 0.02
                        the_end_temp_1 = true;
                    end
                else
                    x_pre_ref_1 = Waypoints_1(end,1);
                    y_pre_ref_1 = Waypoints_1(end,2);
                end
        end

        % 参考航向由当前所在边决定
        idx_theta_1 = min(max(current_idx_temp_1-1,1),4);
        theta_pre_ref_1 = theta_ref_value(idx_theta_1);

        history_1_pre_ref.x(pre_num_1) = x_pre_ref_1;
        history_1_pre_ref.y(pre_num_1) = y_pre_ref_1;
        history_1_pre_ref.theta(pre_num_1) = theta_pre_ref_1;
        history_1_pre_ref.current_idx(pre_num_1) = current_idx_temp_1;
        history_1_pre_ref.the_end(pre_num_1) = the_end_temp_1;

        pre_num_1 = pre_num_1+1;
    end

    % 只把当前时刻（预测序列第1点）对应的段号带到下一拍，
    % 不让未来预测过程提前改变真正的 current_idx。
    current_idx_1 = history_1_pre_ref.current_idx(1);
    the_end_1 = history_1_pre_ref.the_end(1);

    x_ref_1 = history_1_pre_ref.x(1);
    y_ref_1 = history_1_pre_ref.y(1);
    theta_ref_1 = history_1_pre_ref.theta(1);
     % 两车之间的安全通道中心

    corridor_vec_1 = 0.5*(upper_lower_edge+y_2)+vehicle_radius;
    % corridor_vec_2 = 0.5*(lower_upper_edge+y_1)-2*vehicle_radius;
    corridor_vec_2 = 0.5*(lower_upper_edge+y_1)-2*vehicle_radius;
    Q_corridor_vec = 500;
    % 
    % corridor_vec_1 = 0.5*(upper_lower_edge+y_2);
    % corridor_vec_2 = 0.5*(lower_upper_edge+y_1);
    % Q_corridor_vec = 500;

    if test_traj == true
    history_1_ref.x(num_1) = x_ref_1;
    history_1_ref.y(num_1) = y_ref_1;
    num_1 = num_1+1;
    else
    %% =====================================================================
    %% 1号车：规划层 + 最终控制参考序列
    %% =====================================================================
    if mode_1 == nornal_1

        % NORMAL：直接跟踪本车全局参考轨迹
        ctrl_ref_x_1 = history_1_pre_ref.x(1:Np+1)';
        ctrl_ref_y_1 = history_1_pre_ref.y(1:Np+1)';
        ctrl_ref_theta_1 = unwrap(history_1_pre_ref.theta(1:Np+1)');

        ctrl_ref_v_1 = V_target*ones(Np+1,1);
        ctrl_ref_omega_1 = zeros(Np+1,1);

    else

        % AVOID / RECOVER：先把 dt 参考重采样到规划器周期T
        t_old_1 = (0:length(history_1_pre_ref.x)-1)'*dt;
        t_new_1 = (1:plan_Np)'*T;

        Xref_1 = interp1(t_old_1,history_1_pre_ref.x(:),t_new_1,'linear');
        Yref_1 = interp1(t_old_1,history_1_pre_ref.y(:),t_new_1,'linear');
        THETAref_1 = interp1(t_old_1,unwrap(history_1_pre_ref.theta(:)),t_new_1,'linear');

        if mode_1 == avoid_1
            %% -------------------- 1号车 AVOID --------------------
            omega_max_1 = 4;
            omega_min_1 = -4;

            lb_plan_1 = omega_min_1*ones(plan_Nc,1);
            ub_plan_1 = omega_max_1*ones(plan_Nc,1);

            options_plan_1 = optimoptions('fmincon', ...
                'Algorithm','sqp', ...
                'Display','none', ...
                'MaxIterations',45, ...
                'MaxFunctionEvaluations',3000);

            % 角速度变化量约束
            delta_omega_plan_1 = 0.35;

            Domega_1 = eye(plan_Nc);
            Domega_1(2:end,1:end-1) = ...
                Domega_1(2:end,1:end-1)-eye(plan_Nc-1);

            omega_prev_vec_1 = [omega_1;zeros(plan_Nc-1,1)];

            Aineq_plan_1 = [Domega_1;-Domega_1];

            bineq_plan_1 = [ ...
                delta_omega_plan_1*ones(plan_Nc,1)+omega_prev_vec_1; ...
                delta_omega_plan_1*ones(plan_Nc,1)-omega_prev_vec_1];

            % 末端角速度回到0，末端航向回到本车参考路径方向
            theta_target_1 = THETAref_1(min(plan_Nc,length(THETAref_1)));

            Aeq_plan_1 = [zeros(1,plan_Nc-1),1];
            beq_plan_1 = 0;

            theta_limit_1 = deg2rad(35);

            nonlcon_1 = @(target)rect_clearance_constraints( ...
                target,x_1,y_1,theta_1,v_plan, ...
                plan_Np,plan_Nc,T,obs_rect,vehicle_radius,safe_margin, ...
                THETAref_1(1),theta_limit_1);

            [omega_opt_1,~,exitflag_plan_1] = fmincon( ...
                @(target)function_avoid_rect( ...
                target,x_1,y_1,theta_1,v_plan, ...
                plan_Np,plan_Nc,T, ...
                Yref_1,Xref_1,THETAref_1, ...
                Q_plan,R_plan,S_avoid,Q_terminal_avoid, ...
                Q_terminal_theta_avoid, ...
                obs_rect,vehicle_radius,d_influence, ...
                obs_risk_weight_1,corridor_vec_1,Q_corridor_vec), ...
                omega_seed_1, ...
                Aineq_plan_1,bineq_plan_1, ...
                Aeq_plan_1,beq_plan_1, ...
                lb_plan_1,ub_plan_1, ...
                nonlcon_1,options_plan_1);%corridor_y_vec

            if exitflag_plan_1 <= 0 || any(~isfinite(omega_opt_1))
                warning('1号车局部避障规划器在 t=%.2f s 未正常收敛，沿用上一拍控制序列。',t);
                omega_opt_1 = omega_seed_1;
            end

            % 滚动热启动
            omega_seed_1 = [omega_opt_1(2:end);0];
            omega_seed_1(end) = 0;

            [x_pre_1,y_pre_1,theta_pre_1,omega_pre_1] = rollout_plan( ...
                omega_opt_1,x_1,y_1,theta_1,v_plan, ...
                plan_Np,plan_Nc,T);

        else
            %% -------------------- 1号车 RECOVER --------------------
            [x_pre_1,y_pre_1,theta_pre_1,omega_pre_1] = ...
                generate_recover_plan( ...
                x_1,y_1,theta_1,v_1, ...
                plan_Np,T,Xref_1,Yref_1,THETAref_1);

            omega_seed_1(:) = 0;
        end

        history_1.plan_x_test{plan_num_1,1} = x_pre_1;
        history_1.plan_y_test{plan_num_1,1} = y_pre_1;
        plan_num_1 = plan_num_1+1;

        % 规划器 T -> 跟踪 MPC dt
        t_plan_1 = [0;(1:plan_Np)'*T];

        x_plan_1 = [x_1;x_pre_1];
        y_plan_1 = [y_1;y_pre_1];
        theta_plan_1 = unwrap([theta_1;theta_pre_1]);

        t_ctrl_1 = (0:Np)'*dt;

        ctrl_ref_x_1 = interp1(t_plan_1,x_plan_1,t_ctrl_1,'linear');
        ctrl_ref_y_1 = interp1(t_plan_1,y_plan_1,t_ctrl_1,'linear');
        ctrl_ref_theta_1 = interp1(t_plan_1,theta_plan_1,t_ctrl_1,'linear');

        ctrl_ref_v_1 = V_target*ones(Np+1,1);

        ctrl_ref_omega_1 = zeros(Np+1,1);
        ctrl_ref_omega_1(1:Np) = diff(ctrl_ref_theta_1)/dt;
        ctrl_ref_omega_1(Np+1) = ctrl_ref_omega_1(Np);
    end

    %% =====================================================================
    %% 1号车：d_k + 时变误差状态 MPC
    %% =====================================================================
    U_1 = [v_1;omega_1];

    d_ref_1 = zeros(Nx,Np);
    delta_u_ref_1 = zeros(Nu,Np);

    for i_mpc_1 = 1:Np
        xr_i_1 = ctrl_ref_x_1(i_mpc_1);
        yr_i_1 = ctrl_ref_y_1(i_mpc_1);
        theta_r_i_1 = ctrl_ref_theta_1(i_mpc_1);

        vr_i_1 = ctrl_ref_v_1(i_mpc_1);
        wr_i_1 = ctrl_ref_omega_1(i_mpc_1);

        xr_model_next_1 = ...
            xr_i_1+dt*vr_i_1*cos(theta_r_i_1);
        yr_model_next_1 = ...
            yr_i_1+dt*vr_i_1*sin(theta_r_i_1);
        theta_model_next_1 = ...
            theta_r_i_1+dt*wr_i_1;

        xr_next_1 = ctrl_ref_x_1(i_mpc_1+1);
        yr_next_1 = ctrl_ref_y_1(i_mpc_1+1);
        theta_r_next_1 = ctrl_ref_theta_1(i_mpc_1+1);

        d_ref_1(:,i_mpc_1) = [ ...
            xr_model_next_1-xr_next_1; ...
            yr_model_next_1-yr_next_1; ...
            wrap_pi(theta_model_next_1-theta_r_next_1)];

        delta_u_ref_1(:,i_mpc_1) = [ ...
            ctrl_ref_v_1(i_mpc_1+1)-ctrl_ref_v_1(i_mpc_1); ...
            ctrl_ref_omega_1(i_mpc_1+1)-ctrl_ref_omega_1(i_mpc_1)];
    end

    % 当前增广误差
    kesi_1(1) = x_1-ctrl_ref_x_1(1);
    kesi_1(2) = y_1-ctrl_ref_y_1(1);
    kesi_1(3) = wrap_pi(theta_1-ctrl_ref_theta_1(1));
    kesi_1(4) = v_1-ctrl_ref_v_1(1);
    kesi_1(5) = omega_1-ctrl_ref_omega_1(1);

    % 时变 A_i, B_i
    A_cell_1 = cell(Np,1);
    B_cell_1 = cell(Np,1);

    C_tilde_1 = [ ...
        1 0 0 0 0; ...
        0 1 0 0 0; ...
        0 0 1 0 0];

    for j_mpc_1 = 1:3
        A_cell_1{j_mpc_1} = [ ...
            1,0,-dt*V_target*sin(theta_1); ...
            0,1, dt*V_target*cos(theta_1); ...
            0,0,1];

        B_cell_1{j_mpc_1} = [ ...
            dt*cos(theta_1),0; ...
            dt*sin(theta_1),0; ...
            0,dt];
    end

    for i_mpc_1 = 4:Np
        theta_r_i_1 = ctrl_ref_theta_1(i_mpc_1);

        A_cell_1{i_mpc_1} = [ ...
            1,0,-dt*V_target*sin(theta_r_i_1); ...
            0,1, dt*V_target*cos(theta_r_i_1); ...
            0,0,1];

        B_cell_1{i_mpc_1} = [ ...
            dt*cos(theta_r_i_1),0; ...
            dt*sin(theta_r_i_1),0; ...
            0,dt];
    end

    % PHI（保留原代码结构，当前目标函数主要使用 error_free）
    PHI_cell_1 = cell(Np,1);

    for i_mpc_1 = 1:Np
        A_prod_1 = eye(Nx+Nu);

        for k_mpc_1 = 1:i_mpc_1
            A_tilde_k_1 = [ ...
                A_cell_1{k_mpc_1},B_cell_1{k_mpc_1}; ...
                O,E];

            A_prod_1 = A_tilde_k_1*A_prod_1;
        end

        PHI_cell_1{i_mpc_1} = C_tilde_1*A_prod_1;
    end

    PHI_1 = cell2mat(PHI_cell_1);

    % THETA
    THETA_cell_1 = cell(Np,Nc);

    for i_mpc_1 = 1:Np
        for j_mpc_1 = 1:Nc

            if j_mpc_1 <= i_mpc_1

                A_prod_1 = eye(Nx+Nu);

                for k_mpc_1 = (j_mpc_1+1):i_mpc_1
                    A_tilde_k_1 = [ ...
                        A_cell_1{k_mpc_1},B_cell_1{k_mpc_1}; ...
                        O,E];

                    A_prod_1 = A_tilde_k_1*A_prod_1;
                end

                B_tilde_j_1 = [B_cell_1{j_mpc_1};E];

                THETA_cell_1{i_mpc_1,j_mpc_1} = ...
                    C_tilde_1*A_prod_1*B_tilde_j_1;

            else
                THETA_cell_1{i_mpc_1,j_mpc_1} = zeros(Nx,Nu);
            end
        end
    end

    THETA_1 = cell2mat(THETA_cell_1);

    % 含 d_k 的自由误差预测
    Xi_free_1 = kesi_1;
    error_free_1 = zeros(Nx*Np,1);

    for i_mpc_1 = 1:Np

        A_tilde_i_1 = [ ...
            A_cell_1{i_mpc_1},B_cell_1{i_mpc_1}; ...
            O,E];

        g_i_1 = [ ...
            d_ref_1(:,i_mpc_1); ...
            -delta_u_ref_1(:,i_mpc_1)];

        Xi_free_1 = A_tilde_i_1*Xi_free_1+g_i_1;

        error_free_1((i_mpc_1-1)*Nx+1:i_mpc_1*Nx) = ...
            C_tilde_1*Xi_free_1;
    end

    %% 1号车 QP目标函数
    H_cell_1 = cell(2,2);
    H_cell_1{1,1} = THETA_1'*Q*THETA_1+R;
    H_cell_1{1,2} = zeros(Nu*Nc,1);
    H_cell_1{2,1} = zeros(1,Nu*Nc);
    H_cell_1{2,2} = Row;

    H_1 = cell2mat(H_cell_1);
    H_1 = (H_1+H_1')/2;

    f_cell_1 = cell(1,2);
    f_cell_1{1,1} = error_free_1'*Q*THETA_1;
    f_cell_1{1,2} = 0;

    f_1 = cell2mat(f_cell_1)';

    %% 1号车控制约束
    A_t_1 = zeros(Nc,Nc);

    for p_1 = 1:Nc
        for q_1 = 1:Nc
            if q_1 <= p_1
                A_t_1(p_1,q_1) = 1;
            else
                A_t_1(p_1,q_1) = 0;
            end
        end
    end

    A_I_1 = kron(A_t_1,eye(Nu));
    Ut_1 = kron(ones(Nc,1),U_1);

    % 目前三个模式沿用相同输入约束，保留原来的模式判断结构
    if mode_1 == recover_1
        umin_1 = [-1.2;-4.0];
        umax_1 = [ 1.2; 4.0];
        delta_umin_1 = [-0.2;-0.7];
        delta_umax_1 = [ 0.2; 0.7];

    elseif mode_1 == avoid_1
        umin_1 = [-1.2;-4.0];
        umax_1 = [ 1.2; 4.0];
        delta_umin_1 = [-0.2;-0.7];
        delta_umax_1 = [ 0.2; 0.7];

    else
        umin_1 = [-1;-4.0];
        umax_1 = [ 1; 4.0];
        delta_umin_1 = [-0.2;-0.7];
        delta_umax_1 = [ 0.2; 0.7];
    end

    Umin_1 = kron(ones(Nc,1),umin_1);
    Umax_1 = kron(ones(Nc,1),umax_1);

    A_cons_cell_1 = { ...
        A_I_1,zeros(Nu*Nc,1); ...
        -A_I_1,zeros(Nu*Nc,1)};

    b_cons_cell_1 = { ...
        Umax_1-Ut_1; ...
        -Umin_1+Ut_1};

    A_cons_1 = cell2mat(A_cons_cell_1);
    b_cons_1 = cell2mat(b_cons_cell_1);

    delta_Umin_1 = kron(ones(Nc,1),delta_umin_1);
    delta_Umax_1 = kron(ones(Nc,1),delta_umax_1);

    lb_1 = [delta_Umin_1;0];
    ub_1 = [delta_Umax_1;10];

    %% 1号车 QP求解
    options_qp_1 = optimset('Algorithm','interior-point-convex','Display','off');

    [X_qp_1,~,exitflag_qp_1] = quadprog( ...
        H_1,f_1,A_cons_1,b_cons_1, ...
        [],[],lb_1,ub_1,[],options_qp_1);

    if isempty(X_qp_1) || exitflag_qp_1 <= 0 || any(~isfinite(X_qp_1))
        warning('1号车下层MPC在 t=%.2f s 求解失败，本拍保持上一控制。',t);
        d_v_1 = 0;
        d_omega_1 = 0;
    else
        d_v_1 = X_qp_1(1);
        d_omega_1 = X_qp_1(2);
    end

    v_1 = v_1+d_v_1;
    omega_1 = omega_1+d_omega_1;

    % 差速运动学
    vr_1 = v_1+omega_1*L/2;
    vl_1 = v_1-omega_1*L/2;

    v_act_1 = (vr_1+vl_1)/2;
    w_act_1 = (vr_1-vl_1)/L;

    x_1 = x_1+v_act_1*cos(theta_1)*dt;
    y_1 = y_1+v_act_1*sin(theta_1)*dt;
    theta_1 = wrap_pi(theta_1+w_act_1*dt);

 
    end
    %% ======================== 2号车 =============================

    %% =====================================================================
    %% 2号车：全局参考轨迹预测序列
    %% =====================================================================
    if mode_2 == nornal_2
        pre_t_2 = (Np+1)*dt;
    else
        pre_t_2 = plan_Np*T;
    end
    t_end_2 = t + pre_t_2;

    % 每一个控制时刻重新生成本车未来参考序列
    history_2_pre_ref.x = [];
    history_2_pre_ref.y = [];
    history_2_pre_ref.theta = [];
    history_2_pre_ref.current_idx = [];
    history_2_pre_ref.the_end = [];

    pre_num_2 = 1;
    current_idx_temp_2 = current_idx_2;
    the_end_temp_2 = the_end_2;

    for i_ref_2 = t:dt:t_end_2

        switch current_idx_temp_2

            case 2
                x_pre_ref_2 = d_traj+Waypoints_2(current_idx_temp_2-1,1) ...
                    + V_target*(i_ref_2-t_arri_2) ...
                    - squa_len_2*(current_idx_temp_2-2);
                y_pre_ref_2 = Waypoints_2(current_idx_temp_2,2);

                ref_err_2 = [ ...
                    x_pre_ref_2-Waypoints_2(current_idx_temp_2,1), ...
                    y_pre_ref_2-Waypoints_2(current_idx_temp_2,2)];
                ref_distance_2 = norm(ref_err_2);

                if ref_distance_2 <= 0.02
                    current_idx_temp_2 = current_idx_temp_2+1;
                end

            case 3
                x_pre_ref_2 = Waypoints_2(current_idx_temp_2,1);
                y_pre_ref_2 = d_traj+Waypoints_2(current_idx_temp_2-1,2) ...
                    + V_target*(i_ref_2-t_arri_2) ...
                    - squa_len_2*(current_idx_temp_2-2);

                ref_err_2 = [ ...
                    x_pre_ref_2-Waypoints_2(current_idx_temp_2,1), ...
                    y_pre_ref_2-Waypoints_2(current_idx_temp_2,2)];
                ref_distance_2 = norm(ref_err_2);

                if ref_distance_2 <= 0.02
                    current_idx_temp_2 = current_idx_temp_2+1;
                end

            case 4
                x_pre_ref_2 = Waypoints_2(current_idx_temp_2-1,1) ...
                    - (V_target*(i_ref_2-t_arri_2) ...
                    - squa_len_2*(current_idx_temp_2-2));
                y_pre_ref_2 = Waypoints_2(current_idx_temp_2,2);

                ref_err_2 = [ ...
                    x_pre_ref_2-Waypoints_2(current_idx_temp_2,1), ...
                    y_pre_ref_2-Waypoints_2(current_idx_temp_2,2)];
                ref_distance_2 = norm(ref_err_2);

                if ref_distance_2 <= 0.02
                    current_idx_temp_2 = current_idx_temp_2+1;
                end

            otherwise
                if ~the_end_temp_2
                    x_pre_ref_2 = Waypoints_2(current_idx_temp_2,1);
                    y_pre_ref_2 = Waypoints_2(current_idx_temp_2-1,2) ...
                        - abs(V_target*(i_ref_2-t_arri_2) ...
                        - squa_len_2*(current_idx_temp_2-2));

                    ref_err_2 = [ ...
                        x_pre_ref_2-Waypoints_2(current_idx_temp_2,1), ...
                        y_pre_ref_2-Waypoints_2(current_idx_temp_2,2)];
                    ref_distance_2 = norm(ref_err_2);

                    if ref_distance_2 <= 0.02
                        the_end_temp_2 = true;
                    end
                else
                    x_pre_ref_2 = Waypoints_2(end,1);
                    y_pre_ref_2 = Waypoints_2(end,2);
                end
        end

        % 参考航向由当前所在边决定
        idx_theta_2 = min(max(current_idx_temp_2-1,1),4);
        theta_pre_ref_2 = theta_ref_value(idx_theta_2);

        history_2_pre_ref.x(pre_num_2) = x_pre_ref_2;
        history_2_pre_ref.y(pre_num_2) = y_pre_ref_2;
        history_2_pre_ref.theta(pre_num_2) = theta_pre_ref_2;
        history_2_pre_ref.current_idx(pre_num_2) = current_idx_temp_2;
        history_2_pre_ref.the_end(pre_num_2) = the_end_temp_2;

        pre_num_2 = pre_num_2+1;
    end

    % 只把当前时刻（预测序列第1点）对应的段号带到下一拍，
    % 不让未来预测过程提前改变真正的 current_idx。
    current_idx_2 = history_2_pre_ref.current_idx(1);
    the_end_2 = history_2_pre_ref.the_end(1);

    x_ref_2 = history_2_pre_ref.x(1);
    y_ref_2 = history_2_pre_ref.y(1);
    theta_ref_2 = history_2_pre_ref.theta(1);
    if test_traj == true

    history_2_ref.x(num_2) = x_ref_2;
    history_2_ref.y(num_2) = y_ref_2;
    num_2 = num_2+1;
    
    else
    %% =====================================================================
    %% 2号车：规划层 + 最终控制参考序列
    %% =====================================================================
    if mode_2 == nornal_2

        % NORMAL：直接跟踪本车全局参考轨迹
        ctrl_ref_x_2 = history_2_pre_ref.x(1:Np+1)';
        ctrl_ref_y_2 = history_2_pre_ref.y(1:Np+1)';
        ctrl_ref_theta_2 = unwrap(history_2_pre_ref.theta(1:Np+1)');

        ctrl_ref_v_2 = V_target*ones(Np+1,1);
        ctrl_ref_omega_2 = zeros(Np+1,1);

    else

        % AVOID / RECOVER：先把 dt 参考重采样到规划器周期 T
        t_old_2 = (0:length(history_2_pre_ref.x)-1)'*dt;
        t_new_2 = (1:plan_Np)'*T;

        Xref_2 = interp1(t_old_2,history_2_pre_ref.x(:),t_new_2,'linear');
        Yref_2 = interp1(t_old_2,history_2_pre_ref.y(:),t_new_2,'linear');
        THETAref_2 = interp1(t_old_2,unwrap(history_2_pre_ref.theta(:)),t_new_2,'linear');

        if mode_2 == avoid_2
            %% -------------------- 2号车 AVOID --------------------
            omega_max_2 = 4;
            omega_min_2 = -4;

            lb_plan_2 = omega_min_2*ones(plan_Nc,1);
            ub_plan_2 = omega_max_2*ones(plan_Nc,1);

            options_plan_2 = optimoptions('fmincon', ...
                'Algorithm','sqp', ...
                'Display','none', ...
                'MaxIterations',45, ...
                'MaxFunctionEvaluations',3000);

            % 角速度变化量约束
            delta_omega_plan_2 = 0.35;

            Domega_2 = eye(plan_Nc);
            Domega_2(2:end,1:end-1) = ...
                Domega_2(2:end,1:end-1)-eye(plan_Nc-1);

            omega_prev_vec_2 = [omega_2;zeros(plan_Nc-1,1)];

            Aineq_plan_2 = [Domega_2;-Domega_2];

            bineq_plan_2 = [ ...
                delta_omega_plan_2*ones(plan_Nc,1)+omega_prev_vec_2; ...
                delta_omega_plan_2*ones(plan_Nc,1)-omega_prev_vec_2];

            % 末端角速度回到0，末端航向回到本车参考路径方向
            theta_target_2 = THETAref_2(min(plan_Nc,length(THETAref_2)));

            Aeq_plan_2 = [zeros(1,plan_Nc-1),1];
            beq_plan_2 = 0;

            theta_limit_2 = deg2rad(35);

            nonlcon_2 = @(target)rect_clearance_constraints( ...
                target,x_2,y_2,theta_2,v_plan, ...
                plan_Np,plan_Nc,T,obs_rect,vehicle_radius,safe_margin, ...
                THETAref_2(1),theta_limit_2);

            [omega_opt_2,~,exitflag_plan_2] = fmincon( ...
                @(target)function_avoid_rect( ...
                target,x_2,y_2,theta_2,v_plan, ...
                plan_Np,plan_Nc,T, ...
                Yref_2,Xref_2,THETAref_2, ...
                Q_plan,R_plan,S_avoid,Q_terminal_avoid, ...
                Q_terminal_theta_avoid, ...
                obs_rect,vehicle_radius,d_influence, ...
                obs_risk_weight_2,corridor_vec_2,Q_corridor_vec), ...
                omega_seed_2, ...
                Aineq_plan_2,bineq_plan_2, ...
                Aeq_plan_2,beq_plan_2, ...
                lb_plan_2,ub_plan_2, ...
                nonlcon_2,options_plan_2);

            if exitflag_plan_2 <= 0 || any(~isfinite(omega_opt_2))
                warning('2号车局部避障规划器在 t=%.2f s 未正常收敛，沿用上一拍控制序列。',t);
                omega_opt_2 = omega_seed_2;
            end

            % 滚动热启动
            omega_seed_2 = [omega_opt_2(2:end);0];
            omega_seed_2(end) = 0;

            [x_pre_2,y_pre_2,theta_pre_2,omega_pre_2] = rollout_plan( ...
                omega_opt_2,x_2,y_2,theta_2,v_plan, ...
                plan_Np,plan_Nc,T);

        else
            %% -------------------- 2号车 RECOVER --------------------
            [x_pre_2,y_pre_2,theta_pre_2,omega_pre_2] = ...
                generate_recover_plan( ...
                x_2,y_2,theta_2,v_2, ...
                plan_Np,T,Xref_2,Yref_2,THETAref_2);

            omega_seed_2(:) = 0;
        end

        history_2.plan_x_test{plan_num_2,1} = x_pre_2;
        history_2.plan_y_test{plan_num_2,1} = y_pre_2;
        plan_num_2 = plan_num_2+1;

        % 规划器 T -> 跟踪 MPC dt
        t_plan_2 = [0;(1:plan_Np)'*T];

        x_plan_2 = [x_2;x_pre_2];
        y_plan_2 = [y_2;y_pre_2];
        theta_plan_2 = unwrap([theta_2;theta_pre_2]);

        t_ctrl_2 = (0:Np)'*dt;

        ctrl_ref_x_2 = interp1(t_plan_2,x_plan_2,t_ctrl_2,'linear');
        ctrl_ref_y_2 = interp1(t_plan_2,y_plan_2,t_ctrl_2,'linear');
        ctrl_ref_theta_2 = interp1(t_plan_2,theta_plan_2,t_ctrl_2,'linear');

        ctrl_ref_v_2 = V_target*ones(Np+1,1);

        ctrl_ref_omega_2 = zeros(Np+1,1);
        ctrl_ref_omega_2(1:Np) = diff(ctrl_ref_theta_2)/dt;
        ctrl_ref_omega_2(Np+1) = ctrl_ref_omega_2(Np);
    end

    %% =====================================================================
    %% 2号车：d_k + 时变误差状态 MPC
    %% =====================================================================
    U_2 = [v_2;omega_2];

    d_ref_2 = zeros(Nx,Np);
    delta_u_ref_2 = zeros(Nu,Np);

    for i_mpc_2 = 1:Np
        xr_i_2 = ctrl_ref_x_2(i_mpc_2);
        yr_i_2 = ctrl_ref_y_2(i_mpc_2);
        theta_r_i_2 = ctrl_ref_theta_2(i_mpc_2);

        vr_i_2 = ctrl_ref_v_2(i_mpc_2);
        wr_i_2 = ctrl_ref_omega_2(i_mpc_2);

        xr_model_next_2 = ...
            xr_i_2+dt*vr_i_2*cos(theta_r_i_2);
        yr_model_next_2 = ...
            yr_i_2+dt*vr_i_2*sin(theta_r_i_2);
        theta_model_next_2 = ...
            theta_r_i_2+dt*wr_i_2;

        xr_next_2 = ctrl_ref_x_2(i_mpc_2+1);
        yr_next_2 = ctrl_ref_y_2(i_mpc_2+1);
        theta_r_next_2 = ctrl_ref_theta_2(i_mpc_2+1);

        d_ref_2(:,i_mpc_2) = [ ...
            xr_model_next_2-xr_next_2; ...
            yr_model_next_2-yr_next_2; ...
            wrap_pi(theta_model_next_2-theta_r_next_2)];

        delta_u_ref_2(:,i_mpc_2) = [ ...
            ctrl_ref_v_2(i_mpc_2+1)-ctrl_ref_v_2(i_mpc_2); ...
            ctrl_ref_omega_2(i_mpc_2+1)-ctrl_ref_omega_2(i_mpc_2)];
    end

    % 当前增广误差
    kesi_2(1) = x_2-ctrl_ref_x_2(1);
    kesi_2(2) = y_2-ctrl_ref_y_2(1);
    kesi_2(3) = wrap_pi(theta_2-ctrl_ref_theta_2(1));
    kesi_2(4) = v_2-ctrl_ref_v_2(1);
    kesi_2(5) = omega_2-ctrl_ref_omega_2(1);

    % 时变 A_i, B_i
    A_cell_2 = cell(Np,1);
    B_cell_2 = cell(Np,1);

    C_tilde_2 = [ ...
        1 0 0 0 0; ...
        0 1 0 0 0; ...
        0 0 1 0 0];

    for j_mpc_2 = 1:3
        A_cell_2{j_mpc_2} = [ ...
            1,0,-dt*V_target*sin(theta_2); ...
            0,1, dt*V_target*cos(theta_2); ...
            0,0,1];

        B_cell_2{j_mpc_2} = [ ...
            dt*cos(theta_2),0; ...
            dt*sin(theta_2),0; ...
            0,dt];
    end

    for i_mpc_2 = 4:Np
        theta_r_i_2 = ctrl_ref_theta_2(i_mpc_2);

        A_cell_2{i_mpc_2} = [ ...
            1,0,-dt*V_target*sin(theta_r_i_2); ...
            0,1, dt*V_target*cos(theta_r_i_2); ...
            0,0,1];

        B_cell_2{i_mpc_2} = [ ...
            dt*cos(theta_r_i_2),0; ...
            dt*sin(theta_r_i_2),0; ...
            0,dt];
    end

    % PHI（保留原代码结构，当前目标函数主要使用 error_free）
    PHI_cell_2 = cell(Np,1);

    for i_mpc_2 = 1:Np
        A_prod_2 = eye(Nx+Nu);

        for k_mpc_2 = 1:i_mpc_2
            A_tilde_k_2 = [ ...
                A_cell_2{k_mpc_2},B_cell_2{k_mpc_2}; ...
                O,E];

            A_prod_2 = A_tilde_k_2*A_prod_2;
        end

        PHI_cell_2{i_mpc_2} = C_tilde_2*A_prod_2;
    end

    PHI_2 = cell2mat(PHI_cell_2);

    % THETA
    THETA_cell_2 = cell(Np,Nc);

    for i_mpc_2 = 1:Np
        for j_mpc_2 = 1:Nc

            if j_mpc_2 <= i_mpc_2

                A_prod_2 = eye(Nx+Nu);

                for k_mpc_2 = (j_mpc_2+1):i_mpc_2
                    A_tilde_k_2 = [ ...
                        A_cell_2{k_mpc_2},B_cell_2{k_mpc_2}; ...
                        O,E];

                    A_prod_2 = A_tilde_k_2*A_prod_2;
                end

                B_tilde_j_2 = [B_cell_2{j_mpc_2};E];

                THETA_cell_2{i_mpc_2,j_mpc_2} = ...
                    C_tilde_2*A_prod_2*B_tilde_j_2;

            else
                THETA_cell_2{i_mpc_2,j_mpc_2} = zeros(Nx,Nu);
            end
        end
    end

    THETA_2 = cell2mat(THETA_cell_2);

    % 含 d_k 的自由误差预测
    Xi_free_2 = kesi_2;
    error_free_2 = zeros(Nx*Np,1);

    for i_mpc_2 = 1:Np

        A_tilde_i_2 = [ ...
            A_cell_2{i_mpc_2},B_cell_2{i_mpc_2}; ...
            O,E];

        g_i_2 = [ ...
            d_ref_2(:,i_mpc_2); ...
            -delta_u_ref_2(:,i_mpc_2)];

        Xi_free_2 = A_tilde_i_2*Xi_free_2+g_i_2;

        error_free_2((i_mpc_2-1)*Nx+1:i_mpc_2*Nx) = ...
            C_tilde_2*Xi_free_2;
    end

    %% 2号车 QP目标函数
    H_cell_2 = cell(2,2);
    H_cell_2{1,1} = THETA_2'*Q*THETA_2+R;
    H_cell_2{1,2} = zeros(Nu*Nc,1);
    H_cell_2{2,1} = zeros(1,Nu*Nc);
    H_cell_2{2,2} = Row;

    H_2 = cell2mat(H_cell_2);
    H_2 = (H_2+H_2')/2;

    f_cell_2 = cell(1,2);
    f_cell_2{1,1} = error_free_2'*Q*THETA_2;
    f_cell_2{1,2} = 0;

    f_2 = cell2mat(f_cell_2)';

    %% 2号车控制约束
    A_t_2 = zeros(Nc,Nc);

    for p_2 = 1:Nc
        for q_2 = 1:Nc
            if q_2 <= p_2
                A_t_2(p_2,q_2) = 1;
            else
                A_t_2(p_2,q_2) = 0;
            end
        end
    end

    A_I_2 = kron(A_t_2,eye(Nu));
    Ut_2 = kron(ones(Nc,1),U_2);

    % 目前三个模式沿用相同输入约束，保留原来的模式判断结构
    if mode_2 == recover_2
        umin_2 = [-1.2;-4.0];
        umax_2 = [ 1.2; 4.0];
        delta_umin_2 = [-0.2;-0.7];
        delta_umax_2 = [ 0.2; 0.7];

    elseif mode_2 == avoid_2
        umin_2 = [-1.2;-4.0];
        umax_2 = [ 1.2; 4.0];
        delta_umin_2 = [-0.2;-0.7];
        delta_umax_2 = [ 0.2; 0.7];

    else
        umin_2 = [-1;-4.0];
        umax_2 = [ 1; 4.0];
        delta_umin_2 = [-0.2;-0.7];
        delta_umax_2 = [ 0.2; 0.7];
    end

    Umin_2 = kron(ones(Nc,1),umin_2);
    Umax_2 = kron(ones(Nc,1),umax_2);

    A_cons_cell_2 = { ...
        A_I_2,zeros(Nu*Nc,1); ...
        -A_I_2,zeros(Nu*Nc,1)};

    b_cons_cell_2 = { ...
        Umax_2-Ut_2; ...
        -Umin_2+Ut_2};

    A_cons_2 = cell2mat(A_cons_cell_2);
    b_cons_2 = cell2mat(b_cons_cell_2);

    delta_Umin_2 = kron(ones(Nc,1),delta_umin_2);
    delta_Umax_2 = kron(ones(Nc,1),delta_umax_2);

    lb_2 = [delta_Umin_2;0];
    ub_2 = [delta_Umax_2;10];

    %% 2号车 QP求解
    options_qp_2 = optimset('Algorithm','interior-point-convex','Display','off');

    [X_qp_2,~,exitflag_qp_2] = quadprog( ...
        H_2,f_2,A_cons_2,b_cons_2, ...
        [],[],lb_2,ub_2,[],options_qp_2);

    if isempty(X_qp_2) || exitflag_qp_2 <= 0 || any(~isfinite(X_qp_2))
        warning('2号车下层MPC在 t=%.2f s 求解失败，本拍保持上一控制。',t);
        d_v_2 = 0;
        d_omega_2 = 0;
    else
        d_v_2 = X_qp_2(1);
        d_omega_2 = X_qp_2(2);
    end

    v_2 = v_2+d_v_2;
    omega_2 = omega_2+d_omega_2;

    % 差速运动学
    vr_2 = v_2+omega_2*L/2;
    vl_2 = v_2-omega_2*L/2;

    v_act_2 = (vr_2+vl_2)/2;
    w_act_2 = (vr_2-vl_2)/L;

    x_2 = x_2+v_act_2*cos(theta_2)*dt;
    y_2 = y_2+v_act_2*sin(theta_2)*dt;
    theta_2 = wrap_pi(theta_2+w_act_2*dt);

    % 1车全局路径 Frenet 切向/法向
    Tangent_global_1 = [cos(theta_ref_1),sin(theta_ref_1)];
    Normal_global_1 = [-sin(theta_ref_1),cos(theta_ref_1)];

    e_vec_1 = [x_1-x_ref_1,y_1-y_ref_1];

    e_long_1 = dot(e_vec_1,Tangent_global_1); 
    e_lat_1 = dot(e_vec_1,Normal_global_1);

    % 2车全局路径 Frenet 切向/法向
    Tangent_global_2 = [cos(theta_ref_2),sin(theta_ref_2)];
    Normal_global_2 = [-sin(theta_ref_2),cos(theta_ref_2)];

    e_vec_2 = [x_2-x_ref_2,y_2-y_ref_2];

    e_long_2 = dot(e_vec_2,Tangent_global_2);
    e_lat_2 = dot(e_vec_2,Normal_global_2);



    vec_relative_1 = [ x_1-x_2, ...
                    y_1-y_2];

    vec_center_1 = ...
            dot(Tangent_global_1,vec_relative_1);

     d_vec_edge_crach_1_1 = sqrt((norm(vec_relative_1))^2 ...
        -(vec_center_1)^2)- 2*vehicle_radius;%车辆横向间隙

    d_vec_edge_crach_1_2 = abs((norm(vec_center_1))) ...
        - 2*vehicle_radius;%车辆距离间隙


    vec_relative_2 = [ x_1-x_2, ...
                    y_1-y_2];

    vec_center_2 = ...
        dot(Tangent_global_1,vec_relative_1);

    d_vec_edge_crach_2_1 = sqrt((norm(vec_relative_1))^2 ...
        -(vec_center_1)^2)- 2*vehicle_radius;%车辆横向间隙

    d_vec_edge_crach_2_2 = abs((norm(vec_center_1))) ...
        - 2*vehicle_radius;%车辆距离间隙
    
    history_1.lat_vec(num_1) = d_vec_edge_crach_1_1;
    history_2.lat_vec(num_2) = d_vec_edge_crach_2_1;

    history_1.distance_vec(num_1) = d_vec_edge_crach_1_2;
    history_2.distance_vec(num_2) = d_vec_edge_crach_2_2;
    %% =====================================================================
    %% 1号车：数据存储 + 矩形避障状态机
    %% =====================================================================
    history_1_ref.x(num_1) = x_ref_1;
    history_1_ref.y(num_1) = y_ref_1;

    history_1.x(num_1) = x_1;
    history_1.y(num_1) = y_1;
    history_1.theta(num_1) = theta_1;
    history_1.v(num_1) = v_1;
    history_1.omega(num_1) = omega_1;

    % 到矩形边界的净距离（扣除车辆半径）
    rect_clearance_1 = zeros(Nobs,1);

    for j_obs_1 = 1:Nobs
        rect_clearance_1(j_obs_1) = ...
            point_rect_signed_distance( ...
            [x_1,y_1],obs_rect(j_obs_1,:)) ...
            - vehicle_radius;
    end

    history_1.d_obs_dot(num_1) = min(rect_clearance_1);


    %% 1号车状态机：NORMAL -> AVOID -> RECOVER
    % =============================================================
    % 调试变量：每一个控制周期重新初始化
    % =============================================================
    should_avoid_1 = false;       % 本拍是否检测到静态障碍避让需求
    should_avoid_vec_1 = false;   % 本拍是否检测到车辆避让需求
    switch mode_1

        case nornal_1
            % should_avoid_1 = false;
            % should_avoid_vec_1 = false;
            for j_obs_1 = 1:Nobs

                cx_1 = obs_rect(j_obs_1,1);
                cy_1 = obs_rect(j_obs_1,2);

                hx_1 = obs_rect(j_obs_1,3)/2;
                hy_1 = obs_rect(j_obs_1,4)/2;
                
                h_1 = min(hx_1,hy_1);

                obs_relative_1 = [ ...
                    cx_1-x_1, ...
                    cy_1-y_1];

                s_center_1 = ...
                    dot(Tangent_global_1,obs_relative_1);
                 
                half_extent_tangent_1 = ...
                    abs(Tangent_global_1(1))*hx_1 ...
                    + abs(Tangent_global_1(2))*hy_1;
                %% ---------------- 1. 障碍物是否仍在前方 ----------------
                obstacle_ahead_1 = ...
                    (s_center_1+half_extent_tangent_1) > 0;
                %% ---------------- 2. 新增横向风险判据 ----------------
                % lateral_clearance_1 = ...
                %     abs(y_1-cy_1) - hy_1 - vehicle_radius;
                %这里假设避障总是发生在直线上的，并且矩形的长边一定是会发生碰撞的边缘。
                lateral_clearance_1 = sqrt((norm(obs_relative_1))^2-(s_center_1)^2)-h_1- vehicle_radius;
                
                if obstacle_ahead_1 && ...
                        rect_clearance_1(j_obs_1) <= d_trigger&& ... 
                        lateral_clearance_1 <= d_lateral_trigger

                    should_avoid_1 = true;
                    avoid_vec_1 = true;
                end
            end
            %% 这个判据还要改
             % d_vec_edge_crach = ...
             %        y_1 -y_2-2*vehicle_radius;
             % 
             % if d_vec_edge_crach <= vehicle_safe_margin
             %     should_avoid_vec_1 = true;
             % end
            %  on_first_segment = ...
            %         current_idx_1 == 2 && ...
            %         current_idx_2 == 2;
            %  if on_first_segment
            % 
            %     d_vec_edge_crach = ...
            %         abs(y_1-y_2) - 2*vehicle_radius;
            % 
            %     if d_vec_edge_crach <= vehicle_safe_margin
            %         should_avoid_vec_1 = true;
            %     end
            % end

            
            % vehicle_lateral_margin = 0.2;
            % vehicle_distance_margin = 0.6;

            if d_vec_edge_crach_1_1 <= vehicle_lateral_margin &&...
                d_vec_edge_crach_1_2 <= vehicle_distance_margin

                    should_avoid_vec_1 = true;
            end


            if should_avoid_1 || should_avoid_vec_1

                mode_1 = avoid_1;

                % 依据安全通道在本车哪一侧，生成S型角速度热启动
                turn_sign_1 = sign(corridor_vec_1-y_1);

                if turn_sign_1 == 0
                    turn_sign_1 = -1;
                end

                omega_seed_1(:) = 0;

                seed_shape_1 = [ ...
                    0.20;0.40;0.55;0.55;0.40;0.20;0; ...
                    -0.20;-0.40;-0.55;-0.55;-0.40;-0.20;0];

                nseed_1 = min(length(seed_shape_1),plan_Nc-1);

                omega_seed_1(1:nseed_1) = ...
                    turn_sign_1*seed_shape_1(1:nseed_1);

                omega_seed_1(end) = 0;
            end

        case avoid_1

            % 如果是避障车，必须完全越过所有矩形 + pass_margin 才进入RECOVER
            % 如果是要避开避障车的车辆，必须等避障车与其的距离足够大才才进入RECOVER
            if avoid_vec_1 ==true
                passed_all_1 = true;
    
                for j_obs_1 = 1:Nobs
    
                    cx_1 = obs_rect(j_obs_1,1);
                    cy_1 = obs_rect(j_obs_1,2);
    
                    hx_1 = obs_rect(j_obs_1,3)/2;
                    hy_1 = obs_rect(j_obs_1,4)/2;
    
                    vehicle_relative_1 = [ ...
                        x_1-cx_1, ...
                        y_1-cy_1];
    
                    s_center_1 = ...
                        dot(Tangent_global_1,vehicle_relative_1);
    
                    half_extent_tangent_1 = ...
                        abs(Tangent_global_1(1))*hx_1 ...
                        + abs(Tangent_global_1(2))*hy_1;
    
                    if s_center_1 < ...
                            half_extent_tangent_1+pass_margin
                        passed_all_1 = false;
                    end
                end
    
                if passed_all_1
                    
                    mode_1 = recover_1;
                    avoid_vec_1 = false;%避障车标志位重置
                    recover_side_1 = sign(e_lat_1);
    
                    if recover_side_1 == 0
                        % 注意：1号车自己的参考线不一定是 y=0，
                        % 因此必须使用 corridor_y-y_ref_1 判断回归侧。
                        recover_side_1 = sign(corridor_vec_1-y_ref_1);
    
                        if recover_side_1 == 0
                            recover_side_1 = -1;
                        end
                    end
    
                    omega_seed_1(:) = 0;
                end
            else               

            % vehicle_lateral_margin = 0.2;
            % vehicle_distance_margin = 0.6;

                if d_vec_edge_crach_1_1 >= vehicle_return_safe_margin 
                % d_vec_edge_crach = ...
                %         abs(y_1-y_2) - 2*vehicle_radius;
                % if  d_vec_edge_crach >= vehicle_return_safe_margin
                    mode_1 = recover_1;
                    should_avoid_vec_1 = false;
                    % 记录本次进入RECOVER时位于参考线哪一侧
                    recover_side_1 = sign(e_lat_1);

                    if recover_side_1 == 0

                        recover_side_1 = ...
                            sign(corridor_vec_1-y_ref_1);

                        if recover_side_1 == 0
                            recover_side_1 = -1;
                        end
                    end
                    omega_seed_1(:) = 0;
                end
            end

        case recover_1

            e_theta_recover_1 = ...
                wrap_pi(theta_1-theta_ref_1);
            %判断横向误差与角度误差是否满足条件
            recover_converged_1 = ...
                abs(e_lat_1) <= 0.07 && ...
                abs(e_theta_recover_1) <= deg2rad(10);
            %判断车辆是否穿过全局参考轨迹
            recover_crossed_1 = ...
                recover_side_1 ~= 0 && ...
                e_lat_1*recover_side_1 <= 0 && ...
                abs(e_theta_recover_1) <= deg2rad(15);

            if recover_converged_1 || recover_crossed_1

                mode_1 = nornal_1;
                omega_seed_1(:) = 0;
                recover_side_1 = 0;
            end
    end

    history_1.mode(num_1) = mode_1;
    history_1.t(num_1) = t;
    
    % =============================================================
    % 状态机调试量
    % =============================================================
    history_1.should_avoid(num_1) = double(should_avoid_1);
    history_1.should_avoid_vec(num_1) = double(should_avoid_vec_1);
    history_1.avoid_vec(num_1) = double(avoid_vec_1);
    
    num_1 = num_1+1;
    %% =====================================================================
    %% 2号车：数据存储 + 矩形避障状态机
    %% =====================================================================
    history_2_ref.x(num_2) = x_ref_2;
    history_2_ref.y(num_2) = y_ref_2;

    history_2.x(num_2) = x_2;
    history_2.y(num_2) = y_2;
    history_2.theta(num_2) = theta_2;
    history_2.v(num_2) = v_2;
    history_2.omega(num_2) = omega_2;

    % 到矩形边界的净距离（扣除车辆半径）
    rect_clearance_2 = zeros(Nobs,1);

    for j_obs_2 = 1:Nobs
        rect_clearance_2(j_obs_2) = ...
            point_rect_signed_distance( ...
            [x_2,y_2],obs_rect(j_obs_2,:)) ...
            - vehicle_radius;
    end

    history_2.d_obs_dot(num_2) = min(rect_clearance_2);


    %% 2号车状态机：NORMAL -> AVOID -> RECOVER
    % =============================================================
    % 调试变量：每一个控制周期重新初始化
    % =============================================================
    should_avoid_2 = false;       % 本拍是否检测到静态障碍避让需求
    should_avoid_vec_2 = false;   % 本拍是否检测到车辆避让需求
    switch mode_2

        case nornal_2
          
            % should_avoid_2 = false;
            % should_avoid_vec_2 = false;
            for j_obs_2 = 1:Nobs

                cx_2 = obs_rect(j_obs_2,1);
                cy_2 = obs_rect(j_obs_2,2);

                hx_2 = obs_rect(j_obs_2,3)/2;
                hy_2 = obs_rect(j_obs_2,4)/2;
                
                h_2 = min(hx_2,hy_2);
                obs_relative_2 = [ ...
                    cx_2-x_2, ...
                    cy_2-y_2];

                s_center_2 = ...
                    dot(Tangent_global_2,obs_relative_2);
                
                half_extent_tangent_2 = ...
                    abs(Tangent_global_2(1))*hx_2 ...
                    + abs(Tangent_global_2(2))*hy_2;
                %% ---------------- 1. 障碍物是否仍在前方 ----------------
                obstacle_ahead_2 = ...
                    (s_center_2+half_extent_tangent_2) > 0;
                  %% ---------------- 2. 新增横向风险判据 ----------------
                % lateral_clearance_1 = ...
                %     abs(y_1-cy_1) - hy_1 - vehicle_radius;
                %这里假设避障总是发生在直线上的，并且矩形的长边一定是会发生碰撞的边缘。
                lateral_clearance_2= sqrt((norm(obs_relative_2))^2-(s_center_2)^2)...
                -h_2- vehicle_radius;

                if obstacle_ahead_2 && ...
                        rect_clearance_2(j_obs_2) <= d_trigger&& ...
                        lateral_clearance_2 <= d_lateral_trigger

                    should_avoid_2 = true;
                    avoid_vec_2 = true;
                end
            end
            % d_vec_edge_crach = ...
            %         y_1 -y_2-2*vehicle_radius;
            % 
            %  if d_vec_edge_crach <= vehicle_safe_margin
            %      should_avoid_vec_2 = true;
            %  end

            % if on_first_segment
            %     d_vec_edge_crach = ...
            %         abs(y_1-y_2) - 2*vehicle_radius;
            % 
            %     if d_vec_edge_crach <= vehicle_safe_margin
            %         should_avoid_vec_2 = true;
            %     end
            % end

             % vehicle_lateral_margin = 0.2;
            % vehicle_distance_margin = 0.2;

            if d_vec_edge_crach_2_1 <= vehicle_lateral_margin &&...
                d_vec_edge_crach_2_2 <= vehicle_distance_margin

                    should_avoid_vec_2 = true;
            end


            if should_avoid_2 || should_avoid_vec_2

                mode_2 = avoid_2;
                % 依据安全通道在本车哪一侧，生成S型角速度热启动
                turn_sign_2 = sign(corridor_vec_2-y_2);

                if turn_sign_2 == 0
                    turn_sign_2 = -1;
                end

                omega_seed_2(:) = 0;

                seed_shape_2 = [ ...
                    0.20;0.40;0.55;0.55;0.40;0.20;0; ...
                    -0.20;-0.40;-0.55;-0.55;-0.40;-0.20;0];

                nseed_2 = min(length(seed_shape_2),plan_Nc-1);

                omega_seed_2(1:nseed_2) = ...
                    turn_sign_2*seed_shape_2(1:nseed_2);

                omega_seed_2(end) = 0;
            end

        case avoid_2
            % 如果是避障车，必须完全越过所有矩形 + pass_margin 才进入RECOVER
            % 如果是要避开避障车的车辆，必须等避障车与其的距离足够大才才进入RECOVER   
            if avoid_vec_2 == true
            passed_all_2 = true;

            for j_obs_2 = 1:Nobs

                cx_2 = obs_rect(j_obs_2,1);
                cy_2 = obs_rect(j_obs_2,2);

                hx_2 = obs_rect(j_obs_2,3)/2;
                hy_2 = obs_rect(j_obs_2,4)/2;

                vehicle_relative_2 = [ ...
                    x_2-cx_2, ...
                    y_2-cy_2];

                s_center_2 = ...
                    dot(Tangent_global_2,vehicle_relative_2);

                half_extent_tangent_2 = ...
                    abs(Tangent_global_2(1))*hx_2 ...
                    + abs(Tangent_global_2(2))*hy_2;
                %判断车辆是否已经出了障碍物
                if s_center_2 < ...
                        half_extent_tangent_2+pass_margin
                    passed_all_2 = false;
                end
            end

            if passed_all_2

                mode_2 = recover_2;
                avoid_vec_2 = false;%避障车标志位重置
                recover_side_2 = sign(e_lat_2);

                if recover_side_2 == 0
                    % 注意：2号车自己的参考线不一定是 y=0，
                    % 因此必须使用 corridor_y-y_ref_2 判断回归侧。
                    recover_side_2 = sign(corridor_vec_2-y_ref_2);

                    if recover_side_2 == 0
                        recover_side_2 = -1;
                    end
                end

                omega_seed_2(:) = 0;
            end
            else
                % d_vec_edge_crach = ...
                %         abs(y_1-y_2) - 2*vehicle_radius;
                
                % if  d_vec_edge_crach >= vehicle_return_safe_margin
                
            % vehicle_lateral_margin = 0.2;
            % vehicle_distance_margin = 0.6;

                if d_vec_edge_crach_2_1 >= vehicle_return_safe_margin 
                    mode_2 = recover_2;
                    should_avoid_vec_2 = false;
                    % 记录本次进入RECOVER时位于参考线哪一侧
                    recover_side_2 = sign(e_lat_2);

                    if recover_side_2 == 0

                        recover_side_2 = ...
                            sign(corridor_vec_2-y_ref_2);

                        if recover_side_2 == 0
                            recover_side_2 = -1;
                        end
                    end
                    omega_seed_2(:) = 0;
                end
            end
        case recover_2

            e_theta_recover_2 = ...
                wrap_pi(theta_2-theta_ref_2);
             %判断横向误差与角度误差是否满足条件
            recover_converged_2 = ...
                abs(e_lat_2) <= 0.07 && ...
                abs(e_theta_recover_2) <= deg2rad(10);
            %判断车辆是否穿过全局参考轨迹
            recover_crossed_2 = ...
                recover_side_2 ~= 0 && ...
                e_lat_2*recover_side_2 <= 0 && ...
                abs(e_theta_recover_2) <= deg2rad(15);

            if recover_converged_2 || recover_crossed_2

                mode_2 = nornal_2;
                omega_seed_2(:) = 0;
                recover_side_2 = 0;
            end
    end

    history_2.mode(num_2) = mode_2;
    history_2.t(num_2) = t;
    
    % =============================================================
    % 状态机调试量
    % =============================================================
    history_2.should_avoid(num_2) = double(should_avoid_2);
    history_2.should_avoid_vec(num_2) = double(should_avoid_vec_2);
    history_2.avoid_vec(num_2) = double(avoid_vec_2);
    
    num_2 = num_2+1;
    end
    %% 当前两车中心距离（先记录，暂时不参与规划）
    d_vehicle_12 = norm([x_1-x_2,y_1-y_2]);
    history_vehicle_distance(num_1-1) = d_vehicle_12; %#ok<SAGROW>

end

%% ========================================================================
%% 10. 绘图
%% ========================================================================
if draw
    figure;
    hold on;

    plot(history_1_ref.x,history_1_ref.y,'--', ...
        'LineWidth',1.0,'DisplayName','1号车参考轨迹');

    plot(history_1.x,history_1.y,'-', ...
        'LineWidth',1.2,'DisplayName','1号车实际轨迹');

    plot(history_2_ref.x,history_2_ref.y,'--', ...
        'LineWidth',1.0,'DisplayName','2号车参考轨迹');

    plot(history_2.x,history_2.y,'-', ...
        'LineWidth',1.2,'DisplayName','2号车实际轨迹');

    if plan_draw
        for i = 1:length(history_1.plan_x_test)
            plot(history_1.plan_x_test{i},history_1.plan_y_test{i}, ...
                'LineWidth',0.5,'HandleVisibility','off');
        end

        for i = 1:length(history_2.plan_x_test)
            plot(history_2.plan_x_test{i},history_2.plan_y_test{i}, ...
                'LineWidth',0.5,'HandleVisibility','off');
        end
    end

    for j = 1:Nobs
        cx = obs_rect(j,1);
        cy = obs_rect(j,2);
        w = obs_rect(j,3);
        h = obs_rect(j,4);

        x_rect = [cx-w/2,cx+w/2,cx+w/2,cx-w/2];
        y_rect = [cy-h/2,cy-h/2,cy+h/2,cy+h/2];

        if j == 1
            patch(x_rect,y_rect,[0.25 0.25 0.25], ...
                'EdgeColor','k','LineWidth',1.2, ...
                'DisplayName','矩形障碍物');
        else
            patch(x_rect,y_rect,[0.25 0.25 0.25], ...
                'EdgeColor','k','LineWidth',1.2, ...
                'HandleVisibility','off');
        end
    end

    grid on;
    axis equal;
    xlabel('X (m)');
    ylabel('Y (m)');
    title('双车轨迹跟踪 + 各自矩形避障');
    legend('Location','best');
    hold off;
end

figure;
plot(history_1.t,history_1.mode,'LineWidth',1.1); hold on;
plot(history_2.t,history_2.mode,'LineWidth',1.1);
grid on;
xlabel('Time (s)');
ylabel('mode');
title('双车状态机');
legend('1号车','2号车');

figure;
plot(history_1.t,history_1.v,'LineWidth',1.1); hold on;
plot(history_2.t,history_2.v,'LineWidth',1.1);
grid on;
xlabel('Time (s)');
ylabel('v (m/s)');
title('双车速度');
legend('1号车','2号车');

figure;
plot(history_1.t,history_1.lat_vec,'LineWidth',1.1); hold on;
plot(history_2.t,history_1.distance_vec,'LineWidth',1.1);
grid on;
xlabel('Time (s)');
ylabel(' m ');
title('双车距离');
legend('横向间隙','纵向间隙')

figure;
plot(history_1.t,history_1.omega,'LineWidth',1.1); hold on;
plot(history_2.t,history_2.omega,'LineWidth',1.1);
grid on;
xlabel('Time (s)');
ylabel('\omega (rad/s)');
title('双车角速度');
legend('1号车','2号车');

figure;
plot(history_1.t,history_vehicle_distance,'LineWidth',1.1);
grid on;
xlabel('Time (s)');
ylabel('d_{12} (m)');
title('双车中心距离（当前仅记录，不参与控制）');

%% ========================================================================
%% 双车状态机标志位绘制
%% ========================================================================
figure;

stairs(history_1.t, ...
       history_1.should_avoid, ...
       'LineWidth',1.3);
hold on;

stairs(history_1.t, ...
       history_1.should_avoid_vec, ...
       'LineWidth',1.3);

stairs(history_1.t, ...
       history_1.avoid_vec, ...
       'LineWidth',1.3);

grid on;

ylim([-0.1,1.1]);

yticks([0 1]);

xlabel('Time (s)');
ylabel('Flag');

title('1号车避障逻辑标志位');

legend( ...
    'should\_avoid\_1：静态障碍触发', ...
    'should\_avoid\_vec\_1：避车触发', ...
    'avoid\_vec\_1：当前AVOID原因', ...
    'Location','best');

figure;

stairs(history_2.t, ...
       history_2.should_avoid, ...
       'LineWidth',1.3);
hold on;

stairs(history_2.t, ...
       history_2.should_avoid_vec, ...
       'LineWidth',1.3);

stairs(history_2.t, ...
       history_2.avoid_vec, ...
       'LineWidth',1.3);

grid on;

ylim([-0.1,1.1]);

yticks([0 1]);

xlabel('Time (s)');
ylabel('Flag');

title('2号车避障逻辑标志位');

legend( ...
    'should\_avoid\_2：静态障碍触发', ...
    'should\_avoid\_vec\_2：避车触发', ...
    'avoid\_vec\_2：当前AVOID原因', ...
    'Location','best');

%% ========================================================================
%% 11. 任务1：矩形障碍物避障规划函数
%% ========================================================================
function cost = function_avoid_rect(target,x0,y0,theta0,v0, ...
    plan_Np,plan_Nc,T,Yref,Xref,THETAref,Q,R,S, ...
    Q_terminal,Q_terminal_theta,obs_rect,vehicle_radius, ...
    d_influence,obs_risk_weight,corridor_y,Q_corridor)

[x_pre,y_pre,theta_pre,omega_pre] = rollout_plan( ...
    target,x0,y0,theta0,v0,plan_Np,plan_Nc,T);

X_error = x_pre-Xref;
Y_error = y_pre-Yref;

J_obst = zeros(plan_Np,1);
J_corridor = zeros(plan_Np,1);
J_corridor_vec = zeros(plan_Np,1);

obs_x_min = min(obs_rect(:,1)-obs_rect(:,3)/2)-0.65;
obs_x_max = max(obs_rect(:,1)+obs_rect(:,3)/2)+0.65;

k_gate = 8;

for i = 1:plan_Np
    p = [x_pre(i),y_pre(i)];

    for j = 1:size(obs_rect,1)

        d_surface = ...
            point_rect_signed_distance(p,obs_rect(j,:));

        clearance = d_surface-vehicle_radius;

        if clearance < d_influence
            z = max((d_influence-clearance)/d_influence,0);

            J_obst(i) = J_obst(i) ...
                + obs_risk_weight(j)*(z^2+0.3*z^4);
        end

        if clearance < 0
            J_obst(i) = J_obst(i) ...
                + 250*(abs(clearance)+0.02)^2;
        end
    end

    w_in = 1/(1+exp(-k_gate*(x_pre(i)-obs_x_min)));
    w_out = 1/(1+exp(k_gate*(x_pre(i)-obs_x_max)));

    w_corridor = w_in*w_out;

    J_corridor(i) = ...
        w_corridor*(y_pre(i)-corridor_y)^2;

    % J_corridor_vec(i) = ...
    %     w_corridor*(y_pre(i)-corridor_vec_y)^2;
end

J_track_x = X_error'*Q*X_error;
J_track_y = 0.12*(Y_error'*Q*Y_error);

theta_terminal_error = ...
    wrap_pi(theta_pre(end)-THETAref(end));

terminal_error = ...
    (x_pre(end)-Xref(end))^2 ...
    + 2.0*(y_pre(end)-corridor_y)^2;

cost = J_track_x ...
    + J_track_y ...
    + omega_pre(1:plan_Nc)'*R*omega_pre(1:plan_Nc) ...
    + S*sum(J_obst) ...
    + Q_corridor*sum(J_corridor) ...
    + Q_terminal*terminal_error ...
    + Q_terminal_theta*theta_terminal_error^2;
% + Q_corridor_vec*sum(J_corridor_vec) ...;

end

%% ========================================================================
%% 12. RECOVER：Frenet三次Hermite解析回归
%% ========================================================================
function [x_pre,y_pre,theta_pre,omega_pre] = ...
    generate_recover_plan( ...
    x0,y0,theta0,v0,plan_Np,T,Xref,Yref,THETAref)

x_pre = zeros(plan_Np,1);
y_pre = zeros(plan_Np,1);
theta_pre = zeros(plan_Np,1);
omega_pre = zeros(plan_Np,1);

theta_path = THETAref(1);

tangent = [cos(theta_path),sin(theta_path)];
normal = [-sin(theta_path),cos(theta_path)];

p_ref = [Xref(1),Yref(1)];
p0 = [x0,y0];

e0 = dot(p0-p_ref,normal);

e_theta0 = wrap_pi(theta0-theta_path);

e_theta_for_slope = ...
    max(min(e_theta0,deg2rad(25)),-deg2rad(25));

m0 = tan(e_theta_for_slope);

p_proj0 = p0-e0*normal;

v_prog = max(abs(v0),0.55);

L_recover = max( ...
    0.65,0.85*v_prog*plan_Np*T);

for i = 1:plan_Np

    s_long = v_prog*i*T;

    tau = min(max(s_long/L_recover,0),1);

    h00 = 2*tau^3-3*tau^2+1;
    h10 = tau^3-2*tau^2+tau;

    e_lat = ...
        h00*e0+h10*L_recover*m0;

    dh00 = 6*tau^2-6*tau;
    dh10 = 3*tau^2-4*tau+1;

    de_ds = ...
        (dh00*e0+dh10*L_recover*m0)/L_recover;

    if tau >= 1
        e_lat = 0;
        de_ds = 0;
    end

    p = ...
        p_proj0+s_long*tangent+e_lat*normal;

    x_pre(i) = p(1);
    y_pre(i) = p(2);

    theta_pre(i) = ...
        theta_path+atan(de_ds);

    if i == 1
        omega_pre(i) = ...
            wrap_pi(theta_pre(i)-theta0)/T;
    else
        omega_pre(i) = ...
            wrap_pi(theta_pre(i)-theta_pre(i-1))/T;
    end
end

end

%% ========================================================================
%% 13. 非线性硬安全约束
%% ========================================================================
function [c,ceq] = rect_clearance_constraints( ...
    target,x0,y0,theta0,v0, ...
    plan_Np,plan_Nc,T,obs_rect, ...
    vehicle_radius,safe_margin,theta_path,theta_limit)

[x_pre,y_pre,theta_pre,~] = rollout_plan( ...
    target,x0,y0,theta0,v0,plan_Np,plan_Nc,T);

Nobs = size(obs_rect,1);

required_clearance = ...
    vehicle_radius+safe_margin;

c = zeros(plan_Np*Nobs+plan_Np,1);

idx = 1;

for i = 1:plan_Np
    for j = 1:Nobs

        d_surface = ...
            point_rect_signed_distance( ...
            [x_pre(i),y_pre(i)],obs_rect(j,:));

        c(idx) = ...
            required_clearance-d_surface;

        idx = idx+1;
    end
end

for i = 1:plan_Np

    heading_error = ...
        wrap_pi(theta_pre(i)-theta_path);

    c(idx) = ...
        abs(heading_error)-theta_limit;

    idx = idx+1;
end

ceq = [];

end

%% ========================================================================
%% 14. 规划器统一滚动预测模型
%% ========================================================================
function [x_pre,y_pre,theta_pre,omega_pre] = ...
    rollout_plan( ...
    target,x0,y0,theta0,v0,plan_Np,plan_Nc,T)

x_pre = zeros(plan_Np,1);
y_pre = zeros(plan_Np,1);
theta_pre = zeros(plan_Np,1);
omega_pre = zeros(plan_Np,1);

for i = 1:plan_Np

    if i <= plan_Nc
        omega_pre(i) = target(i);
    else
        omega_pre(i) = 0;
    end

    if i == 1

        theta_pre(i) = ...
            theta0+omega_pre(i)*T;

        x_pre(i) = ...
            x0+v0*cos(theta_pre(i))*T;

        y_pre(i) = ...
            y0+v0*sin(theta_pre(i))*T;

    else

        theta_pre(i) = ...
            theta_pre(i-1)+omega_pre(i)*T;

        x_pre(i) = ...
            x_pre(i-1)+v0*cos(theta_pre(i))*T;

        y_pre(i) = ...
            y_pre(i-1)+v0*sin(theta_pre(i))*T;
    end
end

end

%% ========================================================================
%% 15. 点到轴对齐矩形的有符号距离
%% ========================================================================
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

%% ========================================================================
%% 16. wrap到[-pi,pi]
%% ========================================================================
function a = wrap_pi(a)
a = mod(a+pi,2*pi)-pi;
end

%% ========================================================================
%% 17. 计算当前障碍物与车辆风险距离的仿真
%% ========================================================================
