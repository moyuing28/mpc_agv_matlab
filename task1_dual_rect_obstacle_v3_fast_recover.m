%这是运用了时变预测矩阵的误差状态控制的mpc控制与静态避障仿真，涉及避障代价算法的改进,和状态机的改进,重写一个轨迹回归的状态机

%这个代码将会优化掉一些无用的错误思路，加强可读性

%参考轨迹为正方形轨迹

% 存在的问题：！！！！现在这个代码是完全没有做泛化的，只验证了其在水平道路上规避矩形障碍物的可行性！！！
% 我把一部分的泛化代码写在了task3_double_rect_vec_obs_avoid_v2里面。不过至少证明这个思路可行的。
% 后续可以应该可以借助AI延续这个思路
% 直接改进一个泛化版本出来
clear; close all; clc;

%% 1. 共有参数设置

L = 0.5;              % 轮距 (m)

dt = 0.01;             % 采样时间 (s)

T = 0.05;

T_max = 20;           % 最大模拟时长 (s)

V_target = 0.8;       % 期望巡航速度 (m/s)

% LookAhead_Dist = 0.2; % 路径点切换阈值 (m) (判断是否到达航点的距离)(mpc不需要）

%在正方形的轨迹里面还是需要的，但是新的mpc代码没有用到它，作者（就是我）用了0.02提升了采样精度

Np = 60;              %预测步长

Nc = 30;              %控制步长

plan_Np = 20;         % 规划器预测步长：20*0.05 = 1.0 s，兼顾避障效果和计算速度

plan_Nc = 10;         % 规划器控制步长：10*0.05 = 0.5 s

%增广矩阵参数

Nx = 3;               %输出的增广维度参数

Nu = 2;               %输入的增广维度列数

% 初始位姿 (与前面算法统一)

x0 = 0; y0 = 0; theta0 = 0;%改动了初始位置，因为在没进入参考轨迹的情况下mpc控制发挥不好

squa_len = 7;

%% mpc控制的设计

% 初始状态设定

x = x0; y = y0; theta = theta0;

%增量矩阵的填充部分

O = zeros(Nu,Nx);

E = eye(Nu);

% 参考轨迹设计

%正方形轨迹：

Waypoints = [0, 0;

             squa_len, 0;

             squa_len, squa_len;

             0, squa_len;

             0, 0];

theta_ref_value =[0;pi/2;pi;-pi/2];

omega_ref = 0;%正方形轨迹的参考角速度，总为0

v_ref = V_target;

%障碍物设计：任务1，两块轴对齐长方形障碍物
% obs_rect 每行格式：[中心x, 中心y, 长度x, 宽度y]
% 第一块更靠近全局参考线 y=0 -> 高碰撞风险
% 第二块更远 -> 低碰撞风险，但仍作为通道另一侧边界
obs_rect = [ ...
   3.0,  0.38, 1.00, 0.36; ...   % 高风险障碍物
    3.0, -0.85, 1.00, 0.36];      % 低风险障碍物

Nobs = size(obs_rect,1);           %获取障碍物数量
vehicle_radius = 0.16;             % 用圆近似AGV外形半径
safe_margin = 0.10;               % 额外安全距离
d_influence = 0.70;               % 障碍物软代价影响距离

% 第一条边是 y=0，计算两个矩形之间的安全通道中心
obs_upper_idx = find(obs_rect(:,2) > 0, 1, 'first');
obs_lower_idx = find(obs_rect(:,2) < 0, 1, 'first');
upper_lower_edge = obs_rect(obs_upper_idx,2) - obs_rect(obs_upper_idx,4)/2;
lower_upper_edge = obs_rect(obs_lower_idx,2) + obs_rect(obs_lower_idx,4)/2;
corridor_y = 0.5*(upper_lower_edge + lower_upper_edge);
Q_corridor = 500;                  % 避障时把局部轨迹拉向两障碍物之间

% 根据障碍物与全局参考线的距离自动计算风险权重
obs_risk_weight = zeros(Nobs,1);
for j = 1:Nobs
    hy = obs_rect(j,4)/2;
    d_path = max(abs(obs_rect(j,2)) - hy - vehicle_radius, 0);
    obs_risk_weight(j) = 1 + 4*exp(-d_path/0.25);
    %d_path一定是非负数，值越小风险权重越大
end

% 双阈值的模式切换

d_trigger = 0.5;   % 距矩形边界的触发距离
% d_trigger = 1.35;   % 距矩形边界的触发距离
pass_margin = 0.20; % 完全越过矩形后再进入recover

draw = true;%绘图标志符，决定代码是否绘图

arri = true;%判断AGV是否到达轨迹起点，可以忽略，当起点处没有参考轨迹时mpc控制的效果是不好的，直接让mpc从车辆起点开始

t_arri = 0;%AGV到达轨迹起点的时间点

the_end = false;%参考的理想小车是否完成运动的标志

current_idx = 2;%正方形轨迹标签

Row=10;%松弛因子（可用可不用）

plan_draw =true;

v = 0;

omega = 0;

omega_seed = zeros(plan_Nc,1);  % fmincon滚动热启动

history.x=[];

history.y=[];

history.theta=[];

history.v=[];

history.omega=[];

history.obs_exist=[];

history.d_obs_dot=[];

history.plan_x_test={};

history.plan_y_test={};

history.mode =[];

history_ref.x=[];

history_ref.y=[];

% 预测序列和标志位的存储

history_pre_ref.x=[];

history_pre_ref.y=[];

history_pre_ref.theta=[];

history_pre_ref.current_idx =[];

history_pre_ref.the_end =[];

time = [];

history.pre_test={};%验证预测序列的

history.t = [];%仿真运行的时间序列

kesi=zeros(Nx+Nu,1);

num = 1;%实际轨迹序列的存储标签

pre_num = 1;%预测轨迹序列的存储标签

test_num =1;%用于测试的序列存储标签

plan_num =1;%测试规划层用的

% Q=10 * eye(Nx*Np,Nx*Np);    

Q_single = diag([15, 10, 5]);

Q = kron(eye(Np), Q_single);

% R=0.5*eye(Nu*Nc);

R_single = diag([0.5, 0.05]);   % Δv 权重 0.5, Δω 权重 0.05

R = kron(eye(Nc), R_single);

% Q_plan = 100*eye(plan_Np,plan_Np);

Q_plan = diag(linspace(40,120,plan_Np));

R_plan = 12*eye(plan_Nc,plan_Nc);

Q_terminal_avoid = 350;%避障终端位置代价

Q_terminal_theta_avoid = 450;%避障终端航向代价

Q_lat = 500 * eye(plan_Np);

Q_theta = 80 * eye(plan_Np);

Q_terminal_recover = 800;

Q_terminal_theta = 120;

%S_avoid，避障函数的权重。5-50效果都还不错

%S_recover，轨迹回归函数的权重

S_avoid = 70;

S_recover = 0;

%状态机参数

nornal = 1;

avoid = 2;

recover = 3;

mode = nornal;
recover_side = 0;   % 进入RECOVER时车辆位于全局轨迹哪一侧，用于防止状态机卡死

%% 参考轨迹生成+求取状态误差量

for t = 0:dt:T_max
    v_plan = V_target;
    if mode == nornal

        pre_t = (Np+1)*dt;%控制的预测时长,这里可能要斟酌一下

    else

        pre_t = plan_Np*T;%规划的预测时长

    end

    t_end = t + pre_t;

%全局参考轨迹的生成与轨迹跟踪那里一样，可以直接看那里的注释

    for i = t:dt:t_end

        U = [v;omega];

        err1 = [x-Waypoints(1,1),y-Waypoints(1,2)];

        distance1 = norm(err1);

        %if arri == false这一段代码仅作保留，实际是误差状态mpc脱离参考轨迹无法快速回正

            if arri == false

                if distance1 >= 0.02

                    x_ref =Waypoints(1,1);

                    y_ref =Waypoints(1,2);

                    history_ref.x(num) = x_ref;

                    history_ref.y(num) = y_ref;

                    theta_ref = theta_ref_value(1);

                    num = num+1;

                    kesi(1)= x-x_ref;

                    kesi(2)= y-y_ref;

                    kesi(3)= theta-theta_ref;

                    kesi(4)= v-V_target;

                    kesi(5)= omega-omega_ref;

                    arri = false;

                else

                    arri = true;

                    t_arri = t;

                end

            elseif  arri == true

                 switch current_idx

                    case  2

                        % 参考位置设立

                        x_pre_ref = V_target*(i-t_arri)-squa_len*(current_idx-2);

                        y_pre_ref = Waypoints(current_idx,2);        

                        % 参考路线转变标志判断，计算目前参考点到转换点的距离来判断是否进行转换

                        ref_err = [x_pre_ref-Waypoints(current_idx,1),y_pre_ref-Waypoints(current_idx,2)];

                        ref_distance = norm(ref_err);

                        history_pre_ref.x(pre_num) = x_pre_ref;

                        history_pre_ref.y(pre_num) = y_pre_ref;

                        history_pre_ref.theta(pre_num) = theta_ref_value(current_idx-1);

                        if ref_distance<=0.02

                            current_idx = current_idx+1;

                        end

                        history_pre_ref.current_idx(pre_num) = current_idx;

                        pre_num = pre_num+1;

                    case 3

                        x_pre_ref = Waypoints(current_idx,1);

                        y_pre_ref = V_target*(i-t_arri)-squa_len*(current_idx-2);

                        ref_err = [x_pre_ref-Waypoints(current_idx,1),y_pre_ref-Waypoints(current_idx,2)];

                        ref_distance = norm(ref_err);

                        history_pre_ref.x(pre_num) = x_pre_ref;

                        history_pre_ref.y(pre_num) = y_pre_ref;

                        history_pre_ref.theta(pre_num) = theta_ref_value(current_idx-1);               

                        if ref_distance<=0.02

                            current_idx = current_idx+1;

                        end

                        history_pre_ref.current_idx(pre_num) = current_idx;

                        pre_num = pre_num+1;

                     case 4

                        x_pre_ref = squa_len-(V_target*(i-t_arri)-squa_len*(current_idx-2));

                        y_pre_ref = Waypoints(current_idx,2);         

                        ref_err = [x_pre_ref-Waypoints(current_idx,1),y_pre_ref-Waypoints(current_idx,2)];

                        ref_distance = norm(ref_err);

                        history_pre_ref.x(pre_num) = x_pre_ref;

                        history_pre_ref.y(pre_num) = y_pre_ref;

                        history_pre_ref.theta(pre_num) = theta_ref_value(current_idx-1); 

                        if ref_distance<=0.02

                            current_idx = current_idx+1;

                        end

                        history_pre_ref.current_idx(pre_num) = current_idx;

                        pre_num = pre_num+1;

                    otherwise%最后一段比较特殊，如果仿真时间较长，参考点在仿真结束之前要会一直停在结束点处

                        if the_end ==false

                            x_pre_ref = Waypoints(current_idx,1);

                            y_pre_ref = squa_len-abs(V_target*(i-t_arri)-squa_len*(current_idx-2));

                            ref_err = [x_pre_ref-Waypoints(current_idx,1),y_pre_ref-Waypoints(current_idx,2)];

                            ref_distance = norm(ref_err);

                            history_pre_ref.x(pre_num) = x_pre_ref;

                            history_pre_ref.y(pre_num) = y_pre_ref;

                            history_pre_ref.theta(pre_num) = theta_ref_value(current_idx-1);

                            if ref_distance<=0.02

                                the_end =true;

                            end

                            history_pre_ref.the_end(pre_num) = the_end;

                            history_pre_ref.current_idx(pre_num) = current_idx;

                            pre_num = pre_num+1;

                        else%假设参考的小车先到了，实际小车的参考点定位原点

                            x_pre_ref = Waypoints(current_idx,1);

                            y_pre_ref = Waypoints(current_idx,2);              

                            history_pre_ref.x(pre_num) = x_pre_ref;

                            history_pre_ref.y(pre_num) = y_pre_ref;

                            history_pre_ref.theta(pre_num) = theta_ref_value(current_idx-1);

                            history_pre_ref.current_idx(pre_num) = current_idx;

                            pre_num = pre_num+1;

                        end

                        the_end = history_pre_ref.the_end(1);%重置轨迹结束的索引

                end

            end

    end

pre_num = 1;%重置预测时域的索引

current_idx = history_pre_ref.current_idx(1);%重置轨迹变换的索引

%% ================== 规划层 + 最终控制参考序列 ==================
% 保留全局参考，用于状态机和绘图
x_ref = history_pre_ref.x(1);
y_ref = history_pre_ref.y(1);
theta_ref = history_pre_ref.theta(1);

if mode == nornal
    % NORMAL：下层MPC直接跟踪全局参考，必须保留 Np+1 个点供 d_k 使用
    ctrl_ref_x = history_pre_ref.x(1:Np+1)';
    ctrl_ref_y = history_pre_ref.y(1:Np+1)';
    ctrl_ref_theta = unwrap(history_pre_ref.theta(1:Np+1)');
    ctrl_ref_v = V_target*ones(Np+1,1);
    ctrl_ref_omega = zeros(Np+1,1);

else
    %% AVOID / RECOVER：生成上层局部参考
    % 两个模式都先把全局参考重采样到规划器周期 T。
    % AVOID 仍使用 fmincon；RECOVER 改用解析的平滑回归轨迹，
    % 这样能显著减少计算量，而且不会因为短预测时域导致一直回不到 NORMAL。

    t_old = (0:length(history_pre_ref.x)-1)'*dt;
    t_new = (1:plan_Np)'*T;
    Xref = interp1(t_old, history_pre_ref.x(:), t_new, 'linear');
    Yref = interp1(t_old, history_pre_ref.y(:), t_new, 'linear');
    THETAref = interp1(t_old, unwrap(history_pre_ref.theta(:)), t_new, 'linear');

    if mode == avoid
        %% -------------------- AVOID：fmincon局部避障 --------------------
        %角速度约束
        omega_max = 4;
        omega_min = -4;
        lb_plan = omega_min*ones(plan_Nc,1);
        ub_plan = omega_max*ones(plan_Nc,1);

        options_plan = optimoptions('fmincon', ...
            'Algorithm','sqp', ...
            'Display','none', ...
            'MaxIterations',45, ...
            'MaxFunctionEvaluations',3000);

        % 角速度变化率约束
        delta_omega_plan = 0.35;
        Domega = eye(plan_Nc);
        Domega(2:end,1:end-1) = Domega(2:end,1:end-1) - eye(plan_Nc-1);

        omega_prev_vec = [omega; zeros(plan_Nc-1,1)];
        Aineq_plan = [ Domega; -Domega];
        bineq_plan = [ delta_omega_plan*ones(plan_Nc,1) + omega_prev_vec; ...
                       delta_omega_plan*ones(plan_Nc,1) - omega_prev_vec];

        % 控制时域末端角速度回到0，且末端航向回到全局路径方向。
        theta_target = THETAref(min(plan_Nc,length(THETAref)));
        % Aeq_plan = [zeros(1,plan_Nc-1), 1; ...
        %             T*ones(1,plan_Nc)];
        % beq_plan = [0; wrapToPi(theta_target-theta)];

        Aeq_plan = [zeros(1,plan_Nc-1),1];
        beq_plan = 0;
        % 
        %  Aeq_plan = [];
        % beq_plan = [];

        theta_limit = deg2rad(35);%将35度转化为弧度
        nonlcon = @(target)rect_clearance_constraints(target,x,y,theta,v_plan, ...
            plan_Np,plan_Nc,T,obs_rect,vehicle_radius,safe_margin, ...
            THETAref(1),theta_limit);

        [omega_opt,~,exitflag_plan] = fmincon( ...
            @(target)function_avoid_rect(target,x,y,theta,v_plan,plan_Np,plan_Nc,T, ...
            Yref,Xref,THETAref,Q_plan,R_plan,S_avoid,Q_terminal_avoid, ...
            Q_terminal_theta_avoid,obs_rect,vehicle_radius,d_influence, ...
            obs_risk_weight,corridor_y,Q_corridor), ...
            omega_seed,Aineq_plan,bineq_plan,Aeq_plan,beq_plan, ...
            lb_plan,ub_plan,nonlcon,options_plan);

        if exitflag_plan <= 0 || any(~isfinite(omega_opt))
            warning('局部避障规划器在 t=%.2f s 未正常收敛，沿用上一拍控制序列。',t);
            omega_opt = omega_seed;
        end

        % 热启动
        omega_seed = [omega_opt(2:end); 0];
        omega_seed(end) = 0;

        [x_pre,y_pre,theta_pre,omega_pre] = rollout_plan( ...
            omega_opt,x,y,theta,v_plan,plan_Np,plan_Nc,T);

    else
        %% -------------------- RECOVER：解析平滑回归 --------------------
        % 任务1的全局轨迹在障碍物附近是直线，因此不需要继续每0.01s跑fmincon。
        % 在Frenet坐标中用三次Hermite曲线把横向误差和航向误差平滑收敛到0。
        [x_pre,y_pre,theta_pre,omega_pre] = generate_recover_plan( ...
            x,y,theta,v,plan_Np,T,Xref,Yref,THETAref);

        % RECOVER不需要fmincon热启动，避免旧的避障控制继续影响回归
        omega_seed(:) = 0;
    end

    history.plan_x_test{plan_num,1} = x_pre;
    history.plan_y_test{plan_num,1} = y_pre;
    plan_num = plan_num+1;

    %% 多速率插值：规划器 T -> 跟踪MPC dt
    t_plan = [0; (1:plan_Np)'*T];
    x_plan = [x; x_pre];
    y_plan = [y; y_pre];
    theta_plan = unwrap([theta; theta_pre]);

    % 下层MPC只需要0~Np*dt = 0.6 s，规划器现在有1.0 s，覆盖范围足够。
    t_ctrl = (0:Np)'*dt;

    ctrl_ref_x = interp1(t_plan,x_plan,t_ctrl,'linear');
    ctrl_ref_y = interp1(t_plan,y_plan,t_ctrl,'linear');
    ctrl_ref_theta = interp1(t_plan,theta_plan,t_ctrl,'linear');
    ctrl_ref_v = V_target*ones(Np+1,1);

    ctrl_ref_omega = zeros(Np+1,1);
    ctrl_ref_omega(1:Np) = diff(ctrl_ref_theta)/dt;
    ctrl_ref_omega(Np+1) = ctrl_ref_omega(Np);
end

%% ================== d_k 与参考输入变化项 ==================
% 完整误差模型：Xi(k+1)=Atilde*Xi+Btilde*DeltaU+g_k
% g_k = [d_k; -Delta(u_ref)]
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

    xr_next = ctrl_ref_x(i+1);
    yr_next = ctrl_ref_y(i+1);
    theta_r_next = ctrl_ref_theta(i+1);

    d_ref(:,i) = [ ...
        xr_model_next - xr_next;
        yr_model_next - yr_next;
        wrapToPi(theta_model_next - theta_r_next)];

    delta_u_ref(:,i) = [ ...
        ctrl_ref_v(i+1)-ctrl_ref_v(i);
        ctrl_ref_omega(i+1)-ctrl_ref_omega(i)];
end

%% 当前增广误差
kesi(1) = x - ctrl_ref_x(1);
kesi(2) = y - ctrl_ref_y(1);
kesi(3) = wrapToPi(theta - ctrl_ref_theta(1));
kesi(4) = v - ctrl_ref_v(1);
kesi(5) = omega - ctrl_ref_omega(1);

%为每个预测步生成时变 A_i, B_i

A_cell = cell(Np, 1);

B_cell = cell(Np, 1);

C_tilde =[ 1 0 0 0 0;

        0 1 0 0 0;

        0 0 1 0 0];

for  j = 1:3

A_cell{j} = [1, 0, -dt*V_target*sin(theta);

             0, 1,  dt*V_target*cos(theta);

             0, 0, 1];

B_cell{j} = [dt*cos(theta), 0;

             dt*sin(theta), 0;

             0,               dt];

end

for i = 4:Np
    theta_r_i = ctrl_ref_theta(i);

    A_cell{i} = [1, 0, -dt*V_target*sin(theta_r_i);
                 0, 1,  dt*V_target*cos(theta_r_i);
                 0, 0, 1];

    B_cell{i} = [dt*cos(theta_r_i), 0;
                 dt*sin(theta_r_i), 0;
                 0,               dt];
end

% 构造时变 PHI (Nx*Np x (Nx+Nu))

PHI_cell = cell(Np, 1);

for i = 1:Np

    A_prod = eye(Nx + Nu);        % 连乘结果初始化为单位阵

    for k = 1:i

        A_tilde_k = [A_cell{k}, B_cell{k}; O, E];

        A_prod = A_tilde_k * A_prod;

    end

    PHI_cell{i} = C_tilde * A_prod;

end

PHI = cell2mat(PHI_cell);

%构造时变 THETA (Nx*Np x Nu*Nc)

THETA_cell = cell(Np, Nc);

for i = 1:Np

    for j = 1:Nc

        if j <= i

            % THETA_{i,j} = C * A_i * A_{i-1} * ... * A_{j+1} * B_j

            A_prod = eye(Nx + Nu);

            for k = (j+1):i

                A_tilde_k = [A_cell{k}, B_cell{k}; O, E];

                A_prod = A_tilde_k * A_prod;

            end

            B_tilde_j = [B_cell{j}; E];

            THETA_cell{i,j} = C_tilde * A_prod * B_tilde_j;

        else

            THETA_cell{i,j} = zeros(Nx, Nu);

        end

    end

end

THETA = cell2mat(THETA_cell);

%% ================== 加入 d_k 的自由误差预测 ==================

% 当前增广误差状态

Xi_free = kesi;

% 保存Np步未来输出误差

error_free = zeros(Nx*Np, 1);

for i = 1:Np

    %% 当前第i步的增广A矩阵

    A_tilde_i = [ ...

        A_cell{i}, B_cell{i};

        O,         E];

    %% 当前第i步的参考变化项

    % 由于目前正方形：

    % v_ref = 常数

    % omega_ref = 0

    %

    % 所以参考输入没有变化：

    % delta_u_ref = 0

    %

    % 增广扰动就是：

    g_i = [ ...
        d_ref(:,i);
        -delta_u_ref(:,i)];

    %% 在假设未来 ΔU = 0 的情况下传播误差

    Xi_free =A_tilde_i*Xi_free + g_i;

    %% 取出位置、航向三个误差

    error_free((i-1)*Nx+1:i*Nx) =  C_tilde*Xi_free;

end

%% 目标函数

% % error=PHI*kesi;...

H_cell=cell(2,2);

H_cell{1,1}=THETA'*Q*THETA+R;

H_cell{1,2}=zeros(Nu*Nc,1);

H_cell{2,1}=zeros(1,Nu*Nc);

H_cell{2,2}=Row;

H=cell2mat(H_cell);

H=(H+H')/2;

error=error_free;

f_cell=cell(1,2);

f_cell{1,1} = (error'*Q*THETA);

f_cell{1,2} = 0;

f=cell2mat(f_cell);

f = f';

%% 约束

A_t=zeros(Nc,Nc);

    for p=1:1:Nc

        for q=1:1:Nc

            if q<=p 

                A_t(p,q)=1;

            else 

                A_t(p,q)=0;

            end

        end 

    end 

A_I=kron(A_t,eye(Nu));

Ut=kron(ones(Nc,1),U);%U即是输入，在误差状态矩阵中作为输入误差向量，在直接状态矩阵中作为输入向量

% umin=[-10.5;  -0.436];%[min_vel, min_steer]维数与控制变量的个数相同

% umax=[10.0;   0.436]; %[max_vel, max_steer],%0.436rad = 25deg

% delta_umin = [-0.5;  -0.082]; % 0.0082rad = 0.47deg

% delta_umax = [0.5;  0.082];

if mode ==recover

umin = [-1.2; -4.0];       % 放宽角速度下限，缩窄速度下限

umax = [1.2;  4.0];        % 放宽两个速度上限，缩窄速度上限

% umin = [-0.8; -10.0];       % 放宽角速度下限，缩窄速度下限

% umax = [0.8;  10.0];        % 放宽两个速度上限，缩窄速度上限

delta_umin = [-0.2; -0.7];  % 放宽增量下限

delta_umax = [0.2;  0.7];   % 放宽增量上限

elseif mode == avoid

% umin = [0.6; -4.0];       % 放宽角速度下限，缩窄速度下限
% 
% umax = [1.0;  4.0];        % 放宽两个速度上限，缩窄速度上限

umin = [-1.2; -4.0];       % 放宽角速度下限，缩窄速度下限

umax = [1.2;  4.0];        % 放宽两个速度上限，缩窄速度上限

% umin = [-0.8; -10.0];       % 放宽角速度下限，缩窄速度下限

% umax = [0.8;  10.0];        % 放宽两个速度上限，缩窄速度上限

delta_umin = [-0.2; -0.7];  % 放宽增量下限

delta_umax = [0.2;  0.7];   % 放宽增量上限 

else

umin = [-1.2; -4.0];       % 放宽角速度下限，缩窄速度下限

umax = [1.2;  4.0];        % 放宽两个速度上限，缩窄速度上限

% umin = [-0.8; -10.0];       % 放宽角速度下限，缩窄速度下限

% umax = [0.8;  10.0];        % 放宽两个速度上限，缩窄速度上限

delta_umin = [-0.2; -0.7];  % 放宽增量下限

delta_umax = [0.2;  0.7];   % 放宽增量上限 

end

Umin=kron(ones(Nc,1),umin);

Umax=kron(ones(Nc,1),umax);

A_cons_cell={A_I zeros(Nu*Nc, 1); -A_I zeros(Nu*Nc, 1)};

b_cons_cell={Umax-Ut;-Umin+Ut};

A_cons=cell2mat(A_cons_cell);%（求解方程）状态量不等式约束增益矩阵，转换为绝对值的取值范围

b_cons=cell2mat(b_cons_cell);%（求解方程）状态量不等式约束的取值

% 状态量约束

delta_Umin = kron(ones(Nc,1),delta_umin);

delta_Umax = kron(ones(Nc,1),delta_umax);

lb = [delta_Umin; 0];%（求解方程）状态量下界

ub = [delta_Umax; 10];%（求解方程）状态量上界

%% qp求解

options = optimset('Algorithm','interior-point-convex'); 

warning off all  % close the warnings during computation

[X, fval,exitflag]=quadprog(H, f, A_cons, b_cons,[], [],lb,ub,[],options);

%% 求得下一时间步的状态量 

d_v = X(1);

d_omega = X(2);

v = v+d_v;

omega = omega+d_omega;

%等价的：

% v = d_v + dd_v + V_target;

% omega = d_omega + dd_omega + omega_ref;

vr = v + omega * L / 2; vl = v - omega * L / 2;

v_act = (vr + vl) / 2; w_act = (vr - vl) / L;

x = x + v_act * cos(theta) * dt; y = y + v_act * sin(theta) * dt; 

theta = wrapToPi(theta + w_act * dt);

%% 数据存储

history_ref.x(num) = x_ref;

history_ref.y(num) = y_ref;

history.x(num) = x;

history.y(num) = y;

history.theta(num) = theta;

history.v(num)=v;

history.omega(num)=omega;

% 到两个矩形障碍物边界的净安全距离（已经扣除车辆半径）
rect_clearance = zeros(Nobs,1);
for j = 1:Nobs
    rect_clearance(j) = point_rect_signed_distance([x,y],obs_rect(j,:)) - vehicle_radius;
end
history.d_obs_dot(num) = min(rect_clearance);

% 全局路径的切向/法向误差，用于状态机判断回归
Tangent_global = [cos(theta_ref), sin(theta_ref)];
Normal_global  = [-sin(theta_ref), cos(theta_ref)];
e_vec = [x-x_ref, y-y_ref];
e_long = dot(e_vec,Tangent_global);
e_lat  = dot(e_vec,Normal_global);

%% 任务1状态机：NORMAL -> AVOID -> RECOVER
switch mode
    case nornal
        % 只考虑车辆前方的障碍物；只要任一矩形进入触发距离就开始规划
        should_avoid = false;
        for j = 1:Nobs
            cx = obs_rect(j,1); cy = obs_rect(j,2);
            hx = obs_rect(j,3)/2; hy = obs_rect(j,4)/2;
            obs_relative = [cx-x, cy-y];
            s_center = dot(Tangent_global,obs_relative);
            half_extent_tangent = abs(Tangent_global(1))*hx + abs(Tangent_global(2))*hy;
            obstacle_ahead = (s_center + half_extent_tangent) > 0;

            if obstacle_ahead && rect_clearance(j) <= d_trigger
                should_avoid = true;
            end
        end
        if should_avoid
            mode = avoid;

            % 给fmincon一个“S形进入通道”的可行初值，避免从全零直线初值掉进局部最优。
            turn_sign = sign(corridor_y-y);
            if turn_sign == 0
                turn_sign = -1;
            end
            omega_seed(:) = 0;
            seed_shape = [0.20;0.40;0.55;0.55;0.40;0.20;0; ...
                         -0.20;-0.40;-0.55;-0.55;-0.40;-0.20;0];
            nseed = min(length(seed_shape),plan_Nc-1);
            omega_seed(1:nseed) = turn_sign*seed_shape(1:nseed);
            omega_seed(end) = 0;
        end

    case avoid
        % 车辆中心必须完全越过两块矩形的前缘，再进入回归模式
        passed_all = true;
        for j = 1:Nobs
            cx = obs_rect(j,1); cy = obs_rect(j,2);
            hx = obs_rect(j,3)/2; hy = obs_rect(j,4)/2;
            vehicle_relative = [x-cx, y-cy];%这里的向量方向和上面的normal的判断向量刚好相反
            s_center = dot(Tangent_global,vehicle_relative);
            half_extent_tangent = abs(Tangent_global(1))*hx + abs(Tangent_global(2))*hy;
            if s_center < half_extent_tangent + pass_margin%车辆质点完全穿过矩形+满足裕度才会判断成功穿过障碍物
                passed_all = false;
            end
        end
        if passed_all
            mode = recover;

            % 记录进入RECOVER时车辆位于全局路径的哪一侧。
            % 后面即使轻微越过参考线，也可以正常退出RECOVER，避免状态机一直卡住。
            recover_side = sign(e_lat);
            if recover_side == 0
                recover_side = sign(corridor_y);
                if recover_side == 0
                    recover_side = -1;
                end
            end

            omega_seed(:) = 0;
        end

    case recover
        e_theta_recover = wrapToPi(theta-theta_ref);

        % 条件1：已经足够接近全局轨迹，并且车头基本摆正。
        recover_converged = ...
            abs(e_lat) <= 0.07 && ...
            abs(e_theta_recover) <= deg2rad(10);

        % 条件2：车辆已经从原来一侧回到/轻微越过参考线，
        % 且航向偏差不大。此时交给NORMAL下层MPC完成最后的贴线，
        % 不再要求omega也必须恰好小于某个很严的阈值。
        recover_crossed = ...
            recover_side ~= 0 && ...
            e_lat*recover_side <= 0 && ...
            abs(e_theta_recover) <= deg2rad(15);

        if recover_converged || recover_crossed
            mode = nornal;
            omega_seed(:) = 0;
            recover_side = 0;
        end
end

%状态机状态变化存储

history.mode(num)=mode;

% time(num) = t;

history.t(num) = t;

num = num+1;

end

pre_test = cell2mat(history.pre_test);

if draw == true

    figure;

    % 绘制第一组曲线（参考轨迹），线宽调细

    plot(history_ref.x, history_ref.y, '-o', 'LineWidth', 0.8, 'MarkerSize', 2, 'MarkerFaceColor', 'b', 'DisplayName', '参考轨迹');

    hold on;  

    % 绘制第二组曲线（新轨迹），线宽同样调细

    plot(history.x, history.y, '-s', 'LineWidth', 0.8, 'MarkerSize', 2, 'MarkerFaceColor', 'r', 'DisplayName', 'mpc轨迹1');

if plan_draw ==true

    % 绘制第三组曲线（规划轨迹），线宽同样调细

     for i = 1:1:length(history.plan_x_test)

    plot(history.plan_x_test{i,1}, history.plan_y_test{i,1}, '-s', 'LineWidth', 0.8,'MarkerSize', 2, 'MarkerFaceColor', 'g');

     end

end

    %% 绘制两块矩形障碍物（patch兼容旧版MATLAB）
    for j = 1:Nobs
        cx = obs_rect(j,1); cy = obs_rect(j,2);
        w = obs_rect(j,3);  h = obs_rect(j,4);
        x_rect = [cx-w/2,cx+w/2,cx+w/2,cx-w/2];
        y_rect = [cy-h/2,cy-h/2,cy+h/2,cy+h/2];
        if j == 1
            patch(x_rect,y_rect,[0.25 0.25 0.25], ...
                'EdgeColor','k','LineWidth',1.2,'DisplayName','矩形障碍物');
        else
            patch(x_rect,y_rect,[0.25 0.25 0.25], ...
                'EdgeColor','k','LineWidth',1.2,'HandleVisibility','off');
        end
    end

    grid on;

    xlabel('X 轴');

    ylabel('Y 轴');

    title('坐标点可视化');
    axis equal;

    legend('Location', 'best');

    hold off;

end

if draw == true

    figure;

    % 绘制第一组曲线（参考轨迹），线宽调细

    plot(history_ref.x, history_ref.y, '-o', 'LineWidth', 0.8, 'MarkerSize', 2, 'Color', [0.00, 0.45, 0.74], 'DisplayName', '参考轨迹');

    hold on;  

    % 绘制第三组曲线（新轨迹），线宽同样调细

    plot(history.x, history.y, '-s', 'LineWidth', 0.8, 'MarkerSize', 2, 'Color', [0.85, 0.33, 0.10], 'DisplayName', 'mpc轨迹');

    %% 绘制两块矩形障碍物（patch兼容旧版MATLAB）
    for j = 1:Nobs
        cx = obs_rect(j,1); cy = obs_rect(j,2);
        w = obs_rect(j,3);  h = obs_rect(j,4);
        x_rect = [cx-w/2,cx+w/2,cx+w/2,cx-w/2];
        y_rect = [cy-h/2,cy-h/2,cy+h/2,cy+h/2];
        if j == 1
            patch(x_rect,y_rect,[0.25 0.25 0.25], ...
                'EdgeColor','k','LineWidth',1.2,'DisplayName','矩形障碍物');
        else
            patch(x_rect,y_rect,[0.25 0.25 0.25], ...
                'EdgeColor','k','LineWidth',1.2,'HandleVisibility','off');
        end
    end

    grid on;

    xlabel('X 轴');

    ylabel('Y 轴');

    title('坐标点可视化');
    axis equal;

    legend('Location', 'best');

    hold off;

end

% pid和mpc的v，omega的对比

figure;

plot(history.t, history.omega, 'b','DisplayName', 'mpc omega'); hold on;

grid on;

xlabel('Time (s)');

ylabel('\omega (rad/s)');

title('\omega 对比');

legend('MPC', 'PID');

figure;

plot(history.t, history.v, 'b','DisplayName', 'mpc omega'); hold on;

grid on;

xlabel('Time (s)');

ylabel('v (m/s)');

title('v 对比');

legend('MPC', 'PID');

% 绘制参考位置x，y的变化曲线。

figure;

plot(history.t, history_ref.x, 'b','DisplayName', '参考x轴'); hold on;

grid on;

plot(history.t, history_ref.y, 'r','DisplayName', '参考y轴'); hold on;

grid on;

plot(history.t, history.x, 'c','DisplayName', '实际x轴'); hold on;

grid on;

plot(history.t, history.y, 'y','DisplayName', '实际y轴'); hold on;

plot(history.t, history.mode,'DisplayName', '状态机数值'); hold on;

grid on;

grid on;

xlabel('Time (s)');

ylabel('m');

title('参考位置');

%% ================== 任务1：矩形障碍物避障规划函数 ==================
function cost = function_avoid_rect(target,x0,y0,theta0,v0,plan_Np,plan_Nc,T, ...
    Yref,Xref,THETAref,Q,R,S,Q_terminal,Q_terminal_theta,obs_rect, ...
    vehicle_radius,d_influence,obs_risk_weight,corridor_y,Q_corridor)

[x_pre,y_pre,theta_pre,omega_pre] = rollout_plan( ...
    target,x0,y0,theta0,v0,plan_Np,plan_Nc,T);

X_error = x_pre-Xref;
Y_error = y_pre-Yref;
J_obst = zeros(plan_Np,1);
J_corridor = zeros(plan_Np,1);

% 障碍区前后各扩展一段，让轨迹平滑地进入/退出安全通道
% 将矩形的中心向x轴正负向移动了obs_rect(:,3)/2+0.65的距离
obs_x_min = min(obs_rect(:,1)-obs_rect(:,3)/2)-0.65;
obs_x_max = max(obs_rect(:,1)+obs_rect(:,3)/2)+0.65;
k_gate = 8;   % sigmoid门函数斜率

for i = 1:plan_Np
    p = [x_pre(i),y_pre(i)];

    for j = 1:size(obs_rect,1)
        d_surface = point_rect_signed_distance(p,obs_rect(j,:));
        clearance = d_surface-vehicle_radius;%车辆边缘到矩形边缘的距离

        % 障碍物软代价：距离越近代价越高
        if clearance < d_influence
            z = max((d_influence-clearance)/d_influence,0);
            J_obst(i) = J_obst(i) + obs_risk_weight(j)*(z^2 + 0.3*z^4);
        end

        % 进入车辆包络内部时额外强惩罚（正常情况下硬约束会先阻止）
        if clearance < 0
            J_obst(i) = J_obst(i) + 250*(abs(clearance)+0.02)^2;
        end
    end

    % 用平滑门函数把局部参考逐渐拉向两矩形之间，而不是突然开/关。
    w_in  = 1/(1+exp(-k_gate*(x_pre(i)-obs_x_min)));
    w_out = 1/(1+exp( k_gate*(x_pre(i)-obs_x_max)));
    w_corridor = w_in*w_out;
    J_corridor(i) = w_corridor*(y_pre(i)-corridor_y)^2;
end

% 避障阶段不应该强迫 y 始终贴着全局线，否则会和“进入通道”目标打架。
% x方向仍跟踪全局进度，y方向只保留较弱的全局跟踪，再由 corridor 项主导横向位置。
J_track_x = X_error'*Q*X_error;
J_track_y = 0.12*(Y_error'*Q*Y_error);

% 终端希望：继续向前、位于安全通道附近、并且车头重新与全局路径平行。
theta_terminal_error = wrapToPi(theta_pre(end)-THETAref(end));
terminal_error = (x_pre(end)-Xref(end))^2 + 2.0*(y_pre(end)-corridor_y)^2;

cost = J_track_x ...
     + J_track_y ...
     + omega_pre(1:plan_Nc)'*R*omega_pre(1:plan_Nc) ...
     + S*sum(J_obst) ...
     + Q_corridor*sum(J_corridor) ...
     + Q_terminal*terminal_error ...
     + Q_terminal_theta*theta_terminal_error^2;
end

%% 回归全局轨迹的规划函数
function cost = function_recover_simple(target,x0,y0,theta0,v0,plan_Np,plan_Nc,T, ...
    Yref,Xref,THETAref,Q_lat,Q_theta,R,Q_terminal,Q_terminal_theta)

[x_pre,y_pre,theta_pre,omega_pre] = rollout_plan(target,x0,y0,theta0,v0,plan_Np,plan_Nc,T);

X_error = x_pre-Xref;
Y_error = y_pre-Yref;
THETA_error = zeros(plan_Np,1);
e_lat_plan = zeros(plan_Np,1);

for i = 1:plan_Np
    THETA_error(i) = wrapToPi(theta_pre(i)-THETAref(i));
    e_lat_plan(i) = -sin(THETAref(i))*X_error(i) + cos(THETAref(i))*Y_error(i);
end

cost = e_lat_plan'*Q_lat*e_lat_plan ...
     + THETA_error'*Q_theta*THETA_error ...
     + omega_pre(1:plan_Nc)'*R*omega_pre(1:plan_Nc) ...
     + Q_terminal*e_lat_plan(end)^2 ...
     + Q_terminal_theta*THETA_error(end)^2;
end

%% RECOVER解析轨迹：Frenet坐标中的三次Hermite平滑回归
function [x_pre,y_pre,theta_pre,omega_pre] = generate_recover_plan( ...
    x0,y0,theta0,v0,plan_Np,T,Xref,Yref,THETAref)

x_pre = zeros(plan_Np,1);
y_pre = zeros(plan_Np,1);
theta_pre = zeros(plan_Np,1);
omega_pre = zeros(plan_Np,1);

% 取当前全局参考方向，任务1障碍物位于第一条直线段，因此这里近似恒定。
theta_path = THETAref(1);
tangent = [cos(theta_path), sin(theta_path)];
normal  = [-sin(theta_path), cos(theta_path)];

% 当前点相对全局参考线的横向误差
p_ref = [Xref(1),Yref(1)];
p0 = [x0,y0];
e0 = dot(p0-p_ref,normal);

% 当前航向相对路径方向的夹角，转成横向曲线的初始斜率 de/ds。
e_theta0 = wrapToPi(theta0-theta_path);
e_theta_for_slope = max(min(e_theta0,deg2rad(25)),-deg2rad(25));%回归斜率？约束在±25°

m0 = tan(e_theta_for_slope);

% 当前点投影到全局路径上的垂足
p_proj0 = p0 - e0*normal;

% 让横向误差在预测时域前约85%完成回归，尾部留一点直线段供下层MPC稳定。
v_prog = max(abs(v0),0.55);
L_recover = max(0.65,0.85*v_prog*plan_Np*T);

for i = 1:plan_Np
    s_long = v_prog*i*T;
    tau = min(max(s_long/L_recover,0),1);

    % 三次Hermite：e(0)=e0, e'(0)=m0, e(1)=0, e'(1)=0
    % 构造一个逐渐衰减的横向偏移量
    h00 =  2*tau^3 - 3*tau^2 + 1;
    h10 =      tau^3 - 2*tau^2 + tau;

    e_lat = h00*e0 + h10*L_recover*m0;

    % de/dtau，再换成 de/ds
    %求偏移曲线的斜率
    dh00 = 6*tau^2 - 6*tau;
    dh10 = 3*tau^2 - 4*tau + 1;
    de_ds = (dh00*e0 + dh10*L_recover*m0)/L_recover;

    if tau >= 1
        e_lat = 0;
        de_ds = 0;
    end

    p = p_proj0 + s_long*tangent + e_lat*normal;
    x_pre(i) = p(1);
    y_pre(i) = p(2);

    theta_pre(i) = theta_path + atan(de_ds);

    if i == 1
        omega_pre(i) = wrapToPi(theta_pre(i)-theta0)/T;
    else
        omega_pre(i) = wrapToPi(theta_pre(i)-theta_pre(i-1))/T;
    end
end

end

%% 非线性硬安全约束：车辆圆形包络不得进入“矩形+安全距离”
function [c,ceq] = rect_clearance_constraints(target,x0,y0,theta0,v0, ...
    plan_Np,plan_Nc,T,obs_rect,vehicle_radius,safe_margin,theta_path,theta_limit)

[x_pre,y_pre,theta_pre,~] = rollout_plan( ...
    target,x0,y0,theta0,v0,plan_Np,plan_Nc,T);%生成规划轨迹序列

Nobs = size(obs_rect,1);
required_clearance = vehicle_radius + safe_margin;

% 前 plan_Np*Nobs 项：障碍物安全距离
% 后 plan_Np 项：航向不能偏离全局直线太多，防止回环/掉头
c = zeros(plan_Np*Nobs + plan_Np,1);
idx = 1;

for i = 1:plan_Np
    for j = 1:Nobs
        d_surface = point_rect_signed_distance([x_pre(i),y_pre(i)],obs_rect(j,:));
        c(idx) = required_clearance-d_surface;  % c<=0
        idx = idx+1;
    end
end

for i = 1:plan_Np
    heading_error = wrapToPi(theta_pre(i)-theta_path);
    c(idx) = abs(heading_error)-theta_limit;    % |e_theta| <= theta_limit
    idx = idx+1;
end

ceq = [];
end

%% 规划器统一滚动预测模型
function [x_pre,y_pre,theta_pre,omega_pre] = rollout_plan(target,x0,y0,theta0,v0,plan_Np,plan_Nc,T)
x_pre = zeros(plan_Np,1);
y_pre = zeros(plan_Np,1);
theta_pre = zeros(plan_Np,1);
omega_pre = zeros(plan_Np,1);

for i = 1:plan_Np
    if i <= plan_Nc
        omega_pre(i) = target(i);
    else
        % 控制时域以后保持直行。不能长期保持最后一个非零角速度，
        % 否则预测尾部会自然形成圆弧甚至回环。
        omega_pre(i) = 0;
    end

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

%% 点到轴对齐矩形的有符号距离
% rect = [cx,cy,Lx,Ly]
% 外部返回正距离，矩形内部返回负距离
function d = point_rect_signed_distance(p,rect)
cx = rect(1); cy = rect(2);
hx = rect(3)/2; hy = rect(4)/2;
q = abs(p-[cx,cy])-[hx,hy];
outside = norm(max(q,0));
inside = min(max(q(1),q(2)),0);
d = outside+inside;
end

