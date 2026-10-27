%这是运用了时变预测矩阵的误差状态控制的mpc控制
%参考轨迹为圆形轨迹
clear; close all; clc;

%% 1. 共有参数设置
L = 0.5;              % 轮距 (m)
dt = 0.01;             % 采样时间 (s)
T_max = 30;           % 最大模拟时长 (s)
V_target = 0.8;       % 期望巡航速度 (m/s)
Np = 60;              %预测步长
Nc = 30;               %控制步长
%增广矩阵参数
Nx = 3;               %输出的增广维度参数
Nu = 2;               %输入的增广维度列数
% 初始位姿 (与前面算法统一)
x0 = 0; y0 = -0; theta0 = 0;%改动了初始位置，因为在没进入参考轨迹的情况下mpc控制发挥不好

% mpc控制的设计
% 初始状态设定
x = x0; y = y0; theta = theta0;
%增量矩阵的填充部分
O = zeros(Nu,Nx);
E = eye(Nu);
% 参考轨迹设计？

%% 圆形参考轨迹
R_circle = 3;             % 圆半径 (m)
circle_center = [0, R_circle];
% 初始位置 (0,0)，初始航向 theta=0
phi0 = -pi/2;
% 圆周运动角速度
omega_ref = V_target / R_circle;
v_ref = V_target;
draw = true;%绘图标志符，决定代码是否绘图
the_end = false;%参考的理想小车是否完成运动的标志
Row=10;%松弛因子（可用可不用）
pre_t = Np*dt;%控制的预测时长

v = 0;
omega = 0;

history.x=[];
history.y=[];
history.theta=[];
history.v=[];
history.omega=[];

history_ref.x=[];
history_ref.y=[];

% 预测序列和标志位的存储
history_pre_ref.x=[];
history_pre_ref.y=[];
history_pre_ref.theta=[];

history_pre_ref.the_end =[];

time = [];
% history.pre_test={};%验证预测序列的
history.t = [];%仿真运行的时间序列
kesi=zeros(Nx+Nu,1);
num = 1;
pre_num = 1;
% test_num =1;
%% 参考轨迹生成+求取状态误差量
for t = 0:dt:T_max
    for i = t:dt:(t+pre_t)
        U = [v;omega];
        %圆周运动角速度
        %当前仿真时间 t 对应的参考轨迹
        phi_ref = phi0 + omega_ref * i;

        x_pre_ref = circle_center(1) + R_circle*cos(phi_ref);
        y_pre_ref = circle_center(2) + R_circle*sin(phi_ref);
        theta_pre_ref = wrapToPi(phi_ref + pi/2);

        history_pre_ref.x(pre_num) = x_pre_ref;
        history_pre_ref.y(pre_num) = y_pre_ref;
        history_pre_ref.theta(pre_num) = theta_pre_ref;

        pre_num = pre_num+1;         
    end                                    
pre_num = 1;%重置预测时域的索引
%d_k的计算
for i =1:1:Np
    xr_i = history_pre_ref.x(i);
    yr_i = history_pre_ref.y(i);
    theta_r_i = history_pre_ref.theta(i);
    %基于模型预测的下一状态
    xr_model_next = ...
    xr_i + dt*v_ref*cos(theta_r_i);

    yr_model_next = ...
        yr_i + dt*v_ref*sin(theta_r_i);
    
    theta_model_next = ...
        theta_r_i + dt*omega_ref;

    xr_next = history_pre_ref.x(i+1);
    yr_next = history_pre_ref.y(i+1);
    theta_r_next = history_pre_ref.theta(i+1);
    
    d_theta =  wrapToPi(theta_model_next - theta_r_next);
    d_ref(:,i) = [ ...
    xr_model_next - xr_next;
    yr_model_next - yr_next;
    d_theta
    ];

 end
x_ref = history_pre_ref.x(1);
y_ref = history_pre_ref.y(1);
theta_ref = history_pre_ref.theta(1);
%% 误差向量的求取
%误差向量求取
kesi(1)= x-x_ref;
kesi(2)= y-y_ref;
kesi(3)= wrapToPi(theta-theta_ref);
% kesi(3)= (theta-theta_ref);
kesi(4)= v-V_target;
kesi(5)= omega-omega_ref;

% Q=10 * eye(Nx*Np,Nx*Np);    
Q_single = diag([15, 10, 5]);
Q = kron(eye(Np), Q_single);
R=0.5*eye(Nu*Nc);
% ---- 为每个预测步生成时变 A_i, B_i ----
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
    theta_r_i =history_pre_ref.theta(i);
    A_cell{i} = [1, 0, -dt*V_target*sin(theta_r_i);
                 0, 1,  dt*V_target*cos(theta_r_i);
                 0, 0, 1];
    B_cell{i} = [dt*cos(theta_r_i), 0;
                 dt*sin(theta_r_i), 0;
                 0,               dt];
end

% % for i = 1:Np
% %     theta_r_i =history_pre_ref.theta(i);
% %     A_cell{i} = [1, 0, -dt*V_target*sin(theta_r_i);
% %                  0, 1,  dt*V_target*cos(theta_r_i);
% %                  0, 0, 1];
% %     B_cell{i} = [dt*cos(theta_r_i), 0;
% %                  dt*sin(theta_r_i), 0;
% %                  0,               dt];
% % end
% 
%构造时变 PHI (Nx*Np x (Nx+Nu))
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
%%加入 d_k 的自由误差预测
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

    % 由于目前圆形：
    % v_ref = 常数
    % omega_ref =常数
    % 所以参考输入没有变化：
    % delta_u_ref = 0
    % 增广扰动就是：
    g_i = [ ...
        d_ref(:,i);
        zeros(Nu,1)];

    %% 在假设未来 ΔU = 0 的情况下传播误差
    Xi_free = ...
        A_tilde_i*Xi_free + g_i;

    %% 取出位置、航向三个误差
    error_free((i-1)*Nx+1:i*Nx) = ...
        C_tilde*Xi_free;

end
%% 目标函数
H_cell=cell(2,2);
H_cell{1,1}=THETA'*Q*THETA+R;
H_cell{1,2}=zeros(Nu*Nc,1);
H_cell{2,1}=zeros(1,Nu*Nc);
H_cell{2,2}=Row;
H=cell2mat(H_cell);
H=(H+H')/2;

% error=PHI*kesi;%???
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
umin = [-1.2; -4.0];       % 放宽角速度下限，缩窄速度下限
umax = [1.2;  4.0];        % 放宽两个速度上限，缩窄速度上限
% umin = [-0.8; -10.0];       % 放宽角速度下限，缩窄速度下限
% umax = [0.8;  10.0];        % 放宽两个速度上限，缩窄速度上限
delta_umin = [-0.2; -0.7];  % 放宽增量下限
delta_umax = [0.2;  0.7];   % 放宽增量上限

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
history_ref.x(num) = x_ref;
history_ref.y(num) = y_ref;
history.x(num) = x;
history.y(num) = y;
history.theta(num) = theta;
history.v(num)=v;
history.omega(num)=omega;
% time(num) = t;
history.t(num) = t;
num = num+1;
end
%% ================== 圆轨迹参数 ==================
R_circle = 3;
circle_center = [0, R_circle];

% +1：逆时针
% -1：顺时针
turn_dir = 1;

% 理论参考角速度
omega_ref = turn_dir * V_target / R_circle;


%% ================== 矢量合成 PID 参数 ==================

% ----- 矢量纠偏权重 PID -----
% 输入：横向误差绝对值
% 输出：切向量与纠偏向量的混合比例 beta
params_vec.Kp = 3.0;
params_vec.Ki = 0.0;

% 原来0.8稍微有些大，因为这里D项控制的是“权重”
% 可以先降低一些，避免beta变化过于剧烈
params_vec.Kd = 0.2;

% 最大纠偏权重
% 不建议达到1，否则可能完全丢掉前进方向
params_vec.MaxOutput = 0.8;


% ----- 航向 PID -----
params_head.Kp = 4.0;
params_head.Ki = 0.0;
params_head.Kd = 0.5;

% 航向反馈角速度最大值
params_head.MaxOutput = 4.0;

% 总角速度限制
omega_max = 4.0;


%% ================== PID状态初始化 ==================

state_vec.integral_sum = 0;
state_vec.error_prev = 0;

state_head_A.integral_sum = 0;
state_head_A.error_prev = 0;


%% ================== 车辆状态 ==================

x = x0;
y = y0;
theta = theta0;


%% ================== 数据记录 ==================

history_vec.x = [];
history_vec.y = [];
history_vec.error = [];
history_vec.pro = [];

history_pid.v = [];
history_pid.omega = [];
history_pid.omega_ff = [];
history_pid.omega_fb = [];
history_pid.t = [];


%% ================== PID仿真 ==================

for t = 0:dt:T_max

    %% 1. 当前车辆相对于圆心的极角

    phi = atan2( ...
        y - circle_center(2), ...
        x - circle_center(1));


    %% 2. 当前车辆在参考圆上的投影点

    x_ref = circle_center(1) ...
          + R_circle*cos(phi);

    y_ref = circle_center(2) ...
          + R_circle*sin(phi);


    %% 3. 圆的切线方向

    % 逆时针：
    % path_angle = phi + pi/2
    %
    % 顺时针：
    % path_angle = phi - pi/2

    path_angle = ...
        phi + turn_dir*pi/2;

    path_angle = wrapToPi(path_angle);

    T_vec = [ ...
        cos(path_angle), ...
        sin(path_angle)];


    %% 4. 横向误差

    distance_to_center = sqrt( ...
        (x-circle_center(1))^2 ...
      + (y-circle_center(2))^2);

    % e_lat > 0 ：车辆在圆外
    % e_lat < 0 ：车辆在圆内
    e_lat = ...
        distance_to_center - R_circle;


    %% 5. 计算径向单位向量

    radial_vec = [ ...
        x-circle_center(1), ...
        y-circle_center(2)];

    radial_norm = norm(radial_vec);

    if radial_norm > 1e-8
        radial_vec = ...
            radial_vec/radial_norm;
    else
        radial_vec = [1,0];
    end


    %% 6. 计算纠偏方向 F_vec

    %
    % radial_vec：从圆心指向车辆，向外
    %
    % 圆外：
    % e_lat > 0
    % F_vec = -radial_vec
    % → 向圆心
    %
    % 圆内：
    % e_lat < 0
    % F_vec = +radial_vec
    % → 向圆外
    %

    if abs(e_lat) > 1e-6

        F_vec = ...
            -sign(e_lat)*radial_vec;

    else

        % 已经在圆上，无需横向修正
        F_vec = [0,0];

    end


    %% 7. PID计算纠偏权重 beta

    %
    % 注意：
    %
    % PID这里只控制“纠偏程度”
    % 因此输入的是 abs(e_lat)
    %
    % 方向已经由 F_vec 决定
    %

    [beta_raw, state_vec] = ...
        PID_Controller( ...
            abs(e_lat), ...
            0, ...
            dt, ...
            params_vec, ...
            state_vec);

    % 限制为非负权重
    beta = max(beta_raw,0);

    % 再进行最大值限制
    beta = min(beta, ...
               params_vec.MaxOutput);

    % 已经非常接近参考轨迹时直接取消纠偏
    if abs(e_lat) < 1e-4
        beta = 0;
    end

    history_vec.pro(end+1) = beta;


    %% 8. 切向运动 + 横向纠偏矢量合成

    %
    % beta = 0：
    % 完全沿切线运动
    %
    % beta增大：
    % 越来越倾向横向纠偏
    %

    E_vec_raw = ...
        (1-beta)*T_vec ...
        + beta*F_vec;


    %% 9. 单位化期望运动方向

    if norm(E_vec_raw) > 1e-8

        E_vec = ...
            E_vec_raw/norm(E_vec_raw);

    else

        E_vec = T_vec;

    end


    %% 10. 得到期望速度方向

    V_cmd_world = ...
        E_vec*V_target;

    theta_cmd = atan2( ...
        V_cmd_world(2), ...
        V_cmd_world(1));


    %% 11. 航向误差

    angle_diff = ...
        wrapToPi(theta_cmd-theta);


    %% 12. 根据航向误差调节速度

    %
    % 航向误差较大时降低前进速度
    %

    v = ...
        max(V_target*cos(angle_diff),0);


    %% 13. 航向反馈PID

    [omega_fb, state_head_A] = ...
        PID_Controller( ...
            angle_diff, ...
            0, ...
            dt, ...
            params_head, ...
            state_head_A);


    %% 14. 圆轨迹曲率前馈

    %
    % 对圆：
    %
    % kappa = 1/R
    %
    % omega_ff = v*kappa
    %

    omega_ff = ...
        turn_dir*v/R_circle;


    %% 15. 总角速度

    omega = ...
        omega_ff + omega_fb;

    omega = max( ...
        min(omega,omega_max), ...
        -omega_max);


    %% 16. 左右轮速度

    vr = v + omega*L/2;
    vl = v - omega*L/2;


    %% 17. 实际车辆速度

    v_act = (vr+vl)/2;

    w_act = (vr-vl)/L;


    %% 18. 更新车辆运动状态

    x = ...
        x + v_act*cos(theta)*dt;

    y = ...
        y + v_act*sin(theta)*dt;

    theta = ...
        wrapToPi(theta+w_act*dt);


    %% 19. 保存数据

    history_vec.x(end+1) = x;
    history_vec.y(end+1) = y;

    history_vec.error(end+1) = ...
        e_lat;

    history_pid.v(end+1) = v;

    history_pid.omega(end+1) = ...
        omega;

    history_pid.omega_ff(end+1) = ...
        omega_ff;

    history_pid.omega_fb(end+1) = ...
        omega_fb;

    history_pid.t(end+1) = t;

end

% %% ================== 圆轨迹参数 ==================
% R_circle = 3;
% circle_center = [0, R_circle];
% 
% % +1：逆时针
% % -1：顺时针
% turn_dir = 1;
% 
% % 理论参考角速度
% omega_ref = turn_dir * V_target / R_circle;
% 
% 
% %% ================== 矢量合成 PID 参数 ==================
% 
% % ----- 矢量纠偏权重 PID -----
% % 输入：横向误差绝对值
% % 输出：切向量与纠偏向量的混合比例 beta
% params_vec.Kp = 3.0;
% params_vec.Ki = 0.0;
% 
% % 原来0.8稍微有些大，因为这里D项控制的是“权重”
% % 可以先降低一些，避免beta变化过于剧烈
% params_vec.Kd = 0.2;
% 
% % 最大纠偏权重
% % 不建议达到1，否则可能完全丢掉前进方向
% params_vec.MaxOutput = 0.8;
% 
% 
% % ----- 航向 PID -----
% params_head.Kp = 4.0;
% params_head.Ki = 0.0;
% params_head.Kd = 0.5;
% 
% % 航向反馈角速度最大值
% params_head.MaxOutput = 4.0;
% 
% % 总角速度限制
% omega_max = 4.0;
% 
% 
% %% ================== PID状态初始化 ==================
% 
% state_vec.integral_sum = 0;
% state_vec.error_prev = 0;
% 
% state_head_A.integral_sum = 0;
% state_head_A.error_prev = 0;
% 
% 
% %% ================== 车辆状态 ==================
% 
% x = x0;
% y = y0;
% theta = theta0;
% 
% 
% %% ================== 数据记录 ==================
% 
% history_vec.x = [];
% history_vec.y = [];
% history_vec.error = [];
% history_vec.pro = [];
% 
% history_pid.v = [];
% history_pid.omega = [];
% history_pid.omega_ff = [];
% history_pid.omega_fb = [];
% history_pid.t = [];
% 
% 
% %% ================== PID仿真 ==================
% 
% for t = 0:dt:T_max
% 
%     %% 1. 当前车辆相对于圆心的极角
% 
%     phi = atan2( ...
%         y - circle_center(2), ...
%         x - circle_center(1));
% 
% 
%     %% 2. 当前车辆在参考圆上的投影点
% 
%     x_ref = circle_center(1) ...
%           + R_circle*cos(phi);
% 
%     y_ref = circle_center(2) ...
%           + R_circle*sin(phi);
% 
% 
%     %% 3. 圆的切线方向
% 
%     % 逆时针：
%     % path_angle = phi + pi/2
%     %
%     % 顺时针：
%     % path_angle = phi - pi/2
% 
%     path_angle = ...
%         phi + turn_dir*pi/2;
% 
%     path_angle = wrapToPi(path_angle);
% 
%     T_vec = [ ...
%         cos(path_angle), ...
%         sin(path_angle)];
% 
% 
%     %% 4. 横向误差
% 
%     distance_to_center = sqrt( ...
%         (x-circle_center(1))^2 ...
%       + (y-circle_center(2))^2);
% 
%     % e_lat > 0 ：车辆在圆外
%     % e_lat < 0 ：车辆在圆内
%     e_lat = ...
%         distance_to_center - R_circle;
% 
% 
%     %% 5. 计算径向单位向量
% 
%     radial_vec = [ ...
%         x-circle_center(1), ...
%         y-circle_center(2)];
% 
%     radial_norm = norm(radial_vec);
% 
%     if radial_norm > 1e-8
%         radial_vec = ...
%             radial_vec/radial_norm;
%     else
%         radial_vec = [1,0];
%     end
% 
% 
%     %% 6. 计算纠偏方向 F_vec
% 
%     %
%     % radial_vec：从圆心指向车辆，向外
%     %
%     % 圆外：
%     % e_lat > 0
%     % F_vec = -radial_vec
%     % → 向圆心
%     %
%     % 圆内：
%     % e_lat < 0
%     % F_vec = +radial_vec
%     % → 向圆外
%     %
% 
%     if abs(e_lat) > 1e-6
% 
%         F_vec = ...
%             -sign(e_lat)*radial_vec;
% 
%     else
% 
%         % 已经在圆上，无需横向修正
%         F_vec = [0,0];
% 
%     end
% 
% 
%     %% 7. PID计算纠偏权重 beta
% 
%     %
%     % 注意：
%     %
%     % PID这里只控制“纠偏程度”
%     % 因此输入的是 abs(e_lat)
%     %
%     % 方向已经由 F_vec 决定
%     %
% 
%     [beta_raw, state_vec] = ...
%         PID_Controller( ...
%             abs(e_lat), ...
%             0, ...
%             dt, ...
%             params_vec, ...
%             state_vec);
% 
%     % 限制为非负权重
%     beta = max(beta_raw,0);
% 
%     % 再进行最大值限制
%     beta = min(beta, ...
%                params_vec.MaxOutput);
% 
%     % 已经非常接近参考轨迹时直接取消纠偏
%     if abs(e_lat) < 1e-4
%         beta = 0;
%     end
% 
%     history_vec.pro(end+1) = beta;
% 
% 
%     %% 8. 切向运动 + 横向纠偏矢量合成
% 
%     %
%     % beta = 0：
%     % 完全沿切线运动
%     %
%     % beta增大：
%     % 越来越倾向横向纠偏
%     %
% 
%     E_vec_raw = ...
%         (1-beta)*T_vec ...
%         + beta*F_vec;
% 
% 
%     %% 9. 单位化期望运动方向
% 
%     if norm(E_vec_raw) > 1e-8
% 
%         E_vec = ...
%             E_vec_raw/norm(E_vec_raw);
% 
%     else
% 
%         E_vec = T_vec;
% 
%     end
% 
% 
%     %% 10. 得到期望速度方向
% 
%     V_cmd_world = ...
%         E_vec*V_target;
% 
%     theta_cmd = atan2( ...
%         V_cmd_world(2), ...
%         V_cmd_world(1));
% 
% 
%     %% 11. 航向误差
% 
%     angle_diff = ...
%         wrapToPi(theta_cmd-theta);
% 
% 
%     %% 12. 根据航向误差调节速度
% 
%     %
%     % 航向误差较大时降低前进速度
%     %
% 
%     v = ...
%         max(V_target*cos(angle_diff),0);
% 
% 
%     %% 13. 航向反馈PID
% 
%     [omega_fb, state_head_A] = ...
%         PID_Controller( ...
%             angle_diff, ...
%             0, ...
%             dt, ...
%             params_head, ...
%             state_head_A);
% 
% 
%     %% 14. 圆轨迹曲率前馈
% 
%     %
%     % 对圆：
%     %
%     % kappa = 1/R
%     %
%     % omega_ff = v*kappa
%     %
% 
%     omega_ff = ...
%         turn_dir*v/R_circle;
% 
% 
%     %% 15. 总角速度
% 
%     omega = ...
%         omega_ff + omega_fb;
% 
%     omega = max( ...
%         min(omega,omega_max), ...
%         -omega_max);
% 
% 
%     %% 16. 左右轮速度
% 
%     vr = v + omega*L/2;
%     vl = v - omega*L/2;
% 
% 
%     %% 17. 实际车辆速度
% 
%     v_act = (vr+vl)/2;
% 
%     w_act = (vr-vl)/L;
% 
% 
%     %% 18. 更新车辆运动状态
% 
%     x = ...
%         x + v_act*cos(theta)*dt;
% 
%     y = ...
%         y + v_act*sin(theta)*dt;
% 
%     theta = ...
%         wrapToPi(theta+w_act*dt);
% 
% 
%     %% 19. 保存数据
% 
%     history_vec.x(end+1) = x;
%     history_vec.y(end+1) = y;
% 
%     history_vec.error(end+1) = ...
%         e_lat;
% 
%     history_pid.v(end+1) = v;
% 
%     history_pid.omega(end+1) = ...
%         omega;
% 
%     history_pid.omega_ff(end+1) = ...
%         omega_ff;
% 
%     history_pid.omega_fb(end+1) = ...
%         omega_fb;
% 
%     history_pid.t(end+1) = t;
% 
% end
% %% ================== 圆轨迹参数 ==================
% R_circle = 3;
% circle_center = [0, R_circle];
% 
% % 理论圆周角速度
% omega_ref = V_target / R_circle;
% 
% %% ================== PID 参数 ==================
% 
% % ---------- 横向误差 PID ----------
% % 输出：航向修正角 angle_offset
% params_pid_lat.Kp = 3.0;
% params_pid_lat.Ki = 0.0;
% params_pid_lat.Kd = 0.8;
% params_pid_lat.MaxOutput = pi/3;   % 最大修正航向 ±60°
% 
% % ---------- 航向 PID ----------
% % 输出：反馈角速度 omega_fb
% params_head.Kp = 4.0;
% params_head.Ki = 0.0;
% params_head.Kd = 0.5;
% params_head.MaxOutput = 4.0;
% 
% %% ================== PID 状态初始化 ==================
% state_lat.integral_sum = 0;
% state_lat.error_prev = 0;
% 
% state_head.integral_sum = 0;
% state_head.error_prev = 0;
% 
% x = x0;
% y = y0;
% theta = theta0;
% 
% %% ================== 数据记录 ==================
% history_pid.x = [];
% history_pid.y = [];
% history_pid.error = [];
% history_pid.v = [];
% history_pid.omega = [];
% history_pid.omega_ff = [];
% history_pid.omega_fb = [];
% history_pid.t = [];
% 
% %% ================== PID 仿真 ==================
% for t = 0:dt:T_max
% 
%     %% 1. 当前车辆相对于圆心的极角
%     phi = atan2( ...
%         y - circle_center(2), ...
%         x - circle_center(1));
% 
%     %% 2. 当前车辆在圆上的几何投影点
%     x_ref = circle_center(1) + R_circle*cos(phi);
%     y_ref = circle_center(2) + R_circle*sin(phi);
% 
%     %% 3. 圆轨迹切线方向
%     % 当前圆为逆时针运动
%     path_angle = wrapToPi(phi + pi/2);
% 
%     %% 4. 横向误差
%     distance_to_center = sqrt( ...
%         (x - circle_center(1))^2 + ...
%         (y - circle_center(2))^2);
% 
%     % > 0：车辆在圆外
%     % < 0：车辆在圆内
%     e_lat = distance_to_center - R_circle;
% 
%     %% 5. 横向 PID -> 航向修正角
%     %
%     % e_lat > 0：
%     % 车辆在圆外，需要向圆心方向偏转
%     %
%     % e_lat < 0：
%     % 车辆在圆内，需要向圆外方向偏转
%     %
%     [angle_offset, state_lat] = ...
%         PID_Controller( ...
%             e_lat, ...
%             0, ...
%             dt, ...
%             params_pid_lat, ...
%             state_lat);
% 
%     %% 6. 最终期望航向
%     theta_cmd = wrapToPi(path_angle + angle_offset);
% 
%     %% 7. 航向误差
%     angle_diff = wrapToPi(theta_cmd - theta);
% 
%     %% 8. 根据航向误差调节前进速度
%     %
%     % 航向偏差越大，速度越低
%     %
%     v = max(V_target*cos(angle_diff), 0);
% 
%     %% 9. 航向反馈 PID
%     [omega_fb, state_head] = ...
%         PID_Controller( ...
%             angle_diff, ...
%             0, ...
%             dt, ...
%             params_head, ...
%             state_head);
% 
%     %% 10. 圆轨迹曲率前馈
%     %
%     % 圆周运动：
%     %
%     %       omega = v / R
%     %
%     % 即使横向误差和航向误差都为0，
%     % 车辆仍然需要保持非零角速度才能沿圆运动
%     %
%     omega_ff = v / R_circle;
% 
%     %% 11. 总角速度
%     omega = omega_ff + omega_fb;
% 
%     % 角速度限幅
%     omega = max(min(omega, ...
%         params_head.MaxOutput), ...
%         -params_head.MaxOutput);
% 
%     %% 12. 差速轮速度
%     vr = v + omega*L/2;
%     vl = v - omega*L/2;
% 
%     %% 13. 实际车体速度
%     v_act = (vr + vl)/2;
%     w_act = (vr - vl)/L;
% 
%     %% 14. 更新车辆状态
%     x = x + v_act*cos(theta)*dt;
%     y = y + v_act*sin(theta)*dt;
% 
%     theta = wrapToPi(theta + w_act*dt);
% 
%     %% 15. 保存数据
%     history_pid.x(end+1) = x;
%     history_pid.y(end+1) = y;
% 
%     history_pid.error(end+1) = e_lat;
% 
%     history_pid.v(end+1) = v;
%     history_pid.omega(end+1) = omega;
% 
%     history_pid.omega_ff(end+1) = omega_ff;
%     history_pid.omega_fb(end+1) = omega_fb;
% 
%     history_pid.t(end+1) = t;
% 
% end
% % if draw == true
% % % 绘图：'-o' 表示实线+圆圈，'linewidth' 加粗线条，'markersize' 加大点
% % plot(history_ref.x, history_ref.y, '-o', 'LineWidth', 1.5, 'MarkerSize', 2, 'MarkerFaceColor', 'b');
% % 
% % 添加网格和标签
% grid on;
% xlabel('X 轴');
% ylabel('Y 轴');
% title('坐标点可视化');
% end
% 
if draw == true
     % 绘制第一组曲线（参考轨迹），线宽调细
    plot(history_ref.x, history_ref.y, '-o', 'LineWidth', 0.8, ...
         'MarkerSize', 2, 'Color', [0.00, 0.45, 0.74], 'DisplayName', '参考轨迹');
    hold on;  
    % 绘制第三组曲线（新轨迹），线宽同样调细
    plot(history.x, history.y, '-s', 'LineWidth', 0.8, ...
         'MarkerSize', 2, 'Color', [0.85, 0.33, 0.10], 'DisplayName', 'mpc轨迹1');
      % 绘制第二组曲线（新轨迹），线宽同样调细
    plot(history_vec.x, history_vec.y, '-s', 'LineWidth', 0.8, ...
         'MarkerSize', 2, 'Color', [0.47, 0.67, 0.19], 'DisplayName', 'pid轨迹');


    grid on;
    xlabel('X 轴');
    ylabel('Y 轴');
    title('坐标点可视化');
    legend('Location', 'best');
    hold off;
end
% pid和mpc的v，omega的对比
figure;
plot(history.t, history.omega, 'b','DisplayName', 'mpc omega'); hold on;
plot(history_pid.t, history_pid.omega, 'r','DisplayName', 'pid omega');
grid on;
xlabel('Time (s)');
ylabel('\omega (rad/s)');
title('\omega');
legend('MPC', 'PID');

figure;
plot(history.t, history.v, 'b','DisplayName', 'mpc omega'); hold on;
plot(history_pid.t, history_pid.v, 'r','DisplayName', 'pid omega');
grid on;
xlabel('Time (s)');
ylabel('v (m/s)');
title('v');
legend('MPC', 'PID');
% 绘制参考位置x，y的变化曲线。
figure;
plot(history.t, history_ref.x, 'b'); hold on;
plot(history.t, history.theta, 'g'); hold on;
grid on;
xlabel('Time (s)');
ylabel('m');
title('参考位置x');
plot(history.t, history_ref.y, 'r'); hold on;
grid on;
xlabel('Time (s)');
ylabel('m');
title('参考位置y');

% 
%% 6. 共用的 PID 控制器函数定义
function [u, pid_state] = PID_Controller(target, current, dt, params, pid_state)
    error = target - current;
    P = params.Kp * error;
    pid_state.integral_sum = pid_state.integral_sum + error * dt;
    limit_I = params.MaxOutput; 
    pid_state.integral_sum = max(min(pid_state.integral_sum, limit_I), -limit_I);
    I = params.Ki * pid_state.integral_sum;
    D = params.Kd * (error - pid_state.error_prev) / dt;
    u = P + I + D;
    u = max(min(u, params.MaxOutput), -params.MaxOutput); 
    pid_state.error_prev = error;
end

% 角度规范化辅助函数 (-pi 到 pi)
function angle = wrapToPi(angle)
    while angle > pi, angle = angle - 2*pi; end
    while angle < -pi, angle = angle + 2*pi; end
end

