%这是运用了时变预测矩阵的误差状态控制的mpc控制与动态避障仿真，涉及避障代价算法的改进,和状态机的改进,重写一个轨迹回归的状态机
%这个代码将会优化掉一些无用的错误思路，加强可读性
%参考轨迹为正方形轨迹
clear; close all; clc;

%% 1. 共有参数设置
L = 0.5;              % 轮距 (m)
dt = 0.01;             % 采样时间 (s)
T = 0.05;
T_max = 50;           % 最大模拟时长 (s)
V_target = 0.8;       % 期望巡航速度 (m/s)
% LookAhead_Dist = 0.2; % 路径点切换阈值 (m) (判断是否到达航点的距离)(mpc不需要）
%在正方形的轨迹里面还是需要的，但是新的mpc代码没有用到它，作者（就是我）用了0.02提升了采样精度
Np = 60;              %预测步长
Nc = 30;              %控制步长

plan_Np =20;          %规划器的预测步长
plan_Nc =7;           %规划器的控制步长，貌似是对于抑制规划轨迹弯道过大的效果比较好地参数了
%增广矩阵参数
Nx = 3;               %输出的增广维度参数
Nu = 2;               %输入的增广维度列数
% 初始位姿 (与前面算法统一)
x0 = 0; y0 = 0; theta0 = 0;%改动了初始位置，因为在没进入参考轨迹的情况下mpc控制发挥不好
squa_len = 5;
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
%障碍物设计
%圆形障碍物
% obs_pos = [5,2.5];
obs_pos = [0,2.5];
R_obs = 0.15;
% 双阈值的模式切换
d_trigger = 1.1;  % 进入避障
d_release = 1.1;   % 退出避障

draw = true;%绘图标志符，决定代码是否绘图
arri = true;%判断AGV是否到达轨迹起点，可以忽略，当起点处没有参考轨迹时mpc控制的效果是不好的，直接让mpc从车辆起点开始
t_arri = 0;%AGV到达轨迹起点的时间点
the_end = false;%参考的理想小车是否完成运动的标志
current_idx = 2;%正方形轨迹标签
Row=10;%松弛因子（可用可不用）
plan_draw =true;

v = 0;
omega = 0;

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
Q_plan = diag(linspace(100,300,plan_Np));
R_plan = 10*eye(plan_Nc,plan_Nc);
Q_terminal_avoid = 500;%终端代价参数

Q_lat = 500 * eye(plan_Np);
Q_theta = 0 * eye(plan_Np);

Q_terminal_recover = 800;
Q_terminal_theta = 0;
%S_avoid，避障函数的权重。5-50效果都还不错
%S_recover，轨迹回归函数的权重
S_avoid = 5;
S_recover = 0;

%状态机参数
nornal = 1;
avoid = 2;
recover = 3;

mode = nornal;
%% 参考轨迹生成+求取状态误差量
for t = 0:dt:T_max
    
    if mode == nornal
        pre_t = Np*dt;%控制的预测时长
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
if  mode == nornal
    x_ref = history_pre_ref.x(1);
    y_ref = history_pre_ref.y(1);
    theta_ref = history_pre_ref.theta(1);
    
elseif mode == avoid
    %避障模式，约束有改进空间
    omega_max = 6;
    omega_min = -6; 
    lb = omega_min * ones(plan_Nc,1);
    ub = omega_max * ones(plan_Nc,1);
    A=[];
    b=[];
    Aeq=[];
    beq=[];
    options = optimset('Algorithm','active-set');
    % Yref = history_pre_ref.y的插值
     % 原参考轨迹：dt = 0.01 s
    t_old = (1:length(history_pre_ref.y)) * dt;
    % 规划器：T = 0.05 s，共 plan_Np 个点
    t_new = (1:plan_Np) * T;
    % 将原参考轨迹重采样到规划器预测时域
    Xref = interp1(t_old, history_pre_ref.x, t_new, 'linear')';
    Yref = interp1(t_old, history_pre_ref.y, t_new, 'linear')';
    fprintf("length(Yref) = %d\n",length(Yref));
    %规划层的约束？现在还没实验完全
    [omega_opt,fval,exitflag]=fmincon(@(target)function_avoid(target,x,y,theta,omega,v,plan_Np, ...
        plan_Nc,T,Yref,Xref,Q_plan,R_plan,S_avoid,Q_terminal_avoid,obs_pos,R_obs),zeros(plan_Nc,1),A,b,Aeq,beq,lb,ub,[],options);
    fprintf('exitflag=%d\n',exitflag);
    %生成重规划轨迹
    % 以下根据计算出的控制量推导所有的状态量
    y_pre=zeros(plan_Np,1);
    x_pre=zeros(plan_Np,1);
    theta_pre=zeros(plan_Np,1);
    for i=1:1:plan_Np
            if i==1 
                omega_pre(i,1)=omega_opt(1);
                 % 以下完成状态量更新，参考原论文代码
                % x = x + v_act * cos(theta) * dt; y = y + v_act * sin(theta) * dt; 
                % theta = wrapToPi(theta + w_act * dt);
                theta_pre(i,1)=wrapToPi(theta+omega_pre(i)*T);
                x_pre(i,1)=x+v*cos(theta_pre(i,1))*T;
                y_pre(i,1)=y+v*sin(theta_pre(i,1))*T;
            elseif i<=plan_Nc 
                omega_pre(i,1)=omega_opt(i);
                theta_pre(i,1)=wrapToPi(theta_pre(i-1,1)+omega_pre(i)*T);
                x_pre(i,1)=x_pre(i-1,1)+v*cos(theta_pre(i,1))*T;
                y_pre(i,1)=y_pre(i-1,1)+v*sin(theta_pre(i,1))*T;
            else
                omega_pre(i,1)=omega_opt(plan_Nc);
                theta_pre(i,1)=wrapToPi(theta_pre(i-1,1)+omega_pre(i)*T);
                x_pre(i,1)=x_pre(i-1,1)+v*cos(theta_pre(i,1))*T;
                y_pre(i,1)=y_pre(i-1,1)+v*sin(theta_pre(i,1))*T;
            end
    end
    %这个元胞用来存储规划轨迹的数据，用于绘制规划路径的
    history.plan_x_test{plan_num,1} = x_pre;
    history.plan_y_test{plan_num,1} = y_pre;
    plan_num = plan_num+1;
   
     % 全插值绘制，因为规划序列能看到1秒的情况，而预测序列只能看到0.6秒，要用插值法把规划序列变得与预测序列长度一致
     %在设计上规划轨迹的mpc应该在时间上“看”得比跟踪mpc一样甚至更远
    t_plan = [0; (1:plan_Np)' * T];
    x_plan = [x; x_pre];
    y_plan = [y; y_pre];
    theta_plan = [theta; theta_pre];
    t_ctrl = (1:Np)' * dt;
    %测试插值的代码
    % size(t_plan)
    % size(x_plan)
    % size(y_plan)
    % size(theta_plan)
    % size(t_ctrl)
    %插值后得到的新的规划轨迹预测序列，用来给跟踪mpc跟踪
    x_pre_ref = interp1(t_plan, x_plan, t_ctrl);
    y_pre_ref = interp1(t_plan, y_plan, t_ctrl);
    theta_pre_ref = interp1(t_plan, theta_plan, t_ctrl);
    
    x_ref = history_pre_ref.x(1);
    y_ref = history_pre_ref.y(1);
    theta_ref = history_pre_ref.theta(1);
    x_ref_obs = x_pre_ref(1);
    y_ref_obs = y_pre_ref(1);
    theta_ref_obs = theta_pre_ref(1);
    omega_ref_obs = omega_pre(1);
elseif mode ==recover
    %回归模式，约束有改进空间
    omega_max = 8;
    omega_min = -8; 
    lb = omega_min * ones(plan_Nc,1);
    ub = omega_max * ones(plan_Nc,1);
    A=[];
    b=[];
    Aeq=[];
    beq=[];
    options = optimset('Algorithm','active-set');
     % 原参考轨迹：dt = 0.01 s
    t_old = (1:length(history_pre_ref.y)) * dt;
    % 规划器：T = 0.05 s，共 plan_Np 个点
    t_new = (1:plan_Np) * T;
    % 将原参考轨迹重采样到规划器预测时域
    Xref = interp1(t_old, history_pre_ref.x, t_new, 'linear')';
    Yref = interp1(t_old, history_pre_ref.y, t_new, 'linear')';
    THETAref = interp1(t_old, history_pre_ref.theta, t_new, 'linear')';
    fprintf("length(Yref) = %d\n",length(Yref));
    %规划层的约束？
    % Q_lat = 300 * eye(plan_Np);
    % Q_theta = 100 * eye(plan_Np);
    % 
    % Q_terminal_recover = 800;
    % Q_terminal_theta = 300;
    [omega_opt,fval,exitflag]=fmincon(@(target)function_recover(target,x,y,theta,omega,v,plan_Np, ...
        plan_Nc,T,Yref,Xref,THETAref,Q_lat,Q_theta,R_plan,S_recover,Q_terminal_recover,Q_terminal_theta,obs_pos,R_obs),zeros(plan_Nc,1),A,b,Aeq,beq,lb,ub,[],options);
    fprintf('exitflag=%d\n',exitflag);
    %生成重规划轨迹
    % 以下根据计算出的控制量推导所有的状态量
    y_pre=zeros(plan_Np,1);
    x_pre=zeros(plan_Np,1);
    theta_pre=zeros(plan_Np,1);
    for i=1:1:plan_Np
            if i==1 
                omega_pre(i,1)=omega_opt(1);
                % theta_pre(i,1)=wrapToPi(theta+omega_pre(i)*T);
                theta_pre(i,1)=(theta+omega_pre(i)*T);
                x_pre(i,1)=x+v*cos(theta_pre(i,1))*T;
                y_pre(i,1)=y+v*sin(theta_pre(i,1))*T;
            elseif i<=plan_Nc 
                omega_pre(i,1)=omega_opt(i);
                % theta_pre(i,1)=wrapToPi(theta_pre(i-1,1)+omega_pre(i)*T);
                theta_pre(i,1)=(theta_pre(i-1,1)+omega_pre(i)*T);
                x_pre(i,1)=x_pre(i-1,1)+v*cos(theta_pre(i,1))*T;
                y_pre(i,1)=y_pre(i-1,1)+v*sin(theta_pre(i,1))*T;
                
            else
                omega_pre(i,1)=omega_opt(plan_Nc);
                % theta_pre(i,1)=wrapToPi(theta_pre(i-1,1)+omega_pre(i)*T);
                theta_pre(i,1)=(theta_pre(i-1,1)+omega_pre(i)*T);
                x_pre(i,1)=x_pre(i-1,1)+v*cos(theta_pre(i,1))*T;
                y_pre(i,1)=y_pre(i-1,1)+v*sin(theta_pre(i,1))*T;
            end
    end
    history.plan_x_test{plan_num,1} = x_pre;
    history.plan_y_test{plan_num,1} = y_pre;
    plan_num = plan_num+1;
   
    t_plan = [0; (1:plan_Np)' * T];
    x_plan = [x; x_pre];
    y_plan = [y; y_pre];
    theta_plan = [theta; theta_pre];
    t_ctrl = (1:Np)' * dt;
    %插值后得到的新的规划轨迹预测序列，用来给跟踪mpc跟踪
    x_pre_ref = interp1(t_plan, x_plan, t_ctrl);
    y_pre_ref = interp1(t_plan, y_plan, t_ctrl);
    theta_pre_ref = interp1(t_plan, theta_plan, t_ctrl);
    
    x_ref = history_pre_ref.x(1);
    y_ref = history_pre_ref.y(1);
    theta_ref = history_pre_ref.theta(1);
    x_ref_obs = x_pre_ref(1);
    y_ref_obs = y_pre_ref(1);
    theta_ref_obs = theta_pre_ref(1);
    omega_ref_obs = omega_pre(1);
end

%% 误差向量的求取
%误差向量求取
if mode ==nornal
kesi(1)= x-x_ref;
kesi(2)= y-y_ref;
kesi(3)= wrapToPi(theta-theta_ref);
% kesi(3)= (theta-theta_ref);
kesi(4)= v-V_target;
kesi(5)= omega-omega_ref;
else
kesi(1)= x-x_ref_obs;
kesi(2)= y-y_ref_obs;
kesi(3)= wrapToPi(theta-theta_ref_obs);
% kesi(3)= (theta-theta_ref);
kesi(4)= v-V_target;
kesi(5)= omega-omega_ref_obs; 
% kesi(5)= omega-omega_ref; 
end

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
    if mode == nornal
        theta_r_i =history_pre_ref.theta(i);
    else
        theta_r_i = theta_pre_ref(i);
    end
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
        zeros(Nu,1)];

    %% 在假设未来 ΔU = 0 的情况下传播误差
    Xi_free = ...
        A_tilde_i*Xi_free + g_i;

    %% 取出位置、航向三个误差
    error_free((i-1)*Nx+1:i*Nx) = ...
        C_tilde*Xi_free;
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
d_obs_dot = [x-obs_pos(1),y-obs_pos(2)];
history.d_obs_dot(num)=norm(d_obs_dot);
tangent = [cos(theta_ref), sin(theta_ref)];

normal  = [-sin(theta_ref), cos(theta_ref)];
e_vec = [x-x_ref, y-y_ref];
e_long = dot(e_vec, tangent);   % 纵向误差
e_lat  = dot(e_vec, normal);    % 横向误差

% %第二版状态机
% switch mode
%     case nornal
%         if norm(d_obs_dot) <= R_obs+d_trigger
%             mode = avoid;
%         end
%     case avoid
%         a =[cos(theta_ref),sin(theta_ref)];
%         b =[x-obs_pos(1),y-obs_pos(2)];
%         if a*b'>= 0.1
%             mode = recover;
%         end
%     case recover
%         if norm(e_lat)<=0.1
%             mode = nornal;
%         end
% end
%% 第三版状态机
switch mode

    case nornal

        tangent = [cos(theta_ref), sin(theta_ref)];

        obs_relative = [obs_pos(1)-x, ...
                        obs_pos(2)-y];

        s_obs = tangent * obs_relative';
        % 必须同时满足：
        % 1. 距离足够近
        % 2. 障碍物仍然在车辆前方
        if norm(d_obs_dot) <= R_obs+d_trigger && s_obs > 0
            mode = avoid;
        end


    case avoid

        tangent = [cos(theta_ref), sin(theta_ref)];

        vehicle_relative = [x-obs_pos(1), ...
                            y-obs_pos(2)];

        s_pass = tangent * vehicle_relative';

        if s_pass >= 0.1
            % break;
            mode = recover;
        end


    case recover

        if abs(e_lat) <= 0.05
            mode = nornal;
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
    plot(history_ref.x, history_ref.y, '-o', 'LineWidth', 0.8, ...
         'MarkerSize', 2, 'MarkerFaceColor', 'b', 'DisplayName', '参考轨迹');

    hold on;  
    % 绘制第二组曲线（新轨迹），线宽同样调细
    plot(history.x, history.y, '-s', 'LineWidth', 0.8, ...
         'MarkerSize', 2, 'MarkerFaceColor', 'r', 'DisplayName', 'mpc轨迹1');
if plan_draw ==true
    % 绘制第三组曲线（规划轨迹），线宽同样调细
     for i = 1:1:length(history.plan_x_test)
    plot(history.plan_x_test{i,1}, history.plan_y_test{i,1}, '-s', 'LineWidth', 0.8, ...
         'MarkerSize', 2, 'MarkerFaceColor', 'g');
     end
end
     
     %% 3. 绘制圆形障碍物
    theta_circle = linspace(0, 2*pi, 200);

    x_obs = obs_pos(1) + R_obs*cos(theta_circle);
    y_obs = obs_pos(2) + R_obs*sin(theta_circle);

    fill(x_obs, y_obs, 'k',...
        'DisplayName', '障碍物');

    grid on;
    xlabel('X 轴');
    ylabel('Y 轴');
    title('坐标点可视化');
    legend('Location', 'best');
    hold off;
end
if draw == true
    figure;
    % 绘制第一组曲线（参考轨迹），线宽调细
    plot(history_ref.x, history_ref.y, '-o', 'LineWidth', 0.8, ...
         'MarkerSize', 2, 'Color', [0.00, 0.45, 0.74], 'DisplayName', '参考轨迹');
    hold on;  
    % 绘制第三组曲线（新轨迹），线宽同样调细
    plot(history.x, history.y, '-s', 'LineWidth', 0.8, ...
         'MarkerSize', 2, 'Color', [0.85, 0.33, 0.10], 'DisplayName', 'mpc轨迹');

     %% 3. 绘制圆形障碍物
    theta_circle = linspace(0, 2*pi, 200);

    x_obs = obs_pos(1) + R_obs*cos(theta_circle);
    y_obs = obs_pos(2) + R_obs*sin(theta_circle);

    fill(x_obs, y_obs, 'k',...
        'DisplayName', '障碍物');

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


%% 用于避障轨迹规划的代价函数
function cost =function_avoid(target,x0,y0,theta0,omega0,v0,plan_Np,plan_Nc,dt,Yref,Xref,Q,R,S,Q_terminal,obs_pos,R_obs) 
%参量初始化
cost=0;
x = x0;
y = y0;
theta = theta0;
omega = omega0;
v = v0;
%预测序列定义
x_pre = zeros(plan_Np,1);
y_pre = zeros(plan_Np,1);
theta_pre = zeros(plan_Np,1);
J_obst = zeros(plan_Np,1);
Y_error = zeros(plan_Np,1);
X_error = zeros(plan_Np,1);
    for i=1:1:plan_Np
            if i==1 
                omega(i,1)=target(1);
                %以下完成状态量更新
                %模型参照：
                % x = x + v_act * cos(theta) * dt; y = y + v_act * sin(theta) * dt; 
                % theta = wrapToPi(theta + w_act * dt);
                theta_pre(i,1)=wrapToPi(theta+omega(i)*dt);
                x_pre(i,1)=x+v*cos(theta_pre(i,1))*dt;
                y_pre(i,1)=y+v*sin(theta_pre(i,1))*dt;
                % for j=1:1:Nobs
                %     J_obst(i,1)=J_obst(i,1)+1/(((X_predict(i,1))-X_obstacle(j,1))^2+(Y_predict(i,1)-Y_obstacle(j,1))^2+0.000001);
                % end
                %先用圆形障碍物，可以简化部分运算
                d = [x_pre(i,1)-obs_pos(1),y_pre(i,1)-obs_pos(2)];
                %计算距离车辆质心最近的障碍物点
                P = obs_pos + R_obs*d / norm(d);
                J_obst(i,1) = 1/(((x_pre(i,1)-P(1))^2+(y_pre(i,1)-P(2))^2)+0.000001);
            elseif i<=plan_Nc 
                omega(i,1)=target(i);
                theta_pre(i,1)=wrapToPi(theta_pre(i-1,1)+omega(i)*dt);
                x_pre(i,1)=x_pre(i-1,1)+v*cos(theta_pre(i,1))*dt;
                y_pre(i,1)=y_pre(i-1,1)+v*sin(theta_pre(i,1))*dt;
                d = [x_pre(i,1)-obs_pos(1),y_pre(i,1)-obs_pos(2)];
                P = obs_pos + R_obs*d / norm(d);
                J_obst(i,1) = 1/(((x_pre(i,1)-P(1))^2+(y_pre(i,1)-P(2))^2)+0.000001);
            else
                omega(i,1)=target(plan_Nc);
                theta_pre(i,1)=wrapToPi(theta_pre(i-1,1)+omega(i)*dt);
                x_pre(i,1)=x_pre(i-1,1)+v*cos(theta_pre(i,1))*dt;
                y_pre(i,1)=y_pre(i-1,1)+v*sin(theta_pre(i,1))*dt;

                d = [x_pre(i,1)-obs_pos(1),y_pre(i,1)-obs_pos(2)];
                P = obs_pos + R_obs*d / norm(d);
                J_obst(i,1) = 1/(((x_pre(i,1)-P(1))^2+(y_pre(i,1)-P(2))^2)+0.000001);
            end
        %J_obst=J_obst+1/(((X_predict(i,1))-X_obstacle(2,1))^2+(Y_predict(i,1)-Y_obstacle(2,1))^2+0.00001);
            Y_error(i,1)=y_pre(i,1)-Yref(i,1);%注意此处Yref与Y_ref要区分开来，Yref是局部期望路径,Y_ref为全局期望路径
            X_error(i,1)=x_pre(i,1)-Xref(i,1);
    end 
    terminal_error = ...
            (x_pre(plan_Nc)-Xref(plan_Nc))^2 + ...
            (y_pre(plan_Nc)-Yref(plan_Nc))^2;
        cost=cost+Y_error'*Q*Y_error+X_error'*Q*X_error+omega(1:plan_Nc)'*R*omega(1:plan_Nc)+...
        S*sum(J_obst(:))+Q_terminal*terminal_error;
end

%% 用于回归轨迹规划的代价函数
function cost =function_recover(target,x0,y0,theta0,omega0,v0,plan_Np,plan_Nc,dt,Yref,Xref,THETAref ...
    ,Q_lat,Q_theta,R,S,Q_terminal,Q_terminal_theta,obs_pos,R_obs) 
%参量初始化
cost=0;
x = x0;
y = y0;
theta = theta0;
omega = omega0;
v = v0;
%预测序列定义
THETA_error= zeros(plan_Np,1);
e_lat_plan = zeros(plan_Np,1);
x_pre = zeros(plan_Np,1);
y_pre = zeros(plan_Np,1);
theta_pre = zeros(plan_Np,1);
J_obst = zeros(plan_Np,1);
Y_error = zeros(plan_Np,1);
X_error = zeros(plan_Np,1);
    for i=1:1:plan_Np
            if i==1 
                omega(i,1)=target(1);
                % 以下完成状态量更新
                %模型参照：
                % x = x + v_act * cos(theta) * dt; y = y + v_act * sin(theta) * dt; 
                % theta = wrapToPi(theta + w_act * dt);
                theta_pre(i,1)=wrapToPi(theta+omega(i)*dt);
                x_pre(i,1)=x+v*cos(theta_pre(i,1))*dt;
                y_pre(i,1)=y+v*sin(theta_pre(i,1))*dt;
                % for j=1:1:Nobs
                %     J_obst(i,1)=J_obst(i,1)+1/(((X_predict(i,1))-X_obstacle(j,1))^2+(Y_predict(i,1)-Y_obstacle(j,1))^2+0.000001);
                % end
                %先用圆形障碍物，可以简化部分运算
                d = [x_pre(i,1)-obs_pos(1),y_pre(i,1)-obs_pos(2)];
                %计算距离车辆质心最近的障碍物点
                P = obs_pos + R_obs*d / norm(d);
                J_obst(i,1) = 1/(((x_pre(i,1)-P(1))^2+(y_pre(i,1)-P(2))^2)+0.000001);
            elseif i<=plan_Nc 
                omega(i,1)=target(i);
                theta_pre(i,1)=wrapToPi(theta_pre(i-1,1)+omega(i)*dt);
                x_pre(i,1)=x_pre(i-1,1)+v*cos(theta_pre(i,1))*dt;
                y_pre(i,1)=y_pre(i-1,1)+v*sin(theta_pre(i,1))*dt;
                d = [x_pre(i,1)-obs_pos(1),y_pre(i,1)-obs_pos(2)];
                P = obs_pos + R_obs*d / norm(d);
                J_obst(i,1) = 1/(((x_pre(i,1)-P(1))^2+(y_pre(i,1)-P(2))^2)+0.000001);
            else
                omega(i,1)=target(plan_Nc);
                theta_pre(i,1)=wrapToPi(theta_pre(i-1,1)+omega(i)*dt);
                x_pre(i,1)=x_pre(i-1,1)+v*cos(theta_pre(i,1))*dt;
                y_pre(i,1)=y_pre(i-1,1)+v*sin(theta_pre(i,1))*dt;

                d = [x_pre(i,1)-obs_pos(1),y_pre(i,1)-obs_pos(2)];
                P = obs_pos + R_obs*d / norm(d);
                J_obst(i,1) = 1/(((x_pre(i,1)-P(1))^2+(y_pre(i,1)-P(2))^2)+0.000001);
            end
        %J_obst=J_obst+1/(((X_predict(i,1))-X_obstacle(2,1))^2+(Y_predict(i,1)-Y_obstacle(2,1))^2+0.00001);
            THETA_error(i,1) = wrapToPi(theta_pre(i,1) - THETAref(i,1));
            Y_error(i,1)=y_pre(i,1)-Yref(i,1);%注意此处Yref与Y_ref要区分开来，Yref是局部期望路径,Y_ref为全局期望路径
            X_error(i,1)=x_pre(i,1)-Xref(i,1);

            e_lat_plan(i,1) = ...
            -sin(THETAref(i,1)) * X_error(i,1) ...
            +cos(THETAref(i,1)) * Y_error(i,1);           
    end 
    terminal_lat = e_lat_plan(plan_Np)^2;
    terminal_theta = THETA_error(plan_Np)^2;
    
    cost = ...
          e_lat_plan' * Q_lat * e_lat_plan ...
        + THETA_error' * Q_theta * THETA_error ...
        + omega(1:plan_Nc)' * R * omega(1:plan_Nc) ...
        + Q_terminal * terminal_lat ...
        + Q_terminal_theta * terminal_theta+S*sum(J_obst(:));
end