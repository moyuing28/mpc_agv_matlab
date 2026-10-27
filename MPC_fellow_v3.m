%这是运用了时变预测矩阵的误差状态控制的mpc控制
%参考轨迹为正方形轨迹,本项目将设法去应对正方形的角度跳变问题。
%% mpc控制基本原理
% 车辆运动学模型：
% [x,y,theta]' = [cos(theta),sin(theta),0]'*V_target + [0,0,1]'*oobs_existmega;
% 泰勒展开线性化：
% A = [0,0,-V_target*sin(theta);
%     0,0,V_target*cos(theta);
%     0,0,0;];
%B  = [cos(theta),0;
%     sin(theta),0;
%     0,1;];
% kesi = [x-x_r;y-y_r;theta-theta_r];
% u = [v-V_target;omega-omega_r];
% kesi_dot = A*Chi + B*u;
% 欧拉近似离散化：
% A = [1,0,-V_target*sin(theta);
%     0,1,V_target*cos(theta);
%     0,0,1;];
%B  = [dt*cos(theta),0;
%     dt*sin(theta),0;
%     0,dt;];
% kesi_dot = A*kesi + B*u;
% 转化为增量模型：
% Xi = {kesi;u};%u是上一时间步的控制量
%增量控制矩阵：
% A_tilde = [A,B;
%     O,E]
% B_tilde = [B;
%     E]
% C_tilde = [0 0 0 1 0;
%             0 0 0 0 1;];
% 生成预测状态矩阵：
% 由A_tilde,B_tilde,C_tilde生成PHI和THETA

% 奖励函数的设计

% 约束的设计
clear; close all; clc;

%% 1. 共有参数设置
L = 0.5;              % 轮距 (m)
dt = 0.01;             % 采样时间 (s)
T_max = 20;           % 最大模拟时长 (s)
V_target = 0.8;       % 期望巡航速度 (m/s)
v_ref = V_target;
% LookAhead_Dist = 0.2; % 路径点切换阈值 (m) (判断是否到达航点的距离)(mpc不需要）
%在正方形的轨迹里面还是需要的，但是新的mpc代码没有用到它，作者（就是我）用了0.02提升了采样精度
Np = 60;              %预测步长
Nc = 30;               %控制步长
%增广矩阵参数
Nx = 3;               %输出的增广维度参数
Nu = 2;               %输入的增广维度列数
% 初始位姿 (与前面算法统一)
% x0 = 0; y0 = -0.75; theta0 = 0;%改动了初始位置，因为在没进入参考轨迹的情况下mpc控制发挥不好
x0 = 0; y0 = 0; theta0 = 0;
%% mpc控制的设计
% 初始状态设定
x = x0; y = y0; theta = theta0;
%增量矩阵的填充部分
O = zeros(Nu,Nx);
E = eye(Nu);
% 参考轨迹设计？
%正方形轨迹：
squa_len = 5;%正方形边长
Waypoints = [0, 0;
             squa_len, 0;
             squa_len, squa_len;
             0, squa_len;
             0, 0];
theta_ref_value =[0;pi/2;pi;-pi/2];
omega_ref = 0;%正方形轨迹的参考角速度，总为0
draw = true;%绘图标志符，决定代码是否绘图
arri = true;%判断AGV是否到达轨迹起点，可以忽略，当起点处没有参考轨迹时mpc控制的效果是不好的，直接让mpc从车辆起点开始
t_arri = 0;%AGV到达轨迹起点的时间点
the_end = false;%参考的理想小车是否完成运动的标志
comparison = false;
current_idx = 2;%正方形状态机标签，正方形的参考轨迹建立其实就是一个状态机了
Row=10;%松弛因子（可用可不用）
pre_t = (Np+1)*dt;%控制的预测时长

v = 0;
omega = 0;

history.x=[];
history.y=[];
history.theta=[];
history.v=[];
history.omega=[];

history_ref.x=[];
history_ref.y=[];
d_k = {};

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
%% 参考轨迹生成+求取状态误差量
for t = 0:dt:T_max
    for i = t:dt:(t+pre_t)
        U = [v;omega];
        err1 = [x-Waypoints(1,1),y-Waypoints(1,2)];
        distance1 = norm(err1);
            %% 预留代码段，判断mpc在无参考轨迹下达到指定轨迹用
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
            %% 真正的参考轨迹生成
            elseif  arri == true
                %参考轨迹测试代码
                % switch current_idx
                %     case  2
                %         x_ref = V_target*(t-t_arri)-5*(current_idx-2);
                %         y_ref = Waypoints(current_idx,2);
                %         err = [x_ref-Waypoints(current_idx,1),y_ref-Waypoints(current_idx,2)];
                %         err_distance = norm(err);
                %         history_ref.x(num) = x_ref;
                %         history_ref.y(num) = y_ref;
                %         num = num+1;
                %         if err_distance<=0.02
                %             current_idx = current_idx+1;
                %         end
                %     case 3
                %         x_ref = Waypoints(current_idx,1);
                %         y_ref = V_target*(t-t_arri)-5*(current_idx-2);
                %         err = [x_ref-Waypoints(current_idx,1),y_ref-Waypoints(current_idx,2)];
                %         err_distance = norm(err);
                %         history_ref.x(num) = x_ref;
                %         history_ref.y(num) = y_ref;
                %         num = num+1;
                %         if err_distance<=0.02
                %             current_idx = current_idx+1;
                %         end
                %      case 4
                %         x_ref = 5-(V_target*(t-t_arri)-5*(current_idx-2));
                %         y_ref = Waypoints(current_idx,2);
                %         err = [x_ref-Waypoints(current_idx,1),y_ref-Waypoints(current_idx,2)];
                %         err_distance = norm(err);
                %         history_ref.x(num) = x_ref;
                %         history_ref.y(num) = y_ref;
                %         num = num+1;
                %         if err_distance<=0.02
                %             current_idx = current_idx+1;
                %         end
                %     otherwise
                %         if the_end ==false
                %             x_ref = Waypoints(current_idx,1);
                %             y_ref = 5-abs(V_target*(t-t_arri)-5*(current_idx-2));
                %             err = [x_ref-Waypoints(current_idx,1),y_ref-Waypoints(current_idx,2)];
                %             err_distance = norm(err);
                %             history_ref.x(num) = x_ref;
                %             history_ref.y(num) = y_ref;
                %             num = num+1;
                %             if err_distance<=0.02
                %              the_end =true;
                %             end
                %         else
                %             x_ref = Waypoints(current_idx,1);
                %             y_ref = Waypoints(current_idx,2);
                %             history_ref.x(num) = x_ref;
                %             history_ref.y(num) = y_ref;
                %             num = num+1;
                %         end
                % 
                % end
                
                 switch current_idx
                    case  2
                        % 参考位置设立
                        x_pre_ref = V_target*(i-t_arri)-squa_len*(current_idx-2);
                        y_pre_ref = Waypoints(current_idx,2);        
                        % 参考路线转变标志判断，计算目前参考点到转换点的距离来判断是否进行转换
                        ref_err = [x_pre_ref-Waypoints(current_idx,1),y_pre_ref-Waypoints(current_idx,2)];
                        ref_distance = norm(ref_err);
                        %预测序列的生成
                        history_pre_ref.x(pre_num) = x_pre_ref;
                        history_pre_ref.y(pre_num) = y_pre_ref;
                        history_pre_ref.theta(pre_num) = theta_ref_value(current_idx-1);
                        if ref_distance<=0.02
                            current_idx = current_idx+1;
                            % omega_ref =10;
                        end
                        %状态机标签的存储，预测序列的生成过程可能会跨状态机，所以生成完预测序列后要回溯状态机的标签
                        history_pre_ref.current_idx(pre_num) = current_idx;
                        pre_num = pre_num+1;
                    case 3
                        x_pre_ref = Waypoints(current_idx,1);
                        y_pre_ref = V_target*(i-t_arri)-squa_len*(current_idx-2);
                         if y_pre_ref <= 0
                            y_pre_ref = 0;
                        end
                        ref_err = [x_pre_ref-Waypoints(current_idx,1),y_pre_ref-Waypoints(current_idx,2)];
                        ref_distance = norm(ref_err);
                        history_pre_ref.x(pre_num) = x_pre_ref;
                        history_pre_ref.y(pre_num) = y_pre_ref;
                        history_pre_ref.theta(pre_num) = theta_ref_value(current_idx-1);               
                        if ref_distance<=0.02
                            current_idx = current_idx+1;
                            % omega_ref =10;
                        end
                        history_pre_ref.current_idx(pre_num) = current_idx;
                        pre_num = pre_num+1;
                     case 4
                        x_pre_ref = squa_len-(V_target*(i-t_arri)-squa_len*(current_idx-2));
                        if x_pre_ref >=5
                            x_pre_ref = 5;
                        end
                        y_pre_ref = Waypoints(current_idx,2);         
                        ref_err = [x_pre_ref-Waypoints(current_idx,1),y_pre_ref-Waypoints(current_idx,2)];
                        ref_distance = norm(ref_err);
                        history_pre_ref.x(pre_num) = x_pre_ref;
                        history_pre_ref.y(pre_num) = y_pre_ref;
                        history_pre_ref.theta(pre_num) = theta_ref_value(current_idx-1); 
                        if ref_distance<=0.02
                            current_idx = current_idx+1;
                            % omega_ref =10;
                        end
                        history_pre_ref.current_idx(pre_num) = current_idx;
                        pre_num = pre_num+1;
                    otherwise%最后一段比较特殊，如果仿真时间较长，参考点在仿真结束之前会一直停在结束点处
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
% %存储第一个弯道的预测序列,用于测试异常
% if 5.5 <=t && t<= 7
%     history.pre_test{test_num,1} = history_pre_ref.theta;
%     test_num = test_num+1;
% end
pre_num = 1;%重置预测序列的索引
current_idx = history_pre_ref.current_idx(1);%重置轨迹状态机的标签
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
%预测序列的第一个数据即是当前的参考值
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
% R=0.5*eye(Nu*Nc);
R_single = diag([0.5, 0.05]);   % Δv 权重 0.5, Δω 权重 0.05
R = kron(eye(Nc), R_single);
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
    theta_r_i =history_pre_ref.theta(i);
    A_cell{i} = [1, 0, -dt*V_target*sin(theta_r_i);
                 0, 1,  dt*V_target*cos(theta_r_i);
                 0, 0, 1];
    B_cell{i} = [dt*cos(theta_r_i), 0;
                 dt*sin(theta_r_i), 0;
                 0,               dt];
end

% for i = 1:Np
%     theta_r_i =history_pre_ref.theta(i);
%     A_cell{i} = [1, 0, -dt*V_target*sin(theta_r_i);
%                  0, 1,  dt*V_target*cos(theta_r_i);
%                  0, 0, 1];
%     B_cell{i} = [dt*cos(theta_r_i), 0;
%                  dt*sin(theta_r_i), 0;
%                  0,               dt];
% end

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
%% 线性时不变mpc设计参考
% %这个用参考的状态量更好？
% A = [1,0,-dt*V_target*sin(theta);
%     0,1,dt*V_target*cos(theta);
%     0,0,1;];
% B  = [dt*cos(theta),0;
%     dt*sin(theta),0;
%     0,dt;];
% % A = [1,0,-dt*V_target*sin(theta_ref);
% %     0,1,dt*V_target*cos(theta_ref);
% %     0,0,1;];
% % B  = [dt*cos(theta_ref),0;
% %     dt*sin(theta_ref),0;
% %     0,dt;];
% %增量控制矩阵：
% A_tilde = [A,B;
%     O,E];
% B_tilde = [B;
%     E];
% C_tilde =[ 1 0 0 0 0;
%         0 1 0 0 0;
%         0 0 1 0 0];
% %预测状态矩阵：
% PHI_cell=cell(Np,1);
% THETA_cell=cell(Np,Nc); 
%     for i = 1:Np
%         PHI_cell{i,1} = C_tilde * A_tilde^i;
%     end
% % disp(puthy);
%     for i =1:Np
%         for j = 1:Nc
%             if j<=i
%              THETA_cell{i,j} = C_tilde*A_tilde^(i-j)*B_tilde;
%             else
%              THETA_cell{i,j} = zeros(Nx,Nu);
%             end
%         end
%     end
% PHI=cell2mat(PHI_cell);%size=[Nx*Np Nx+Nu]
% THETA=cell2mat(THETA_cell);%size=[Nx*Np Nu*Nc]
% %%
% % disp(xita);
% 

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
% error=PHI*kesi;...
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
% umin = [-1.2; -13.0];       % 放宽角速度下限，缩窄速度下限
% umax = [1.2;  13.0];        % 放宽两个速度上限，缩窄速度上限
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
% 
% % dd_v = X(1);
% % dd_omega = X(2);
% % d_v = kesi(4);
% % d_omega = kesi(5);
% % v = d_v + dd_v + V_target;
% % omega = d_omega + dd_omega + omega_ref;
% 
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
% theta = wrapToPi(theta + w_act * dt);
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



%% 1. 共有参数设置
L = 0.5;              % 轮距 (m)
% dt = 0.01;             % 采样时间 (s)
T_max = 50;           % 最大模拟时长 (s)
V_target = 0.8;       % 期望巡航速度 (m/s)
LookAhead_Dist = 0.2; % 路径点切换阈值 (m) (判断是否到达航点的距离)

% --- 共用的底层航向控制 PID 参数 (用于跟踪目标角度) ---
params_head.Kp = 4.0;
params_head.Ki = 0.0;
params_head.Kd = 0.5;
params_head.MaxOutput = 12.0; % 最大角速度限制 (rad/s)

%% 2. 算法特有参数设置
% --- [算法 A] 矢量合成 PID 参数 ---
params_vec.Kp = 3.0;
params_vec.Ki = 0.0;
params_vec.Kd = 0.8;
params_vec.MaxOutput = 1.0;  % 权重最大限制为1

% --- [算法 B] 纯PID跟踪参数 (计算偏角) ---
params_pid_lat.Kp = 3.0;
params_pid_lat.Ki = 0.0;
params_pid_lat.Kd = 0.8;
params_pid_lat.MaxOutput = pi/2; % 最大纠偏角度限制为 90 度 

% --- [算法 C] Pure Pursuit 参数 ---
Ld_pp = 0.4; % 前瞻距离 Lookahead Distance (m)

% --- [算法 D] Stanley 参数 ---
k_stanley = 1.5; % 增益系数，控制对横向误差的敏感度

%% 3. 定义路径点 (Waypoints)
Waypoints = [0, 0;
             squa_len, 0;
             squa_len, squa_len;
             0, squa_len;
             0, 0];
[num_points, ~] = size(Waypoints);

% 初始位姿 (所有算法统一)
% x0 = 0; y0 = 0; theta0 = 0; 

%% ================== 仿真1：矢量合成算法 ==================
state_vec.integral_sum = 0; state_vec.error_prev = 0;
state_head_A.integral_sum = 0; state_head_A.error_prev = 0;
x = x0; y = y0; theta = theta0; current_idx = 2;          
history_vec.x = []; history_vec.y = []; history_vec.error = []; history_vec.pro = [];
history_pid.v=[]; history_pid.omega=[];history_pid.t =[];
for t = 0:dt:T_max
    Last_Point = Waypoints(current_idx-1, :); Target_Point = Waypoints(current_idx, :);
    if norm([x, y] - Target_Point) < LookAhead_Dist
        if current_idx < num_points, current_idx = current_idx + 1;
            Last_Point = Waypoints(current_idx-1, :); Target_Point = Waypoints(current_idx, :);
        else, break; end
    end
    
    dx = Target_Point(1) - Last_Point(1); dy = Target_Point(2) - Last_Point(2);
    if abs(dx) < 1e-6, k = 1e4; path_angle = pi/2 * sign(dy); 
    else, k = dy / dx; path_angle = atan2(dy, dx); end
    b_line = Last_Point(2) - k * Last_Point(1);
    
    if abs(k) >= 1e4, error = abs(x - Last_Point(1));
    % 横向误差计算
    else, error = abs(k*x - y + b_line) / sqrt(k^2 + 1); end
    % 前馈方向向量
    T_vec = [cos(path_angle), sin(path_angle)];
    %二维叉积（标量值）,可以判定点是在
    cross_product = (x - Last_Point(1)) * dy - (y - Last_Point(2)) * dx;
    if cross_product > 0, norm_angle = path_angle + pi/2; % 左侧，向右纠
    else, norm_angle = path_angle - pi/2; end             % 右侧，向左纠
    %结论好像有问题，cross_product > 0点应该是在参考直线的右侧的，向左纠。另一边同理
    % 反馈纠正方向向量
    F_vec = [cos(norm_angle), sin(norm_angle)];
    
    [pid_out_vec, state_vec] = PID_Controller(error, 0, dt, params_vec, state_vec);
    Proportion = min(max(pid_out_vec, 0), 1); %论文里的beta
    history_vec.pro(end+1) = Proportion;
    E_vec_raw = Proportion * F_vec + (1 - Proportion) * T_vec;%最终的控制方向向量
    E_vec = E_vec_raw / norm(E_vec_raw);
    V_cmd_world = E_vec * V_target; 
    
    theta_cmd = atan2(V_cmd_world(2), V_cmd_world(1));%航向角，控制量
    angle_diff = wrapToPi(theta_cmd - theta);
    
    [omega, state_head_A] = PID_Controller(angle_diff, 0, dt, params_head, state_head_A);

    v = max(V_target * cos(angle_diff), 0); 
    % 它通过在航向误差较大时降低前进速度，
    % 让车辆有更多时间通过 omega 调整自身航向，从而减少高速前进过程中产生的横向偏差。
    vr = v + omega * L / 2; vl = v - omega * L / 2;
    v_act = (vr + vl) / 2; w_act = (vr - vl) / L;
    x = x + v_act * cos(theta) * dt; y = y + v_act * sin(theta) * dt; 
    theta = wrapToPi(theta + w_act * dt);
    history_pid.v(end+1) = v;
    history_pid.omega(end+1) = omega;
    history_pid.t(end+1) = t; 
    signed_error = sign(-cross_product) * error; if cross_product == 0, signed_error = error; end
    history_vec.x(end+1) = x; history_vec.y(end+1) = y; history_vec.error(end+1) = signed_error;
end
% if draw == true
% % 绘图：'-o' 表示实线+圆圈，'linewidth' 加粗线条，'markersize' 加大点
% plot(history_ref.x, history_ref.y, '-o', 'LineWidth', 1.5, 'MarkerSize', 2, 'MarkerFaceColor', 'b');
% 
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
    % 绘制第二组曲线（新轨迹），线宽同样调细
    plot(history_vec.x, history_vec.y, '-s', 'LineWidth', 0.8, ...
         'MarkerSize', 2, 'Color', [0.47, 0.67, 0.19], 'DisplayName', 'pid轨迹');  
    % 绘制第三组曲线（新轨迹），线宽同样调细
    plot(history.x, history.y, '-s', 'LineWidth', 0.8, ...
         'MarkerSize', 2, 'Color', [0.85, 0.33, 0.10], 'DisplayName', 'mpc轨迹1');
     

 
if comparison == true

    % 绘制第四组曲线（新轨迹），线宽同样调细
    plot(history.x1, history.y1, '-s', 'LineWidth', 0.8, ...
         'MarkerSize', 2, 'MarkerFaceColor', 'y', 'DisplayName', 'mpc轨迹2');

     % 绘制第五组曲线（新轨迹），线宽同样调细
    plot(history.x2, history.y2, '-s', 'LineWidth', 0.8, ...
         'MarkerSize', 2, 'MarkerFaceColor', 'c', 'DisplayName', 'mpc轨迹3');
end

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
title('\omega 对比');
legend('MPC', 'PID');

figure;
plot(history.t, history.v, 'b','DisplayName', 'mpc omega'); hold on;
plot(history_pid.t, history_pid.v, 'r','DisplayName', 'pid omega');
grid on;
xlabel('Time (s)');
ylabel('v (m/s)');
title('v 对比');
legend('MPC', 'PID');
% 绘制参考位置x，y的变化曲线。
figure;
plot(history.t, history_ref.x, 'b','DisplayName', '参考x轴'); hold on;
plot(history.t, history.theta, 'g'); hold on;
grid on;
xlabel('Time (s)');
ylabel('m');
title('参考位置x');
plot(history.t, history_ref.y, 'r','DisplayName', '参考y轴'); hold on;
grid on;
plot(history.t, history.x, 'c','DisplayName', '实际x轴'); hold on;
grid on;
plot(history.t, history.y, 'y','DisplayName', '实际y轴'); hold on;
grid on;
xlabel('Time (s)');
ylabel('m');
title('参考位置y');


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

