%% 1. 数据初始化
clear; clc;

% 基本参数
n_customers = 50;       % 客户点数量
Q = 50;                 % 车辆载重上限（单位：吨）
v = 40;                 % 车辆速度（公里/小时）
C = 2;                  % 单位距离行驶成本（元/公里）
F = 100;                % 车辆固定使用成本（元/辆）
mu = 0.01;              % 单位时间货损系数
G = 5;                  % 温湿度设备单位时间成本（元/小时）
alpha = 1;              % 信息素启发因子
beta = 5;               % 期望启发因子
rho = 0.8;              % 信息素挥发因子
m = round(2 * n_customers); % 蚂蚁数量
max_iter = 100;         % 最大迭代次数
Q_aco = 100;            % 信息素增强系数

% 随机生成数据
rng(0); % 固定随机种子便于复现
depot = [0, 0];        % 配送中心坐标
customer_coord = rand(n_customers, 2) * 100;  % 客户坐标（0 - 100公里）
customer_demand = rand(n_customers, 1) * 10;  % 客户需求（0 - 10吨）
customer_time_window = [rand(n_customers, 1) * 2, rand(n_customers, 1) * 8 + 2]; % 时间窗[2,10]小时

% 构建距离矩阵和时间矩阵
distance = zeros(n_customers + 1, n_customers + 1); % 0为配送中心
travel_time = zeros(n_customers + 1, n_customers + 1);
for i = 0:n_customers
    for j = 0:n_customers
        if i ~= j
            if i == 0
                coord_i = depot;
            else
                coord_i = customer_coord(i, :);
            end
            if j == 0
                coord_j = depot;
            else
                coord_j = customer_coord(j, :);
            end
            distance(i + 1, j + 1) = norm(coord_i - coord_j);
            travel_time(i + 1, j + 1) = distance(i + 1, j + 1) / v;
        end
    end
end

% 其他参数
p = 100;                % 卷烟单位价值（元/吨）
penalty_coeff = rand(n_customers, 1) * 5; % 延迟惩罚系数（随机生成）

%% 2. 蚁群算法参数初始化
nodes = 1:n_customers + 1; % 节点编号，1为配送中心，2~n+1为客户
tau = ones(n_customers + 1, n_customers + 1); % 初始信息素矩阵
best_cost = inf;
best_routes = [];

%% 3. 迭代过程
cost_history = [];
for iter = 1:max_iter
    routes = cell(m, 1); % 存储每只蚂蚁的路径
    costs = zeros(m, 1);
    
    for k = 1:m
        tabu = {}; % 每辆车的路径集合
        load = []; % 车辆载重
        time = []; % 车辆当前时间
        current_node = 1; % 从配送中心出发（节点1）
        vehicle_idx = 1;
        tabu{vehicle_idx} = [1]; % 初始位置为配送中心
        load(vehicle_idx) = 0;
        time(vehicle_idx) = 0;
        
        while true
            allowed = setdiff(nodes, [tabu{:}]);
            allowed(allowed == 1) = []; % 配送中心只能作为起点/终点
            
            % 计算转移概率
            if isempty(allowed)
                % 返回配送中心
                tabu{vehicle_idx} = [tabu{vehicle_idx}, 1];
                break
            end
            
            prob = zeros(1, length(allowed));
            for j = 1:length(allowed)
                ij = [current_node, allowed(j)];
                tau_ij = tau(current_node, allowed(j));
                eta_ij = 1 / distance(current_node, allowed(j));
                numerator = tau_ij^alpha * eta_ij^beta;
                denominator = sum(tau(current_node, allowed).^alpha .* (1./distance(current_node, allowed)).^beta);
                prob(j) = numerator / denominator;
            end
            
            % 轮盘赌选择下一个节点（替代 randsample）
            r = rand;
            cum_prob = cumsum(prob);
            for j = 1:length(allowed)
                if r <= cum_prob(j)
                    next_node = allowed(j);
                    break;
                end
            end
            
            % 检查载重约束
            if load(vehicle_idx) + customer_demand(next_node - 1) > Q
                % 换车
                vehicle_idx = vehicle_idx + 1;
                tabu{vehicle_idx} = [1, next_node]; % 新车辆从配送中心出发
                load(vehicle_idx) = customer_demand(next_node - 1);
                time(vehicle_idx) = travel_time(current_node, next_node);
            else
                tabu{vehicle_idx} = [tabu{vehicle_idx}, next_node];
                load(vehicle_idx) = load(vehicle_idx) + customer_demand(next_node - 1);
                time(vehicle_idx) = time(vehicle_idx) + travel_time(current_node, next_node);
            end
            
            current_node = next_node;
        end
        
        % 转换为完整路径并计算成本
        full_routes = tabu;
        route_cost = calculate_route_cost(full_routes, distance, travel_time, customer_demand, customer_time_window, penalty_coeff, F, C, G, mu, p);
        routes{k} = full_routes;
        costs(k) = route_cost;
        
        % 更新最优解
        if route_cost < best_cost
            best_cost = route_cost;
            best_routes = full_routes;
        end
    end
    
    % 信息素更新
    delta_tau = zeros(n_customers + 1, n_customers + 1);
    for k = 1:m
        route = routes{k};
        for v_idx = 1:length(route)
            path = route{v_idx};
            for i = 1:length(path) - 1
                i_node = path(i);
                j_node = path(i + 1);
                delta_tau(i_node, j_node) = delta_tau(i_node, j_node) + Q_aco / costs(k);
            end
        end
    end
    tau = (1 - rho) * tau + delta_tau;
    
    cost_history = [cost_history; best_cost];
    disp(['迭代', num2str(iter), '，最优成本：', num2str(best_cost)]);
end

%% 4. 结果可视化
% 绘制迭代图
figure;
plot(1:max_iter, cost_history, 'b-', 'LineWidth', 1.5);
xlabel('迭代次数');
ylabel('总成本（元）');
title('蚁群算法迭代曲线');
grid on;

% 绘制路径图
figure;
scatter(depot(1), depot(2), 100, 'r', 'filled', 'MarkerEdgeColor', 'k');
text(depot(1) + 2, depot(2), '配送中心', 'FontSize', 10);
hold on;
for v_idx = 1:length(best_routes)
    path = best_routes{v_idx};
    coords = zeros(length(path), 2);
    for i = 1:length(path)
        if path(i) == 1
            coords(i, :) = depot;
        else
            coords(i, :) = customer_coord(path(i) - 1, :);
        end
    end
    plot(coords(:, 1), coords(:, 2), 'b-o', 'LineWidth', 1, 'MarkerSize', 6);
    for i = 1:length(path)
        if path(i) > 1
            text(customer_coord(path(i) - 1, 1) + 2, customer_coord(path(i) - 1, 2), ...
                 ['C', num2str(path(i) - 1)], 'FontSize', 8);
        end
    end
end
xlabel('X坐标（公里）');
ylabel('Y坐标（公里）');
title('最优配送路径图');
axis equal;
grid on;