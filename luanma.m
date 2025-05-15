% main_cvrptw_aco.m

% --- 0. 清理工作区并设置随机种子 ---
clear; clc; close all;
rng('default'); % 为了结果可复现

% --- 1. 问题定义 ---
% 节点1是仓库 (Depot)。节点2到N+1是客户。
% Coords = [x_depot, y_depot; x_cust1, y_cust1; ...];
% Demand = [0; demand_cust1; ...]; % 仓库需求为0
% ServiceTime = [0; service_cust1; ...]; % 仓库服务时间为0
% TimeWindow = [ES_depot, LS_depot; ES_cust1, LS_cust1; ...]; % [最早开始服务时间, 最晚开始服务时间]

% 示例数据 (请替换为您的实际数据)
num_customers = 20; % 客户数量
depot_coord = [50 50]; % 仓库坐标
customer_coords = rand(num_customers, 2) * 100; % 客户坐标 (0-100范围)

Nodes.Coords = [depot_coord; customer_coords];
Nodes.Demand = [0; randi([10, 30], num_customers, 1)]; % 客户需求
Nodes.ServiceTime = [0; randi([5, 10], num_customers, 1)]; % 客户服务时长

% 时间窗设置: ES_i = 随机, LS_i = ES_i + 随机裕量
Nodes.TimeWindow = [0, 1000; % 仓库的时间窗 (通常很大)
                   randi([0, 200], num_customers, 1), zeros(num_customers, 1)]; 
for i = 1:num_customers
    Nodes.TimeWindow(i+1, 2) = Nodes.TimeWindow(i+1, 1) + randi([50, 150]); % LS = ES + 裕量
end

Vehicle.Capacity = 150; % 车辆容量
Vehicle.Speed = 1;     % 车辆速度 (距离单位/时间单位)

% 成本参数
Cost.PerDistance = 1.0;     % 每单位距离的成本
Cost.PerWaitTime = 0.5;     % 每单位等待时间的成本
Cost.PerTardinessTime = 2.0;% 每单位晚到时间的惩罚成本 (服务开始时间晚于LS)

num_nodes = size(Nodes.Coords, 1); % 总节点数 (仓库 + 客户)

% --- 2. ACO 参数 ---
ACO.NumAnts = 30;        % 蚂蚁数量
ACO.MaxIter = 150;       % 最大迭代次数
ACO.Alpha = 1;           % 信息素启发因子 (Pheromone influence)
ACO.Beta = 3;            % 期望启发因子 (Heuristic influence)
ACO.Rho = 0.1;           % 信息素蒸发率 (Pheromone evaporation rate)
ACO.Q = 100;             % 信息素强度常数 (Pheromone deposit constant)
ACO.Tau0 = 0.1;          % 初始信息素水平 (Initial pheromone level)

% --- 3. 初始化 ---
% 计算距离矩阵
DistMatrix = zeros(num_nodes, num_nodes);
for i = 1:num_nodes
    for j = i+1:num_nodes
        DistMatrix(i,j) = norm(Nodes.Coords(i,:) - Nodes.Coords(j,:));
        DistMatrix(j,i) = DistMatrix(i,j);
    end
end

% 初始化信息素矩阵
Pheromone = ones(num_nodes, num_nodes) * ACO.Tau0;

% 计算启发信息矩阵 (通常是距离的倒数)
Heuristic = zeros(num_nodes, num_nodes);
for i = 1:num_nodes
    for j = 1:num_nodes
        if DistMatrix(i,j) > 0
            Heuristic(i,j) = 1 / DistMatrix(i,j); 
        end
    end
end

% 存储最优解
BestSol.Routes = {};
BestSol.Cost = inf;
BestSol.CostHistory = zeros(ACO.MaxIter, 1); % 记录每次迭代的最优成本

% --- 4. 主ACO循环 ---
fprintf('开始ACO求解CVRPTW...\n');
for iter = 1:ACO.MaxIter
    AntSolutions = cell(ACO.NumAnts, 1); % 存储每只蚂蚁的解
    AntCosts = zeros(ACO.NumAnts, 1);    % 存储每只蚂蚁解的总成本

    for k = 1:ACO.NumAnts % 对每只蚂蚁
        CurrentAnt.Routes = {};    % 当前蚂蚁构建的路径集合
        CurrentAnt.TotalCost = 0;  % 当前蚂蚁的总成本
        
        CustomersToVisit = true(1, num_customers); % 标记客户是否已被访问 (true表示未访问)
        num_unvisited = num_customers; % 未访问客户数量

        while num_unvisited > 0 % 当仍有客户未被服务时
            % 为当前蚂蚁开始一条新路径 (一辆新车)
            current_route = [1]; % 从仓库出发 (节点索引1)
            current_load = 0;    % 当前车辆载重
            current_time = Nodes.TimeWindow(1,1); % 车辆从仓库出发的时间 (通常是仓库的ES, e.g., 0)
            current_route_cost = 0; % 当前路径的成本

            % 构建当前路径
            while true 
                last_node_idx = current_route(end); % 当前路径中的最后一个节点
                
                possible_next_nodes_indices = []; % 存储候选下一个节点的实际索引
                selection_probabilities = [];     % 存储选择各候选节点的概率

                % 遍历所有客户，寻找可行的下一个客户
                for cust_idx = 1:num_customers % cust_idx 是客户的编号 (1 到 num_customers)
                    actual_node_idx = cust_idx + 1; % actual_node_idx 是客户在Nodes结构中的索引

                    if CustomersToVisit(cust_idx) % 如果客户尚未被访问
                        % 1. 检查容量约束
                        if current_load + Nodes.Demand(actual_node_idx) <= Vehicle.Capacity
                            % 2. 计算选择该节点的启发信息和信息素浓度
                            % (软时间窗意味着总是"可行"，但成本会变化，此处选择概率不直接惩罚，成本在之后计算)
                            prob_value = (Pheromone(last_node_idx, actual_node_idx) ^ ACO.Alpha) * ...
                                         (Heuristic(last_node_idx, actual_node_idx) ^ ACO.Beta);
                            
                            if prob_value == 0 || isnan(prob_value) || isinf(prob_value)
                                prob_value = 1e-9; % 防止概率为0或无效值
                            end
                            
                            possible_next_nodes_indices(end+1) = actual_node_idx;
                            selection_probabilities(end+1) = prob_value;
                        end
                    end
                end

                if isempty(possible_next_nodes_indices)
                    break; % 没有可选择的下一个客户 (或所有客户已服务完毕)，结束当前路径构建
                end

                % 轮盘赌选择下一个客户
                selection_probabilities = selection_probabilities / sum(selection_probabilities);
                if any(isnan(selection_probabilities)) || any(isinf(selection_probabilities)) || sum(selection_probabilities)==0
                    rand_idx = randi(length(possible_next_nodes_indices)); % 随机选一个
                else
                    r = rand();
                    cum_probs = cumsum(selection_probabilities);
                    chosen_idx_in_possible = find(r <= cum_probs, 1, 'first');
                    if isempty(chosen_idx_in_possible) % 以防万一
                         rand_idx = randi(length(possible_next_nodes_indices));
                         chosen_idx_in_possible = rand_idx;
                    end
                end
                chosen_node_actual_idx = possible_next_nodes_indices(chosen_idx_in_possible);
                
                % --- 将选中的客户加入当前路径，并更新路径状态 ---
                % a. 计算行驶成本和时间
                travel_dist = DistMatrix(last_node_idx, chosen_node_actual_idx);
                travel_time = travel_dist / Vehicle.Speed;
                current_route_cost = current_route_cost + travel_dist * Cost.PerDistance;
                
                arrival_at_chosen = current_time + travel_time; % 到达选中客户的时间
                
                % b. 计算服务时间、等待成本、晚到成本
                service_start_chosen = max(arrival_at_chosen, Nodes.TimeWindow(chosen_node_actual_idx, 1));
                
                wait_time_chosen = service_start_chosen - arrival_at_chosen;
                current_route_cost = current_route_cost + wait_time_chosen * Cost.PerWaitTime;
                
                tardiness_chosen = max(0, service_start_chosen - Nodes.TimeWindow(chosen_node_actual_idx, 2));
                current_route_cost = current_route_cost + tardiness_chosen * Cost.PerTardinessTime;
                
                departure_time_chosen = service_start_chosen + Nodes.ServiceTime(chosen_node_actual_idx);
                
                % c. 更新路径状态
                current_route(end+1) = chosen_node_actual_idx; % 添加到路径
                current_load = current_load + Nodes.Demand(chosen_node_actual_idx); % 更新载重
                current_time = departure_time_chosen; % 更新当前时间为离开此客户的时间
                
                % d. 标记客户为已访问
                CustomersToVisit(chosen_node_actual_idx - 1) = false; % -1 是因为CustomersToVisit索引对应客户编号
                num_unvisited = num_unvisited - 1;
            end % 结束当前路径的客户添加 (inner while true)
            
            % 当前路径构建完毕 (没有更多客户可服务 或 车辆容量不足以服务任何剩余客户)
            % 如果路径服务了至少一个客户，则返回仓库
            if length(current_route) > 1 
                last_customer_in_route_idx = current_route(end);
                
                % a. 计算返回仓库的行驶成本和时间
                travel_dist_to_depot = DistMatrix(last_customer_in_route_idx, 1);
                travel_time_to_depot = travel_dist_to_depot / Vehicle.Speed;
                current_route_cost = current_route_cost + travel_dist_to_depot * Cost.PerDistance;
                
                arrival_at_depot = current_time + travel_time_to_depot;
                
                % b. 检查返回仓库是否晚点 (如果仓库有最晚返回时间LS)
                depot_tardiness = max(0, arrival_at_depot - Nodes.TimeWindow(1,2));
                current_route_cost = current_route_cost + depot_tardiness * Cost.PerTardinessTime;
                
                current_route(end+1) = 1; % 路径末尾添加仓库
                
                CurrentAnt.Routes{end+1} = current_route; % 将此路径存入蚂蚁的解中
                CurrentAnt.TotalCost = CurrentAnt.TotalCost + current_route_cost; %累加总成本
            elseif num_unvisited > 0 && isempty(CurrentAnt.Routes) && length(current_route) == 1
                % 特殊情况: 蚂蚁无法服务任何客户 (可能是参数或问题定义问题)
                % 对于软时间窗，这应该很少发生，除非需求远超容量或时间窗极度不合理
                CurrentAnt.TotalCost = inf; % 给予极高成本
                break; % 退出当前蚂蚁的路径构建
            end

        end % 结束当前蚂蚁的所有路径构建 (while num_unvisited > 0)
        
        AntSolutions{k} = CurrentAnt; % 存储当前蚂蚁的完整解
        AntCosts(k) = CurrentAnt.TotalCost; % 存储当前蚂蚁解的总成本

        % 更新全局最优解
        if CurrentAnt.TotalCost < BestSol.Cost
            BestSol.Cost = CurrentAnt.TotalCost;
            BestSol.Routes = CurrentAnt.Routes;
            fprintf('迭代 %d, 蚂蚁 %d: 新最优成本 = %.2f\n', iter, k, BestSol.Cost);
        end
    end % 结束蚂蚁循环

    % --- 信息素更新 ---
    % 1. 信息素蒸发
    Pheromone = (1 - ACO.Rho) * Pheromone;

    % 2. 信息素沉积
    for k_ant = 1:ACO.NumAnts
        if isinf(AntCosts(k_ant)); continue; end % 如果蚂蚁解无效，则跳过
        
        ant_solution_routes = AntSolutions{k_ant}.Routes;
        delta_tau = ACO.Q / AntCosts(k_ant); % Q / Lk
        
        for r = 1:length(ant_solution_routes) % 遍历蚂蚁的每条路径
            route = ant_solution_routes{r};
            for i = 1:length(route)-1 % 遍历路径上的每条边
                node_from = route(i);
                node_to = route(i+1);
                Pheromone(node_from, node_to) = Pheromone(node_from, node_to) + delta_tau;
                Pheromone(node_to, node_from) = Pheromone(node_to, node_from) + delta_tau; % 对称问题
            end
        end
    end
    
    % 确保信息素不低于某个阈值，防止过早停滞
    Pheromone(Pheromone < ACO.Tau0/10) = ACO.Tau0/10; % 例如，初始信息素的1/10

    BestSol.CostHistory(iter) = BestSol.Cost; % 记录本次迭代的最优成本
    if mod(iter, 10) == 0 || iter == 1 || iter == ACO.MaxIter
      fprintf('迭代 %d: 当前最优成本 = %.2f\n', iter, BestSol.Cost);
    end
end % 结束迭代循环

fprintf('ACO求解完成。\n找到的最优成本: %.2f\n', BestSol.Cost);

% --- 5. 结果与可视化 ---

% 迭代图
figure;
plot(1:ACO.MaxIter, BestSol.CostHistory, 'LineWidth', 2);
xlabel('迭代次数');
ylabel('总成本');
title('ACO收敛曲线 (CVRPTW)');
grid on;

% 路径图
figure;
hold on;
% 绘制仓库
plot(Nodes.Coords(1,1), Nodes.Coords(1,2), 'ks', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'DisplayName', '仓库');
% 绘制客户点
plot(Nodes.Coords(2:end,1), Nodes.Coords(2:end,2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', '客户');

% 为每条路径使用不同颜色
route_colors = lines(length(BestSol.Routes)); 
for r = 1:length(BestSol.Routes)
    route_nodes_indices = BestSol.Routes{r}; % 获取路径中的节点索引
    route_coords = Nodes.Coords(route_nodes_indices, :); % 获取这些节点的坐标
    plot(route_coords(:,1), route_coords(:,2), 'LineWidth', 1.5, 'Color', route_colors(r,:), 'Marker', 'none');
    
    % 为客户节点添加文本标签 (客户编号, 需求, 时间窗)
    for i = 1:length(route_nodes_indices)
        node_idx = route_nodes_indices(i);
        if node_idx > 1 % 如果是客户节点
            customer_label_num = node_idx - 1;
            text_str = sprintf('C%d (D:%d, TW:[%d,%d])', ...
                               customer_label_num, ...
                               Nodes.Demand(node_idx), ...
                               Nodes.TimeWindow(node_idx,1), ...
                               Nodes.TimeWindow(node_idx,2));
            text(Nodes.Coords(node_idx,1)+0.8, Nodes.Coords(node_idx,2)+0.8, text_str, 'FontSize', 7);
        end
    end
end

title(['找到的最优路径方案 (总成本: ' num2str(BestSol.Cost, '%.2f') ')']);
xlabel('X 坐标');
ylabel('Y 坐标');
legend('show', 'Location', 'bestoutside');
axis equal;
grid on;
hold off;

% 在控制台输出最优路径
disp('找到的最优路径:');
for r = 1:length(BestSol.Routes)
    % 将节点索引映射回客户编号 (0 代表仓库)
    route_display = BestSol.Routes{r} - 1; 
    fprintf('路径 %d: %s\n', r, mat2str(route_display));
end

% (可选) 详细验证最优解的成本计算
TotalVerifiedCost = 0;
disp('最优解的详细成本分解:');
for r_idx = 1:length(BestSol.Routes)
    route = BestSol.Routes{r_idx};
    route_cost_detail = 0;
    current_time_verify = Nodes.TimeWindow(1,1); % 从仓库出发时间
    current_load_verify = 0;
    
    fprintf('\n路径 %d: 仓库 -> ', r_idx);
    
    for i = 1:(length(route)-1) % 遍历路径的每一段
        from_node = route(i);
        to_node = route(i+1);
        
        % 1. 行驶成本
        dist = DistMatrix(from_node, to_node);
        travel_c = dist * Cost.PerDistance;
        route_cost_detail = route_cost_detail + travel_c;
        
        travel_t = dist / Vehicle.Speed;
        arrival_at_to = current_time_verify + travel_t;
        
        if to_node ~= 1 % 如果目标是客户
            fprintf('C%d -> ', to_node-1);
            
            % 2. 等待成本
            service_start = max(arrival_at_to, Nodes.TimeWindow(to_node,1));
            wait_t = service_start - arrival_at_to;
            wait_c = wait_t * Cost.PerWaitTime;
            route_cost_detail = route_cost_detail + wait_c;
            
            % 3. 晚到惩罚成本
            tardiness_t = max(0, service_start - Nodes.TimeWindow(to_node,2));
            tardiness_c = tardiness_t * Cost.PerTardinessTime;
            route_cost_detail = route_cost_detail + tardiness_c;
            
            % 更新时间和负载
            departure_time = service_start + Nodes.ServiceTime(to_node);
            current_time_verify = departure_time;
            current_load_verify = current_load_verify + Nodes.Demand(to_node);

            fprintf('(到达:%.1f, 等待:%.1f, 服务开始:%.1f, 晚到:%.1f, 离开:%.1f, 载重:%d) - ', ...
                arrival_at_to, wait_t, service_start, tardiness_t, departure_time, current_load_verify);

        else % 如果目标是仓库 (路径结束)
            fprintf('仓库 ');
            current_time_verify = arrival_at_to; % 到达仓库时间
            % 检查返回仓库是否晚点
            depot_tardiness_t = max(0, current_time_verify - Nodes.TimeWindow(1,2));
            depot_tardiness_c = depot_tardiness_t * Cost.PerTardinessTime;
            route_cost_detail = route_cost_detail + depot_tardiness_c;
            fprintf('(到达仓库:%.1f, 仓库晚到:%.1f) ', current_time_verify, depot_tardiness_t);
        end
    end
    fprintf('\n路径 %d 成本: %.2f. 最终载重: %d (车辆容量: %d)\n', r_idx, route_cost_detail, current_load_verify, Vehicle.Capacity);
    TotalVerifiedCost = TotalVerifiedCost + route_cost_detail;
end
fprintf('\n总验证成本 (来自详细分解): %.2f\n', TotalVerifiedCost);
fprintf('ACO找到的最优成本: %.2f\n', BestSol.Cost);
if abs(TotalVerifiedCost - BestSol.Cost) < 1e-3
    fprintf('成本验证成功。\n');
else
    fprintf('警告: 成本验证失败。差额: %.4f\n', abs(TotalVerifiedCost - BestSol.Cost));
end