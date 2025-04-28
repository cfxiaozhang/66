clc;
clear;

%% 参数设置
numCustomers = 50; % 客户数量
maxIter = 50; % 最大迭代次数
numAnts = 100; % 蚂蚁数量
alpha = 1; % 信息素重要性因子
beta = 5; % 启发式因子重要性
rho = 0.8; % 信息素挥发系数
Q = 200; % 信息素总量

% 成本系数
F = 300; % 每辆车的固定成本
C = 2; % 每单位距离的运输成本
penalty_coeff = randi([5, 15], numCustomers, 1); % 每个客户的时间窗惩罚系数
p = 100; % 货损成本系数
mu = 1; % 每单位时间的货损成本
G = 0.5; % 温湿度设备成本系数

% 客户点和时间窗数据
depot = [0, 0]; % 配送中心坐标
customerPos = randi([0, 100], numCustomers, 2); % 随机生成客户点
timeWindows = sort(randi([1, 100], numCustomers, 2), 2); % 时间窗 [start, end]
demand = randi([1, 10], numCustomers, 1); % 每个客户的需求
vehicleCapacity = 50; % 车辆最大载量

% 计算距离矩阵和时间矩阵
numNodes = numCustomers + 1; % 节点总数（含配送中心）
nodes = [depot; customerPos]; % 所有节点坐标
distMatrix = squareform(pdist(nodes)); % 距离矩阵
travelTime = distMatrix / 30; % 假设平均速度为30单位距离/时间

%% 初始化
pheromone = ones(numNodes); % 信息素矩阵
bestRoute = [];
bestTotalCost = inf;
costHistory = zeros(maxIter, 1); % 用于保存每次迭代的最优总成本

%% 蚁群算法主循环
for iter = 1:maxIter
    routes = cell(numAnts, 1);
    totalCosts = zeros(numAnts, 1);
    
    % 每只蚂蚁构造解
    for k = 1:numAnts
        route = [1]; % 从配送中心开始
        load = 0; % 当前载量
        time = 0; % 当前时间
        totalCost = F; % 初始总成本包含车辆固定成本
        
        while length(route) < numCustomers + 1
            currentNode = route(end);
            unvisited = setdiff(2:numNodes, route); % 未访问节点
            
            if isempty(unvisited)
                route = [route, 1]; % 返回配送中心
                break;
            end
            
            % 计算转移概率
            prob = zeros(size(unvisited));
            for j = 1:length(unvisited)
                nextNode = unvisited(j);
                tau = pheromone(currentNode, nextNode)^alpha; % 信息素
                eta = (1 / distMatrix(currentNode, nextNode))^beta; % 启发式信息
                prob(j) = tau * eta;
            end
            prob = prob / sum(prob);
            
            % 按概率选择下一个节点
            nextNode = unvisited(rouletteWheelSelection(prob));
            
            % 判断是否满足载量与时间窗约束
            travelTime_ij = travelTime(currentNode, nextNode);
            arrivalTime = time + travelTime_ij; % 到达下一个节点的时间
            timeViolation = max(0, arrivalTime - timeWindows(nextNode - 1, 2)); % 时间窗约束违规
            
            if load + demand(nextNode - 1) <= vehicleCapacity
                route = [route, nextNode];
                load = load + demand(nextNode - 1);
                time = max(arrivalTime, timeWindows(nextNode - 1, 1)); % 更新时间，考虑时间窗约束
                
                % 累加行驶成本
                totalCost = totalCost + C * distMatrix(currentNode, nextNode);
                % 累加时间窗惩罚成本
                totalCost = totalCost + penalty_coeff(nextNode - 1) * timeViolation;
                % 累加货损成本
                totalCost = totalCost + p * demand(nextNode - 1) * mu * travelTime_ij;
                % 累加温湿度设备成本
                totalCost = totalCost + G * travelTime_ij;
            else
                route = [route, 1]; % 返回配送中心
                load = 0;
                time = 0;
            end
        end
        
        routes{k} = route;
        totalCosts(k) = totalCost;
    end
    
    % 更新全局最优解
    [minTotalCost, minIndex] = min(totalCosts);
    if minTotalCost < bestTotalCost
        bestTotalCost = minTotalCost;
        bestRoute = routes{minIndex};
    end
    
    % 保存当前迭代的最优总成本
    costHistory(iter) = bestTotalCost;

    % 打印当前迭代的总成本到命令行窗口
    fprintf('迭代 %d: 最优总成本 = %.2f\n', iter, bestTotalCost);

    % 更新信息素
    pheromone = (1 - rho) * pheromone;
    for k = 1:numAnts
        route = routes{k};
        for j = 1:length(route) - 1
            pheromone(route(j), route(j + 1)) = pheromone(route(j), route(j + 1)) + Q / totalCosts(k);
        end
    end
end

%% 绘制迭代图（折线图）
figure(1);
plot(1:maxIter, costHistory, '-o', 'LineWidth', 1.5);
title('迭代过程');
xlabel('迭代次数');
ylabel('总成本');
grid on;

%% 绘制最优路径
figure(2);
plot(nodes(:, 1), nodes(:, 2), 'bo');
hold on;
plot(nodes(1, 1), nodes(1, 2), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k'); % 配送中心
for i = 2:numCustomers + 1
    text(nodes(i, 1), nodes(i, 2), [' ', num2str(i - 1)]);
end
for i = 1:length(bestRoute) - 1
    plot([nodes(bestRoute(i), 1), nodes(bestRoute(i + 1), 1)], ...
         [nodes(bestRoute(i), 2), nodes(bestRoute(i + 1), 2)], 'r-');
end
title('最优路径');
xlabel('X坐标');
ylabel('Y坐标');
grid on;