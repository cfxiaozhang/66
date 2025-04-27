%% 辅助函数：计算路径成本
function total_cost = calculate_route_cost(routes, distance, travel_time, customer_demand, customer_time_window, penalty_coeff, F, C, G, mu, p)
total_cost = 0;
for v_idx = 1:length(routes)
    path = routes{v_idx};
    if isempty(path)
        continue;
    end
    
    % 车辆固定成本
    total_cost = total_cost + F;
    
    % 行驶成本
    for i = 1:length(path) - 1
        i_node = path(i);
        j_node = path(i + 1);
        total_cost = total_cost + C * distance(i_node, j_node);
    end
    
    % 时间窗惩罚成本
    current_time = 0;
    for i = 2:length(path) - 1 % 跳过起点和终点（配送中心）
        customer_idx = path(i) - 1;
        service_time = travel_time(path(i - 1), path(i));
        current_time = current_time + service_time;
        if current_time > customer_time_window(customer_idx, 2)
            delay = current_time - customer_time_window(customer_idx, 2);
            total_cost = total_cost + penalty_coeff(customer_idx) * delay;
        end
    end
    
    % 货损成本和温湿度设备成本
    for i = 1:length(path) - 1
        i_node = path(i);
        j_node = path(i + 1);
        travel_time_ij = travel_time(i_node, j_node);
        % 货损成本（假设每个客户点的需求在运输中产生损耗）
        if i_node > 1
            demand = customer_demand(i_node - 1);
            total_cost = total_cost + p * demand * mu * travel_time_ij;
        end
        % 温湿度设备成本
        total_cost = total_cost + G * travel_time_ij;
    end
end
end    