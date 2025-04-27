%% ��������������·���ɱ�
function total_cost = calculate_route_cost(routes, distance, travel_time, customer_demand, customer_time_window, penalty_coeff, F, C, G, mu, p)
total_cost = 0;
for v_idx = 1:length(routes)
    path = routes{v_idx};
    if isempty(path)
        continue;
    end
    
    % �����̶��ɱ�
    total_cost = total_cost + F;
    
    % ��ʻ�ɱ�
    for i = 1:length(path) - 1
        i_node = path(i);
        j_node = path(i + 1);
        total_cost = total_cost + C * distance(i_node, j_node);
    end
    
    % ʱ�䴰�ͷ��ɱ�
    current_time = 0;
    for i = 2:length(path) - 1 % ���������յ㣨�������ģ�
        customer_idx = path(i) - 1;
        service_time = travel_time(path(i - 1), path(i));
        current_time = current_time + service_time;
        if current_time > customer_time_window(customer_idx, 2)
            delay = current_time - customer_time_window(customer_idx, 2);
            total_cost = total_cost + penalty_coeff(customer_idx) * delay;
        end
    end
    
    % ����ɱ�����ʪ���豸�ɱ�
    for i = 1:length(path) - 1
        i_node = path(i);
        j_node = path(i + 1);
        travel_time_ij = travel_time(i_node, j_node);
        % ����ɱ�������ÿ���ͻ���������������в�����ģ�
        if i_node > 1
            demand = customer_demand(i_node - 1);
            total_cost = total_cost + p * demand * mu * travel_time_ij;
        end
        % ��ʪ���豸�ɱ�
        total_cost = total_cost + G * travel_time_ij;
    end
end
end    