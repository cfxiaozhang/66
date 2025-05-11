%% �Ľ���Ⱥ�㷨����ʱ�䴰�ĳ���·������ (VRPTW)
clear;
clc;
close all;

%% ��������
n = 50;                 % �ͻ�����
m = 20;                 % ��������
capacity = 100;          % ��������
depot = [0, 0];         % �ֿ�����
F = 300;                 % �����̶�ʹ�óɱ�
C = 2;                  % ��λ��������ɱ�
mu = 0.001;               % ��λʱ�����ϵ��
P = 10000;                 % ��λ�����ֵ
G = 0.5;                  % �����豸��λʱ��ɱ�
late_penalty = 2;       % �ٵ��ͷ�ϵ��
speed = 1;              % �����ٶ�

%% ������ɿͻ�����
rng(0);
customers.x = rand(n, 1) * 100;
customers.y = rand(n, 1) * 100;
customers.demand = randi([1, 10], n, 1);
customers.ready_time = randi([0, 50], n, 1);
customers.due_time = customers.ready_time + randi([20, 50], n, 1);
customers.service_time = randi([5, 10], n, 1);

%% ��Ⱥ�㷨����
num_ants = 100;
max_iterations = 300;
alpha0 = 1;     % ��ʼ��Ϣ����Ҫ�̶�
beta0 = 5;      % ��ʼ����������Ҫ�̶�
rho0 = 0.8;     % ��ʼ��Ϣ�ػӷ�ϵ��
Q = 100;
elitist_weight = 10;  % ��Ӣ���ϼ�Ȩ����

tau_max_init = 1;     % ��ʼ�����Ϣ��
tau_min_init = 1e-4;  % ��ʼ��С��Ϣ��

alpha = alpha0;
beta = beta0;
rho = rho0;
tau_max = tau_max_init;
tau_min = tau_min_init;

%% �������
distance = zeros(n+1, n+1);
for i = 1:n+1
    for j = 1:n+1
        if i ~= j
            if i == 1
                x1 = depot(1);
                y1 = depot(2);
            else
                x1 = customers.x(i-1);
                y1 = customers.y(i-1);
            end
            if j == 1
                x2 = depot(1);
                y2 = depot(2);
            else
                x2 = customers.x(j-1);
                y2 = customers.y(j-1);
            end
            distance(i, j) = sqrt((x1 - x2)^2 + (y1 - y2)^2);
        end
    end
end

%% ��ʼ����Ϣ�ؾ���
pheromone = tau_max * ones(n+1, n+1);
for i = 1:n+1
    pheromone(i, i) = 0;
end

%% ��¼����
best_cost = inf;
best_routes = [];
best_iteration = 0;
iter_best_cost = zeros(max_iterations, 1);

%% ��ѭ��
for iter = 1:max_iterations
    % *** ����Ӧ�������� ***
    alpha = alpha0 + (3 - alpha0) * iter / max_iterations;  % alpha����
    beta  = beta0  - (beta0 - 1) * iter / max_iterations;   % beta����
    rho   = rho0   * (1 - iter / (2*max_iterations));        % �𲽼�Сrho
    
    % *** �����С��Ϣ������ ***
    if iter == 1
        tau_max = tau_max_init;
        tau_min = tau_min_init;
    else
        tau_max = 1 / (rho * best_cost);
        tau_min = tau_max / (2 * n);
    end

    all_routes = cell(num_ants, 1);
    all_costs = zeros(num_ants, 1);
    
    % �����������
    for ant = 1:num_ants
        unvisited = 2:n+1;
        routes = cell(m, 1);
        vehicle_load = zeros(m, 1);
        vehicle_time = zeros(m, 1);
        route_index = ones(m, 1);
        for k = 1:m
            current_node = 1;
            routes{k}(route_index(k)) = current_node;
            while ~isempty(unvisited) && vehicle_load(k) < capacity
                probabilities = zeros(size(unvisited));
                for i = 1:length(unvisited)
                    j = unvisited(i);
                    if vehicle_load(k) + customers.demand(j-1) <= capacity
                        travel_time = distance(current_node, j) / speed;
                        arrival_time = vehicle_time(k) + travel_time;
                        wait_time = max(0, customers.ready_time(j-1) - arrival_time);
                        departure_time = arrival_time + wait_time + customers.service_time(j-1);
                        late_penalty_cost = late_penalty * max(0, departure_time - customers.due_time(j-1));
                        eta = 1 / (distance(current_node, j) + late_penalty_cost + 0.001);
                        probabilities(i) = pheromone(current_node, j)^alpha * eta^beta;
                    end
                end
                if sum(probabilities) == 0
                    break;
                end
                probabilities = probabilities / sum(probabilities);
                r = rand;
                cumulative = 0;
                next_index = 1;
                for i = 1:length(probabilities)
                    cumulative = cumulative + probabilities(i);
                    if r <= cumulative
                        next_index = i; break;
                    end
                end
                next_node = unvisited(next_index);
                vehicle_load(k) = vehicle_load(k) + customers.demand(next_node-1);
                travel_time = distance(current_node, next_node) / speed;
                arrival_time = vehicle_time(k) + travel_time;
                wait_time = max(0, customers.ready_time(next_node-1) - arrival_time);
                vehicle_time(k) = arrival_time + wait_time + customers.service_time(next_node-1);
                route_index(k) = route_index(k) + 1;
                routes{k}(route_index(k)) = next_node;
                unvisited(next_index) = [];
                current_node = next_node;
            end
            route_index(k) = route_index(k) + 1;
            routes{k}(route_index(k)) = 1;
        end
        
        % *** ����������(2-opt)�ֲ����� ***
        for k = 1:m
            route = routes{k};
            if length(route) > 3
                routes{k} = two_opt(route, distance);
            end
        end
        
        % �����ܴ���
        total_cost = 0;
        for k = 1:m
            route = routes{k};
            if length(route) > 2
                vehicle_time = 0;
                vehicle_distance = 0;
                damage_cost = 0;
                equipment_cost = 0;
                for i = 1:length(route)-1
                    current_node = route(i);
                    next_node = route(i+1);
                    vehicle_distance = vehicle_distance + distance(current_node, next_node);
                    travel_time = distance(current_node, next_node) / speed;
                    equipment_cost = equipment_cost + travel_time * G;
                    if next_node ~= 1
                        j = next_node - 1;
                        damage_cost = damage_cost + travel_time * mu * P * customers.demand(j);
                        vehicle_time = vehicle_time + travel_time;
                        arrival_time = vehicle_time;
                        wait_time = max(0, customers.ready_time(j) - arrival_time);
                        departure_time = arrival_time + wait_time + customers.service_time(j);
                        late_cost = late_penalty * max(0, departure_time - customers.due_time(j));
                        total_cost = total_cost + late_cost;
                        vehicle_time = departure_time;
                    else
                        vehicle_time = vehicle_time + travel_time;
                    end
                end
                total_cost = total_cost + vehicle_distance * C + F;
                total_cost = total_cost + damage_cost + equipment_cost;
            end
        end
        
        all_routes{ant} = routes;
        all_costs(ant) = total_cost;
        if total_cost < best_cost
            best_cost = total_cost;
            best_routes = routes;
            best_iteration = iter;
        end
    end
    
    % ��¼��������
    [min_cost, idx] = min(all_costs);
    iter_best_cost(iter) = min_cost;

    % ��Ϣ�ػӷ�
    pheromone = (1 - rho) * pheromone;
    
    % ����·������
    for ant = 1:num_ants
        routes = all_routes{ant};
        for k = 1:m
            route = routes{k};
            for i = 1:length(route)-1
                pheromone(route(i), route(i+1)) = pheromone(route(i), route(i+1)) + Q / all_costs(ant);
            end
        end
    end
    % *** ��Ӣ����ϵͳ������ʷ����·������ǿ����Ϣ�� ***
    for k = 1:m
        route = best_routes{k};
        for i = 1:length(route)-1
            pheromone(route(i), route(i+1)) = pheromone(route(i), route(i+1)) + elitist_weight * Q / best_cost;
        end
    end
    % *** ���-��С��Ϣ�����ƣ�MMAS�� ***
    pheromone(pheromone > tau_max) = tau_max;
    pheromone(pheromone < tau_min) = tau_min;

    % ��ʾ����
    if mod(iter, 10) == 0
        fprintf('��������: %d, ��ǰ�����ܴ���: %.2f\n', iter, best_cost);
    end
end

%% ������
fprintf('���Ž�����ڵ� %d ����\n', best_iteration);
fprintf('�ܴ���: %.2f\n', best_cost);

fprintf('·������:\n');
for k = 1:m
    route = best_routes{k};
    if length(route) > 2
        fprintf('���� %d: ', k);
        for i = 1:length(route)
            if route(i) == 1
                fprintf('�ֿ�');
            else
                fprintf('�ͻ� %d', route(i)-1);
            end
            if i < length(route)
                fprintf(' -> ');
            end
        end
        fprintf('\n');
    end
end

% ����ɱ�ͳ��
fprintf('\n�ɱ��ֽ�:\n');
total_distance = 0;
total_fixed_cost = 0;
total_late_penalty = 0;
total_damage_cost = 0;
total_equipment_cost = 0;
for k = 1:m
    route = best_routes{k};
    if length(route) > 2
        vehicle_time = 0;
        vehicle_distance = 0;
        damage_cost = 0;
        equipment_cost = 0;
        late_cost = 0;
        for i = 1:length(route)-1
            current_node = route(i);
            next_node = route(i+1);
            vehicle_distance = vehicle_distance + distance(current_node, next_node);
            travel_time = distance(current_node, next_node) / speed;
            equipment_cost = equipment_cost + travel_time * G;
            if next_node ~= 1
                j = next_node - 1;
                damage_cost = damage_cost + travel_time * mu * P * customers.demand(j);
                vehicle_time = vehicle_time + travel_time;
                arrival_time = vehicle_time;
                wait_time = max(0, customers.ready_time(j) - arrival_time);
                departure_time = arrival_time + wait_time + customers.service_time(j);
                late_cost = late_cost + late_penalty * max(0, departure_time - customers.due_time(j));
                vehicle_time = departure_time;
            else
                vehicle_time = vehicle_time + travel_time;
            end
        end
        total_distance = total_distance + vehicle_distance;
        total_fixed_cost = total_fixed_cost + F;
        total_late_penalty = total_late_penalty + late_cost;
        total_damage_cost = total_damage_cost + damage_cost;
        total_equipment_cost = total_equipment_cost + equipment_cost;
    end
end
fprintf('����ɱ�(����*C): %.2f\n', total_distance * C);
fprintf('�����̶��ɱ�: %.2f\n', total_fixed_cost);
fprintf('�ٵ��ͷ�: %.2f\n', total_late_penalty);
fprintf('��ĳɱ�: %.2f\n', total_damage_cost);
fprintf('�����豸�ɱ�: %.2f\n', total_equipment_cost);
fprintf('�ܴ���: %.2f\n', total_distance * C + total_fixed_cost + total_late_penalty + total_damage_cost + total_equipment_cost);

%% ��������
figure;
plot(1:max_iterations, iter_best_cost, 'LineWidth', 2);
title('��Ⱥ�㷨��������');
xlabel('��������');
ylabel('�ܴ���');
grid on;

%% ·�����ӻ�
figure; hold on;
plot(depot(1), depot(2), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
text(depot(1)+1, depot(2)+1, '�ֿ�');
for i = 1:n
    plot(customers.x(i), customers.y(i), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    text(customers.x(i)+1, customers.y(i)+1, num2str(i));
end
colors = lines(m);
for k = 1:m
    route = best_routes{k};
    if length(route) > 2
        x = zeros(length(route), 1);
        y = zeros(length(route), 1);
        for i = 1:length(route)
            if route(i) == 1
                x(i) = depot(1);
                y(i) = depot(2);
            else
                x(i) = customers.x(route(i)-1);
                y(i) = customers.y(route(i)-1);
            end
        end
        plot(x, y, 'Color', colors(k,:), 'LineWidth', 1.5);
        text(mean(x), mean(y), ['���� ', num2str(k)], 'Color', colors(k,:));
    end
end
title('���ų���·��');
xlabel('X����');
ylabel('Y����');
grid on;
axis equal;
legend('�ֿ�', '�ͻ�', '����·��');
hold off;