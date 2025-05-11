%% ��Ⱥ�㷨������ʱ�䴰�ĳ���·������(VRPTW)
clear;
clc;
close all;

%% �����������
n = 50;                 % �ͻ�����
m = 20;                  % ��������
capacity = 100;          % ��������
depot = [0, 0];         % ������������
F = 300;                 % �����̶�ʹ�óɱ����ĵ��е�F��
C = 2;                  % ��λ������ʻ�ɱ����ĵ��е�C��
mu = 0.001;               % ��λʱ�����ϵ�����ĵ��еĦ̣�
P = 10000;                 % ��λ���̼�ֵ���ĵ��е�P��
G = 0.5;                  % ��ʪ���豸��λʱ��ɱ����ĵ��е�G��
late_penalty = 2;       % ���ͷ�ϵ�����ĵ��н����гͷ���
speed = 30;              % ������ʻ�ٶ�

%% ��������ͻ�����
rng(0); % ����������ӣ���֤���������
customers.x = rand(n, 1) * 100;      % �ͻ�x����
customers.y = rand(n, 1) * 100;      % �ͻ�y����
customers.demand = randi([1, 10], n, 1); % �ͻ�����
customers.ready_time = randi([0, 50], n, 1); % �ͻ��������ʱ��
customers.due_time = customers.ready_time + randi([20, 50], n, 1); % �ͻ��������ʱ��
customers.service_time = randi([5, 10], n, 1); % �ͻ�����ʱ��

%% ��Ⱥ�㷨��������
num_ants = 100;          % ��������
max_iterations = 300;   % ����������
alpha = 1;              % ��Ϣ����Ҫ�̶�����
beta = 5;               % ����ʽ����
rho = 0.8;              % ��Ϣ�ػӷ�ϵ��
Q = 100;                % ��Ϣ�ظ���ǿ��

%% ����������
distance = zeros(n+1, n+1); % ������������
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
pheromone = ones(n+1, n+1); % ��ʼ��Ϣ��
for i = 1:n+1
    pheromone(i, i) = 0; % �����������Ϣ��Ϊ0
end

%% ��ʼ����ѽ�
best_cost = inf;
best_routes = [];
best_iteration = 0;
iter_best_cost = zeros(max_iterations, 1); % ��¼ÿ�ε��������Ž�

%% ��Ⱥ�㷨��ѭ��
for iter = 1:max_iterations
    % ��ʼ���������ϵ�·���ͳɱ�
    all_routes = cell(num_ants, 1);
    all_costs = zeros(num_ants, 1);
    
    % ÿֻ���Ϲ���·��
    for ant = 1:num_ants
        % ��ʼ��
        unvisited = 2:n+1; % δ���ʵĿͻ���2��n+1��1���������ģ�
        routes = cell(m, 1); % m������·��
        vehicle_load = zeros(m, 1); % ÿ�����ĵ�ǰ����
        vehicle_time = zeros(m, 1); % ÿ�����ĵ�ǰʱ��
        route_index = ones(m, 1); % ÿ������ǰ·��������
        
        % Ϊÿ��������·��
        for k = 1:m
            current_node = 1; % ���������Ŀ�ʼ
            routes{k}(route_index(k)) = current_node;
            
            % ������δ���ʵĿͻ��ҳ���δ����ʱ����
            while ~isempty(unvisited) && vehicle_load(k) < capacity
                % ����ת�Ƹ���
                probabilities = zeros(size(unvisited));
                for i = 1:length(unvisited)
                    j = unvisited(i);
                    if vehicle_load(k) + customers.demand(j-1) <= capacity
                        % ����ʱ����ز���
                        travel_time = distance(current_node, j) / speed;
                        arrival_time = vehicle_time(k) + travel_time;
                        wait_time = max(0, customers.ready_time(j-1) - arrival_time);
                        departure_time = arrival_time + wait_time + customers.service_time(j-1);
                        late_penalty_cost = late_penalty * max(0, departure_time - customers.due_time(j-1));
                        
                        % ��������ʽ��Ϣ
                        eta = 1 / (distance(current_node, j) + late_penalty_cost + 0.001);
                        
                        % ����ת�Ƹ���
                        probabilities(i) = pheromone(current_node, j)^alpha * eta^beta;
                    end
                end
                
                % ���û�п��еĿͻ���������ǰ������·��
                if sum(probabilities) == 0
                    break;
                end
                
                % ��һ������
                probabilities = probabilities / sum(probabilities);
                
                % ���̶�ѡ����һ���ͻ�
                r = rand;
                cumulative = 0;
                next_index = 1;
                for i = 1:length(probabilities)
                    cumulative = cumulative + probabilities(i);
                    if r <= cumulative
                        next_index = i;
                        break;
                    end
                end
                next_node = unvisited(next_index);
                
                % ���³���״̬
                vehicle_load(k) = vehicle_load(k) + customers.demand(next_node-1);
                travel_time = distance(current_node, next_node) / speed;
                arrival_time = vehicle_time(k) + travel_time;
                wait_time = max(0, customers.ready_time(next_node-1) - arrival_time);
                vehicle_time(k) = arrival_time + wait_time + customers.service_time(next_node-1);
                
                % ����·��
                route_index(k) = route_index(k) + 1;
                routes{k}(route_index(k)) = next_node;
                
                % ��δ���ʿͻ����Ƴ�
                unvisited(next_index) = [];
                
                % ���µ�ǰ�ڵ�
                current_node = next_node;
            end
            
            % ����������������
            route_index(k) = route_index(k) + 1;
            routes{k}(route_index(k)) = 1;
        end
        
        % �����ܳɱ���ʹ���޸ĺ��Ŀ�꺯����
        total_cost = 0;
        vehicle_used = 0;       % ��¼ʹ�õĳ�����
        for k = 1:m
            route = routes{k};
            if length(route) > 2 % �������������ĺ�һ���ͻ�
                vehicle_used = vehicle_used + 1;
                vehicle_time = 0;
                vehicle_distance = 0;
                damage_cost = 0;
                equipment_cost = 0;
                
                for i = 1:length(route)-1
                    current_node = route(i);
                    next_node = route(i+1);
                    
                    % ��ʻ�ɱ�
                    vehicle_distance = vehicle_distance + distance(current_node, next_node);
                    
                    % ����ʱ�䣨���ڻ�����豸�ɱ���
                    travel_time = distance(current_node, next_node) / speed;
                    equipment_cost = equipment_cost + travel_time * G;
                    
                    if next_node ~= 1 % ���������Ľڵ㣨�ͻ��ڵ㣩
                        j = next_node - 1; % �ͻ�����
                        damage_cost = damage_cost + travel_time * mu * P * customers.demand(j);
                        
                        % ʱ�䴰�ͷ��������ͷ���
                        vehicle_time = vehicle_time + travel_time;
                        arrival_time = vehicle_time;
                        wait_time = max(0, customers.ready_time(j) - arrival_time);
                        departure_time = arrival_time + wait_time + customers.service_time(j);
                        late_cost = late_penalty * max(0, departure_time - customers.due_time(j));
                        total_cost = total_cost + late_cost;
                        
                        vehicle_time = departure_time;
                    else
                        vehicle_time = vehicle_time + travel_time; % �����������ĵ�ʱ��
                    end
                end
                
                % �ۼӹ̶��ɱ�����ʻ�ɱ�
                total_cost = total_cost + vehicle_distance * C + F;
                total_cost = total_cost + damage_cost + equipment_cost;
            end
        end
        
        % ���浱ǰ���ϵ�·���ͳɱ�
        all_routes{ant} = routes;
        all_costs(ant) = total_cost;
        
        % ����ȫ�����Ž�
        if total_cost < best_cost
            best_cost = total_cost;
            best_routes = routes;
            best_iteration = iter;
        end
    end
    
    % ��¼��ǰ���������Ž�
    [min_cost, ~] = min(all_costs);
    iter_best_cost(iter) = min_cost;
    
    % ��Ϣ�ظ���
    pheromone = (1 - rho) * pheromone; % ��Ϣ�ػӷ�
    
    % ÿֻ���ϸ�����Ϣ��
    for ant = 1:num_ants
        routes = all_routes{ant};
        for k = 1:m
            route = routes{k};
            for i = 1:length(route)-1
                pheromone(route(i), route(i+1)) = pheromone(route(i), route(i+1)) + Q / all_costs(ant);
            end
        end
    end
    
    % ��ʾ����
    if mod(iter, 10) == 0
        fprintf('��������: %d, ���ųɱ�: %.2f\n', iter, best_cost);
    end
end

%% ������
fprintf('���Ž��ڵ� %d �ε����ҵ�\n', best_iteration);
fprintf('�ܳɱ�: %.2f\n', best_cost);

% ��ʾ����·��
fprintf('����·��:\n');
for k = 1:m
    route = best_routes{k};
    if length(route) > 2 % �������������ĺ�һ���ͻ�
        fprintf('���� %d: ', k);
        for i = 1:length(route)
            if route(i) == 1
                fprintf('��������');
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

% �ֽ���ʾ����ɱ�
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

fprintf('��ʻ�ɱ� (�����C): %.2f\n', total_distance * C);
fprintf('�����̶��ɱ�: %.2f\n', total_fixed_cost);
fprintf('���ͷ��ɱ�: %.2f\n', total_late_penalty);
fprintf('����ɱ�: %.2f\n', total_damage_cost);
fprintf('��ʪ���豸�ɱ�: %.2f\n', total_equipment_cost);
fprintf('�ܳɱ�: %.2f\n', total_distance * C + total_fixed_cost + total_late_penalty + total_damage_cost + total_equipment_cost);

%% ���ӻ�
% ��������ͼ
figure;
plot(1:max_iterations, iter_best_cost, 'LineWidth', 2);
title('��Ⱥ�㷨��������');
xlabel('��������');
ylabel('�ܳɱ�');
grid on;

% ����·��ͼ
figure;
hold on;
% ������������
plot(depot(1), depot(2), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
text(depot(1)+1, depot(2)+1, '��������');

% ���ƿͻ���
for i = 1:n
    plot(customers.x(i), customers.y(i), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    text(customers.x(i)+1, customers.y(i)+1, num2str(i));
end

% ����·��
colors = lines(m);
for k = 1:m
    route = best_routes{k};
    if length(route) > 2 % �������������ĺ�һ���ͻ�
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

title('���ų���·������');
xlabel('X����');
ylabel('Y����');
grid on;
axis equal;
legend('��������', '�ͻ�', '����·��');
hold off;    