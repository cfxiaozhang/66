%% 1. ���ݳ�ʼ��
clear; clc;

% ��������
n_customers = 20;       % �ͻ�������
Q = 50;                 % �����������ޣ���λ���֣�
v = 40;                 % �����ٶȣ�����/Сʱ��
C = 2;                  % ��λ������ʻ�ɱ���Ԫ/���
F = 100;                % �����̶�ʹ�óɱ���Ԫ/����
mu = 0.01;              % ��λʱ�����ϵ��
G = 5;                  % ��ʪ���豸��λʱ��ɱ���Ԫ/Сʱ��
alpha = 1;              % ��Ϣ����������
beta = 5;               % ������������
rho = 0.8;              % ��Ϣ�ػӷ�����
m = round(2 * n_customers); % �����������ĵ�����1.5 - 2.5���ͻ�����
max_iter = 300;         % ����������
Q_aco = 100;            % ��Ϣ����ǿϵ��

% �����������
rng(0); % �̶�������ӱ��ڸ���
depot = [0, 0];        % ������������
customer_coord = rand(n_customers, 2) * 100;  % �ͻ����꣨0 - 100���
customer_demand = rand(n_customers, 1) * 10;  % �ͻ�����0 - 10�֣�
customer_time_window = [rand(n_customers, 1) * 2, rand(n_customers, 1) * 8 + 2]; % ʱ�䴰[2,10]Сʱ

% ������������ʱ�����
distance = zeros(n_customers + 1, n_customers + 1); % 0Ϊ��������
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

% ��������
p = 100;                % ���̵�λ��ֵ��Ԫ/�֣�
penalty_coeff = rand(n_customers, 1) * 5; % �ӳٳͷ�ϵ����������ɣ�

%% 2. ��Ⱥ�㷨������ʼ��
nodes = 1:n_customers + 1; % �ڵ��ţ�1Ϊ�������ģ�2~n+1Ϊ�ͻ�
tau = ones(n_customers + 1, n_customers + 1); % ��ʼ��Ϣ�ؾ���
best_cost = inf;
best_routes = [];

%% 3. ��������
cost_history = [];
for iter = 1:max_iter
    routes = cell(m, 1); % �洢ÿֻ���ϵ�·��
    costs = zeros(m, 1);
    
    for k = 1:m
        tabu = {}; % ÿ������·������
        load = []; % ��������
        time = []; % ������ǰʱ��
        current_node = 1; % ���������ĳ������ڵ�1��
        vehicle_idx = 1;
        tabu{vehicle_idx} = [1]; % ��ʼλ��Ϊ��������
        load(vehicle_idx) = 0;
        time(vehicle_idx) = 0;
        
        while true
            allowed = setdiff(nodes, [tabu{:}]);
            allowed(allowed == 1) = []; % ��������ֻ����Ϊ���/�յ�
            
            % ����ת�Ƹ���
            if isempty(allowed)
                % ������������
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
            
            % ���̶�ѡ����һ���ڵ㣨��� randsample��
            r = rand;
            cum_prob = cumsum(prob);
            for j = 1:length(allowed)
                if r <= cum_prob(j)
                    next_node = allowed(j);
                    break;
                end
            end
            
            % �������Լ��
            if load(vehicle_idx) + customer_demand(next_node - 1) > Q
                % ����
                vehicle_idx = vehicle_idx + 1;
                tabu{vehicle_idx} = [1, next_node]; % �³������������ĳ���
                load(vehicle_idx) = customer_demand(next_node - 1);
                time(vehicle_idx) = travel_time(current_node, next_node);
            else
                tabu{vehicle_idx} = [tabu{vehicle_idx}, next_node];
                load(vehicle_idx) = load(vehicle_idx) + customer_demand(next_node - 1);
                time(vehicle_idx) = time(vehicle_idx) + travel_time(current_node, next_node);
            end
            
            current_node = next_node;
        end
        
        % ת��Ϊ����·��������ɱ�
        full_routes = tabu;
        route_cost = calculate_route_cost(full_routes, distance, travel_time, customer_demand, customer_time_window, penalty_coeff, F, C, G, mu, p);
        routes{k} = full_routes;
        costs(k) = route_cost;
        
        % �������Ž�
        if route_cost < best_cost
            best_cost = route_cost;
            best_routes = full_routes;
        end
    end
    
    % ��Ϣ�ظ���
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
    disp(['����', num2str(iter), '�����ųɱ���', num2str(best_cost)]);
end

%% 4. ������ӻ�
% ���Ƶ���ͼ
figure;
plot(1:max_iter, cost_history, 'b-', 'LineWidth', 1.5);
xlabel('��������');
ylabel('�ܳɱ���Ԫ��');
title('��Ⱥ�㷨��������');
grid on;

% ����·��ͼ
figure;
scatter(depot(1), depot(2), 100, 'r', 'filled', 'MarkerEdgeColor', 'k');
text(depot(1) + 2, depot(2), '��������', 'FontSize', 10);
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
xlabel('X���꣨���');
ylabel('Y���꣨���');
title('��������·��ͼ');
axis equal;
grid on;

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