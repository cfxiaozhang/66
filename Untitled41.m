clc;
clear;

%% ��������
numCustomers = 50; % �ͻ�����
maxIter = 50; % ����������
numAnts = 100; % ��������
alpha = 1; % ��Ϣ������ʽ����
beta = 5; % ��������ʽ����
rho = 0.8; % ��Ϣ�ػӷ�����
Q = 100; % ��Ϣ����ǿϵ��

% �ɱ�ϵ��
F = 500; % ÿ�����Ĺ̶��ɱ�
C = 2; % ÿ��������ɱ�
penalty_coeff = randi([5, 15], numCustomers, 1); % ÿ���ͻ�ʱ�䴰��ϵ��
mu = 1; % ÿ��λʱ��Ļ���ɱ�
G = 0.5; % ��ʪ���豸�ɱ�ϵ��

% �ͻ�λ�ú�������������
depot = [0, 0]; % ������������
customerPos = randi([0, 100], numCustomers, 2); % ������ɿͻ�����
timeWindows = sort(randi([1, 100], numCustomers, 2), 2); % ʱ�䴰 [start, end]
demand = randi([1, 10], numCustomers, 1); % ÿ���ͻ�������
vehicleCapacity = 50; % �����������

% ������ʱ�����
nodes = [depot; customerPos]; % ���нڵ�����
distMatrix = squareform(pdist(nodes)); % �ڵ��������
travelTime = distMatrix / 30; % ƽ���ٶȼ���Ϊ30��λ/Сʱ

% ʹ��������㷨��ó�ʼ��
numNodes = numCustomers + 1;
bestRoute = nearestNeighborInit(numNodes, distMatrix);
bestTotalCost = inf;

%% ��Ⱥ�㷨����
pheromone = ones(numNodes); % ��Ϣ�ؾ���
costHistory = zeros(maxIter, 1);

for iter = 1:maxIter
    routes = cell(numAnts, 1);
    totalCosts = zeros(numAnts, 1);
    
    for k = 1:numAnts
        route = [1]; % ���������ĳ���
        load = 0; % ��ǰ����
        time = 0; % ��ǰʱ��
        totalCost = F; % ��ʼ�ܳɱ������̶��ɱ�
        
        while length(route) < numCustomers + 1
            currentNode = route(end);
            unvisited = setdiff(2:numNodes, route); % δ���ʽڵ�
            
            if isempty(unvisited)
                route = [route, 1];
                break;
            end
            
            % ת�Ƹ��ʼ���
            prob = zeros(size(unvisited));
            for j = 1:length(unvisited)
                nextNode = unvisited(j);
                tau = pheromone(currentNode, nextNode)^alpha; % ��Ϣ��
                eta = (1 / distMatrix(currentNode, nextNode))^beta; % ����ʽ��Ϣ
                prob(j) = tau * eta;
            end
            prob = prob / sum(prob);
            
            % ���̶�ѡ����һ�ڵ�
            nextNode = unvisited(rouletteWheelSelection(prob));
            
            % �ж��Ƿ�����ʱ�䴰Լ��������Լ��
            travelTime_ij = travelTime(currentNode, nextNode);
            arrivalTime = time + travelTime_ij; % ������һ���ڵ��ʱ��
            timeViolation = max(0, arrivalTime - timeWindows(nextNode - 1, 2));
            
            if load + demand(nextNode - 1) <= vehicleCapacity
                route = [route, nextNode];
                load = load + demand(nextNode - 1);
                time = max(arrivalTime, timeWindows(nextNode - 1, 1)); % ����ʱ��
                
                % �ۼ�����ɱ�
                totalCost = totalCost + C * distMatrix(currentNode, nextNode);
                totalCost = totalCost + penalty_coeff(nextNode - 1) * timeViolation;
                totalCost = totalCost + demand(nextNode - 1) * mu * travelTime_ij;
                totalCost = totalCost + G * travelTime_ij;
            else
                route = [route, 1]; % ������������
                load = 0;
                time = 0;
            end
        end
        
        routes{k} = route;
        totalCosts(k) = totalCost;
    end
    
    % ����ȫ�����Ž�
    [minTotalCost, minIndex] = min(totalCosts);
    if minTotalCost < bestTotalCost
        bestTotalCost = minTotalCost;
        bestRoute = routes{minIndex};
    end
    
    % �����������Ż�·��
    bestRoute = twoOpt(bestRoute, distMatrix);
    
    % ��¼��ʷ���Ž�
    costHistory(iter) = bestTotalCost;
    fprintf('���� %d: �����ܳɱ� = %.2f\n', iter, bestTotalCost);
    
    % ��Ϣ�ظ���
    pheromone = (1 - rho) * pheromone;
    for k = 1:numAnts
        route = routes{k};
        for j = 1:length(route) - 1
            pheromone(route(j), route(j + 1)) = pheromone(route(j), route(j + 1)) + Q / totalCosts(k);
        end
    end
end

%% �����ͼ
figure(1);
plot(1:maxIter, costHistory, '-', 'LineWidth', 1.5);
title('��������');
xlabel('��������');
ylabel('�ܳɱ�');
grid on;

figure(2);
plot(nodes(:, 1), nodes(:, 2), 'bo');
hold on;
plot(nodes(1, 1), nodes(1, 2), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k'); % ��������
for i = 2:numCustomers + 1
    text(nodes(i, 1), nodes(i, 2), [' ', num2str(i - 1)]);
end
for i = 1:length(bestRoute) - 1
    plot([nodes(bestRoute(i), 1), nodes(bestRoute(i + 1), 1)], ...
         [nodes(bestRoute(i), 2), nodes(bestRoute(i + 1), 2)], 'r-');
end
title('����·��');
xlabel('X����');
ylabel('Y����');
grid on;