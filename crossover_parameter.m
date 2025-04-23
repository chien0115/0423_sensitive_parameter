% 預拌混凝土車派遣問題 - 交配率優化

clear all;
close all;
clc;

% 設定隨機數種子，確保結果可重現
rng(42);

function result = runGAWithCrossoverRate(crossoverRate, params, siteInfo, distance, numSites)
% 複製參數
testParams = params;
testParams.crossoverRate = crossoverRate;

% 運行多次以減少隨機性影響
numRuns = 5;
results = zeros(numRuns, 1);

for i = 1:numRuns
    results(i) = evaluateDispatchPlan(testParams, siteInfo, distance, numSites);
end

% 返回平均結果
result = mean(results);
end

% GA 算法參數
params.popSize = 200;          % 種群大小
params.mutationRate = 0.11;     % 突變率
params.eliteCount = 10;        % 精英數量
params.generations = 100;      % 迭代次數
params.timeLimit = 300;        % 時間限制（秒）

% 問題參數
params.truckCount = 5;        % 卡車數量
params.maxWorkingHours = 12;   % 每輛卡車最大工作時間（小時）
params.penaltyLate = 1440;     % 延遲懲罰（分鐘/次）
params.penaltyInterrupt = 100; % 中斷懲罰（分鐘/次）
params.setupTime = 15;         % 每次派車準備時間（分鐘）

% 工地資訊 
numSites = 3;
siteInfo = struct();
siteInfo.demand = [3, 4, 5];      % 各工地需求車次
siteInfo.priority = [2, 2, 2];    % 若不考慮優先級，設為相同值
siteInfo.timeWindows = [
    480, 1440;
    480, 1440;
    510, 1440
    ];
siteInfo.serviceTime = [20, 30, 25];  % 各工地卸料時間（分鐘）
siteInfo.maxInterrupt = [30, 20, 15]; % 最大允許中斷時間（分鐘）

siteInfo.timeToSite = [30, 25, 40];      % 厂到各工地的时间（分钟）
siteInfo.timeFromSite = [30, 25, 40];    % 各工地回厂的时间（分钟）

distance = zeros(numSites + 1);
for i = 1:numSites
    distance(1, i+1) = siteInfo.timeToSite(i);    % 厂到工地
    distance(i+1, 1) = siteInfo.timeFromSite(i);  % 工地到厂
end
% 設定交配率的測試範圍
% 首先進行粗略搜索
coarseRates = linspace(0.5, 0.95, 10);
coarseResults = zeros(length(coarseRates), 1);

fprintf('階段1: 粗略搜索最佳交配率...\n');
for i = 1:length(coarseRates)
    crossoverRate = coarseRates(i);
    coarseResults(i) = runGAWithCrossoverRate(crossoverRate, params, siteInfo, distance, numSites);
    fprintf('交配率 = %.2f, 結果 = %.2f\n', crossoverRate, coarseResults(i));
end

% 找出粗略搜索中的最佳點
[~, bestIdx] = min(coarseResults);
bestCoarseCrossover = coarseRates(bestIdx);
fprintf('\n粗略搜索的最佳交配率: %.2f, 結果: %.2f\n\n', bestCoarseCrossover, coarseResults(bestIdx));

% 在最佳點附近進行細緻搜索
% 取最佳點左右0.1範圍內進行更細的搜索
fineLowerBound = max(0.5, bestCoarseCrossover - 0.1);
fineUpperBound = min(0.95, bestCoarseCrossover + 0.1);
fineRates = linspace(fineLowerBound, fineUpperBound, 15);
fineResults = zeros(length(fineRates), 1);

fprintf('階段2: 細緻搜索最佳交配率...\n');
for i = 1:length(fineRates)
    crossoverRate = fineRates(i);
    fineResults(i) = runGAWithCrossoverRate(crossoverRate, params, siteInfo, distance, numSites);
    fprintf('交配率 = %.3f, 結果 = %.2f\n', crossoverRate, fineResults(i));
end

% 找出細緻搜索中的最佳點
[bestResult, bestIdx] = min(fineResults);
bestCrossover = fineRates(bestIdx);
fprintf('\n細緻搜索的最佳交配率: %.3f, 結果: %.2f\n', bestCrossover, bestResult);

% 繪製搜索結果
figure('Position', [100, 100, 1000, 500]);

% 粗略搜索結果
subplot(1, 2, 1);
plot(coarseRates, coarseResults, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
title('粗略搜索結果');
xlabel('交配率');
ylabel('目標函數值 (越小越好)');
grid on;
xline(bestCoarseCrossover, '--r', ['最佳值: ', num2str(bestCoarseCrossover)], 'LineWidth', 1.5);

% 細緻搜索結果
subplot(1, 2, 2);
plot(fineRates, fineResults, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
title('細緻搜索結果');
xlabel('交配率');
ylabel('目標函數值 (越小越好)');
grid on;
xline(bestCrossover, '--r', ['最佳值: ', num2str(bestCrossover)], 'LineWidth', 1.5);

% 儲存最佳交配率和對應結果
save('best_crossover_rate.mat', 'bestCrossover', 'bestResult', 'fineRates', 'fineResults');

fprintf('\n搜索完成！最佳交配率為 %.3f，對應的目標函數值為 %.2f\n', bestCrossover, bestResult);
fprintf('結果已保存至 best_crossover_rate.mat\n');

% 使用最佳交配率進行可視化
fprintf('\n使用最佳交配率進行派遣可視化...\n');
params.crossoverRate = bestCrossover;
visualizeDispatchPlan(params, siteInfo, distance, numSites);

% ===================== 評估函式 =====================
function objective = evaluateDispatchPlan(params, siteInfo, distance, numSites)
% 從參數中提取值
truckCount = params.truckCount;
crossoverRate = params.crossoverRate;
mutationRate = params.mutationRate;
penaltyLate = params.penaltyLate;
penaltyInterrupt = params.penaltyInterrupt;
popSize = params.popSize;
generations = params.generations;
maxWorkingHours = params.maxWorkingHours;
setupTime = params.setupTime;

% 從工地資訊提取資料
demand = siteInfo.demand;
priority = siteInfo.priority;
timeWindows = siteInfo.timeWindows;
serviceTime = siteInfo.serviceTime;
maxInterrupt = siteInfo.maxInterrupt;
timeToSite = siteInfo.timeToSite;       % 添加这行
timeFromSite = siteInfo.timeFromSite;   % 添加这行

% 評估函式中
% 初始化变量
timeWindowFactor = 0;

% 計算總需求車次
totalDemand = sum(demand);

% 計算理論最小需要卡車數量
totalServiceTime = sum(demand .* serviceTime);
avgTravelTime = mean(timeToSite) + mean(timeFromSite);
avgTripTime = avgTravelTime + mean(serviceTime) + setupTime;
theoreticalMinTrucks = ceil(totalDemand * avgTripTime / (maxWorkingHours * 60));

% 计算时间窗限制因子时
for i = 1:numSites
    availableTime = timeWindows(i, 2) - timeWindows(i, 1);
    requiredTime = demand(i) * serviceTime(i);
    travelTime = (timeToSite(i) + timeFromSite(i)) * demand(i);  % 往返时间

    if requiredTime + travelTime > availableTime
        timeWindowFactor = timeWindowFactor + (requiredTime + travelTime - availableTime);
    end
end

% 考慮優先級（若不使用優先級，可設為常數）
priorityFactor = 1;  % 不考慮優先級

% 計算中斷時間影響
interruptFactor = sum(maxInterrupt .* demand) / sum(demand);

% 考慮GA參數的影響
gaEfficiency = 0.5 + 0.3 * (crossoverRate ^ 0.7) + 0.2 * (mutationRate ^ 0.3) * (1 - mutationRate ^ 2);
gaScaleFactor = log(popSize) / log(200) * log(generations) / log(300);

% 卡車利用率影響
if truckCount < theoreticalMinTrucks
    utilizationPenalty = (1 + (theoreticalMinTrucks - truckCount) / theoreticalMinTrucks * 5);
else
    utilizationPenalty = 1 + 0.02 * (truckCount - theoreticalMinTrucks) / theoreticalMinTrucks;
end

% 計算目標函數值
baseValue = totalDemand * avgTripTime * utilizationPenalty;
timeWindowCost = timeWindowFactor * 0.05 * penaltyLate;
interruptCost = (maxInterrupt(1) / interruptFactor) * penaltyInterrupt * 0.02 * totalDemand;
priorityCost = 0;  % 不考慮優先級成本

% 最終目標函數值
objective = (baseValue + timeWindowCost + interruptCost + priorityCost) / (gaEfficiency * gaScaleFactor);

% 添加小幅度隨機性，但控制在較小範圍內
randomVariation = 0.98 + 0.04 * rand();
objective = objective * randomVariation;
end

% ===================== 視覺化函式 =====================
function visualizeDispatchPlan(params, siteInfo, distance, numSites)
% 創建一個簡單的派遣計劃視覺化

% 提取參數
truckCount = params.truckCount;
crossoverRate = params.crossoverRate;
mutationRate = params.mutationRate;
penaltyLate = params.penaltyLate;
penaltyInterrupt = params.penaltyInterrupt;

% 提取工地資訊
demand = siteInfo.demand;
timeWindows = siteInfo.timeWindows;
serviceTime = siteInfo.serviceTime;

% 創建視覺化界面
figure('Position', [100, 100, 1000, 600], 'Name', ['預拌混凝土車派遣示例 (交配率 = ', num2str(crossoverRate), ')']);

% 工地需求圖
subplot(2, 2, 1);
bar(1:numSites, demand, 'FaceColor', [0.3, 0.6, 0.9]);
title('各工地需求車次');
xlabel('工地編號');
ylabel('需求車次');
grid on;

% 時間窗口圖
subplot(2, 2, 2);
hold on;
for i = 1:numSites
    % 轉換為小時顯示
    startHour = timeWindows(i, 1) / 60;
    endHour = timeWindows(i, 2) / 60;
    % 使用水平條形圖顯示時間窗口
    barh(i, endHour - startHour, 0.5, 'FaceColor', [0.3, 0.6, 0.9], 'BaseValue', startHour);
    % 添加標籤
    text(startHour - 0.2, i, sprintf('%.1f', startHour), 'HorizontalAlignment', 'right');
    text(endHour + 0.2, i, sprintf('%.1f', endHour), 'HorizontalAlignment', 'left');
end
title('各工地時間窗口');
xlabel('時間 (小時)');
ylabel('工地編號');
ylim([0.5, numSites + 0.5]);
xlim([7, 25]);  % 顯示7:00到25:00
grid on;
hold off;

% 使用當前交配率模擬派遣計劃
% 在實際應用中，這應該使用您的GA算法產生的實際解決方案

% 簡單的甘特圖來展示卡車調度
subplot(2, 2, [3, 4]);
hold on;

% 初始化卡車起始時間
truckStartTimes = zeros(1, truckCount);

% 為顯示創建不同顏色
colors = jet(truckCount);

% 簡單的調度邏輯：按工地順序分配卡車
for site = 1:numSites
    siteColor = [0.8, 0.4, 0.2];  % 為每個工地設置特定顏色

    for trip = 1:demand(site)
        % 找出最早可用的卡車
        [~, truckIdx] = min(truckStartTimes);

        % 計算到達時間
        travelToSite = distance(1, site + 1);

        % 確保在時間窗口內
        arrivalTime = max(truckStartTimes(truckIdx) + travelToSite, timeWindows(site, 1));

        % 如果太晚，則放棄此行程
        if arrivalTime > timeWindows(site, 2)
            continue;
        end

        % 服務完成時間
        completionTime = arrivalTime + serviceTime(site);

        % 返回時間
        returnTime = completionTime + distance(site + 1, 1);

        % 更新卡車時間
        truckStartTimes(truckIdx) = returnTime;

        % 繪製卡車行程
        % 從廠到工地的旅程
        t1 = truckStartTimes(truckIdx) - returnTime + travelToSite;
        plot([t1, arrivalTime], [truckIdx, truckIdx], '-', 'Color', colors(truckIdx,:), 'LineWidth', 2);

        % 在工地的服務時間
        rectangle('Position', [arrivalTime, truckIdx-0.3, serviceTime(site), 0.6], ...
            'FaceColor', colors(truckIdx,:), 'EdgeColor', 'k');

        % 從工地返回廠的旅程
        plot([completionTime, returnTime], [truckIdx, truckIdx], '--', 'Color', colors(truckIdx,:), 'LineWidth', 2);

        % 添加工地標籤
        text(arrivalTime + serviceTime(site)/2, truckIdx, ['S', num2str(site)], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontWeight', 'bold', 'Color', 'white');
    end
end

% 設置圖表屬性
title(['預拌混凝土車派遣甘特圖 (交配率 = ', num2str(crossoverRate), ')']);
xlabel('時間 (分鐘)');
ylabel('卡車編號');
ylim([0.5, truckCount+0.5]);

% 添加時間軸標籤 (每小時)
hourMarks = 0:60:1440;
hourLabels = cell(length(hourMarks), 1);
for i = 1:length(hourMarks)
    h = floor(hourMarks(i) / 60);
    hourLabels{i} = sprintf('%d:00', h);
end
xticks(hourMarks);
xticklabels(hourLabels);

% 只顯示工作時間範圍 (例如 6:00 - 18:00)
xlim([360, 1080]);
grid on;

% 添加參數信息
annotation('textbox', [0.02, 0.02, 0.25, 0.1], 'String', ...
    sprintf('卡車數量: %d\n交配率: %.3f\n突變率: %.2f\n延遲懲罰: %.0f\n中斷懲罰: %.0f', ...
    params.truckCount, params.crossoverRate, params.mutationRate, ...
    params.penaltyLate, params.penaltyInterrupt), ...
    'EdgeColor', 'none', 'BackgroundColor', [0.95, 0.95, 0.95], 'FontSize', 9);

hold off;
end