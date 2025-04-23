clear all;
close all;
clc;

% 設定基本參數
rng(42);  % 設定隨機數種子，確保結果可重現

% GA 算法參數
params.popSize = 200;          % 種群大小
params.crossoverRate = 0.8;    % 交配率
params.mutationRate = 0.1;     % 突變率
params.eliteCount = 10;        % 精英數量
params.generations = 300;      % 迭代次數
params.timeLimit = 300;        % 時間限制（秒）

% 問題參數
params.truckCount = 5;        % 卡車數量
params.maxWorkingHours = 12;   % 每輛卡車最大工作時間（小時）
params.penaltyLate = 1440;     % 延遲懲罰（分鐘/次）
params.penaltyInterrupt = 100; % 中斷懲罰（分鐘/次）
params.setupTime = 3;         % 每次派車準備時間（分鐘）

% 工地資訊
numSites = 3;
siteInfo = struct();
siteInfo.demand = [3,4,5];      % 各工地需求車次
siteInfo.priority = [2, 2, 2];    % 優先級 (1高-3低)  設置為相同優先級
siteInfo.timeWindows = [480, 1440;
    480, 1440;
    510, 1440;
    ];

siteInfo.serviceTime = [20, 30, 25];  % 各工地卸料時間（分鐘）
siteInfo.maxInterrupt = [30,20,15];    % 最大允許中斷時間（分鐘）

% 距離矩陣（廠到各工地及工地間的距離，單位：分鐘）
% 第一行/列為工廠，其餘為工地
distance = [
    0, 30, 25, 40, 15, 35, 45, 20, 25;   % 工廠到各工地
    30, 0, 15, 25, 20, 30, 35, 15, 20;   % 工地1到其他點
    25, 15, 0, 20, 15, 25, 30, 10, 15;   % 工地2到其他點
    40, 25, 20, 0, 30, 20, 15, 25, 30;   % 工地3到其他點
    15, 20, 15, 30, 0, 25, 35, 15, 20;   % 工地4到其他點
    35, 30, 25, 20, 25, 0, 20, 25, 30;   % 工地5到其他點
    45, 35, 30, 15, 35, 20, 0, 30, 35;   % 工地6到其他點
    20, 15, 10, 25, 15, 25, 30, 0, 10;   % 工地7到其他點
    25, 20, 15, 30, 20, 30, 35, 10, 0    % 工地8到其他點
];

% 定義要進行敏感度分析的參數
% 格式：[參數名稱，基準值，變化範圍百分比，變化點數]
sensitivityParams = {
    {'卡車數量', params.truckCount, 0.5, 9},              % 卡車數量變化範圍 ±50%
    {'交配率', params.crossoverRate, 0.5, 11},            % 交配率變化範圍 ±50%
    {'突變率', params.mutationRate, 0.8, 11},             % 突變率變化範圍 ±80%
    {'延遲懲罰', params.penaltyLate, 0.5, 9},             % 延遲懲罰變化範圍 ±50%
    {'中斷懲罰', params.penaltyInterrupt, 0.6, 9},        % 中斷懲罰變化範圍 ±60%
    {'種群大小', params.popSize, 0.5, 9}                  % 種群大小變化範圍 ±50%
};

numParams = length(sensitivityParams);
results = cell(numParams, 1);
sensitivityCoeff = zeros(numParams, 1);

% 執行基準案例評估
fprintf('執行基準案例評估...\n');
baseResult = evaluateDispatchPlan(params, siteInfo, distance, numSites);
fprintf('基準案例目標函數值: %.2f\n\n', baseResult);

% 進行敏感度分析
figure('Position', [100, 100, 1200, 800]);

for i = 1:numParams
    paramName = sensitivityParams{i}{1};
    baseValue = sensitivityParams{i}{2};
    variationPct = sensitivityParams{i}{3};
    numPoints = sensitivityParams{i}{4};
    
    fprintf('分析參數 %s 的敏感度...\n', paramName);
    
    % 創建變化範圍
    if strcmp(paramName, '交配率') || strcmp(paramName, '突變率')
        % 對於0-1範圍的參數，使用絕對值變化
        variationRange = linspace(-baseValue * variationPct, (1-baseValue) * variationPct, numPoints);
    else
        % 對於其他參數，使用百分比變化
        variationRange = linspace(-variationPct, variationPct, numPoints);
    end
    
    paramValues = zeros(numPoints, 1);
    resultValues = zeros(numPoints, 1);
    
    for j = 1:numPoints
        % 複製參數
        testParams = params;
        
        % 設定當前參數值
        if strcmp(paramName, '卡車數量')
            adjustedValue = round(baseValue * (1 + variationRange(j)));
            testParams.truckCount = max(1, adjustedValue);
            paramValues(j) = testParams.truckCount;
        elseif strcmp(paramName, '交配率')
            adjustedValue = baseValue + variationRange(j);
            testParams.crossoverRate = max(0.01, min(0.99, adjustedValue));
            paramValues(j) = testParams.crossoverRate;
        elseif strcmp(paramName, '突變率')
            adjustedValue = baseValue + variationRange(j);
            testParams.mutationRate = max(0.01, min(0.99, adjustedValue));
            paramValues(j) = testParams.mutationRate;
        elseif strcmp(paramName, '延遲懲罰')
            adjustedValue = baseValue * (1 + variationRange(j));
            testParams.penaltyLate = max(1, adjustedValue);
            paramValues(j) = testParams.penaltyLate;
        elseif strcmp(paramName, '中斷懲罰')
            adjustedValue = baseValue * (1 + variationRange(j));
            testParams.penaltyInterrupt = max(1, adjustedValue);
            paramValues(j) = testParams.penaltyInterrupt;
        elseif strcmp(paramName, '種群大小')
            adjustedValue = round(baseValue * (1 + variationRange(j)));
            testParams.popSize = max(10, adjustedValue);
            paramValues(j) = testParams.popSize;
        end
        
        % 評估當前參數組合
        resultValues(j) = evaluateDispatchPlan(testParams, siteInfo, distance, numSites);
        fprintf('  %s = %.2f, 目標函數值 = %.2f\n', paramName, paramValues(j), resultValues(j));
    end
    
    % 存儲結果
    results{i} = [paramValues, resultValues];
    
    % 計算敏感度係數（使用線性迴歸近似）
    if ~(strcmp(paramName, '交配率') || strcmp(paramName, '突變率'))
        % 對於百分比變化的參數
        normalizedParams = (paramValues - baseValue) / baseValue;
    else
        % 對於0-1範圍的參數
        normalizedParams = (paramValues - baseValue);
    end
    normalizedResults = (resultValues - baseResult) / baseResult;
    
    % 線性迴歸計算斜率（敏感度係數）
    p = polyfit(normalizedParams, normalizedResults, 1);
    sensitivityCoeff(i) = p(1);
    
    % 繪製敏感度曲線
    subplot(2, 3, i);
    plot(paramValues, resultValues, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
    title(['參數: ', paramName, ' 的敏感度分析']);
    xlabel(paramName);
    ylabel('目標函數值');
    grid on;
    
    % 繪製擬合線
    hold on;
    xRange = linspace(min(paramValues), max(paramValues), 100);
    if ~(strcmp(paramName, '交配率') || strcmp(paramName, '突變率'))
        yFit = baseResult * (1 + p(1) * ((xRange - baseValue) / baseValue));
    else
        yFit = baseResult * (1 + p(1) * (xRange - baseValue));
    end
    plot(xRange, yFit, '--r', 'LineWidth', 1.5);
    hold off;
end

% 顯示敏感度係數
figure;
barNames = cell(numParams, 1);
for i = 1:numParams
    barNames{i} = sensitivityParams{i}{1};
end

% 按敏感度絕對值排序
[~, sortIdx] = sort(abs(sensitivityCoeff), 'descend');
sortedCoeff = sensitivityCoeff(sortIdx);
sortedNames = barNames(sortIdx);


%絕對數值越大對適應度值影響越大
barh(sortedCoeff);
set(gca, 'YTick', 1:numParams, 'YTickLabel', sortedNames);
title('參數敏感度係數');
xlabel('敏感度係數（正值表示正相關，負值表示負相關）');
grid on;

% 輸出敏感度係數
fprintf('\n敏感度係數（按絕對值排序）:\n');
for i = 1:numParams
    fprintf('參數 %s: %.4f\n', sortedNames{i}, sortedCoeff(i));
end

% 視覺化基準案例的調度計劃
visualizeDispatchPlan(params, siteInfo, distance, numSites);

% ===================== 評估函式 =====================
function objective = evaluateDispatchPlan(params, siteInfo, distance, numSites)
    % 在實際應用中，此函數應該調用完整的GA算法
    % 這裡使用一個更合理的簡化模型來模擬GA的效果
    
    % 提取參數
    truckCount = params.truckCount;
    crossoverRate = params.crossoverRate;
    mutationRate = params.mutationRate;
    penaltyLate = params.penaltyLate;
    penaltyInterrupt = params.penaltyInterrupt;
    popSize = params.popSize;
    generations = params.generations;
    maxWorkingHours = params.maxWorkingHours;
    setupTime = params.setupTime;
    
    % 提取工地資訊
    demand = siteInfo.demand;
    priority = siteInfo.priority;
    timeWindows = siteInfo.timeWindows;
    serviceTime = siteInfo.serviceTime;
    maxInterrupt = siteInfo.maxInterrupt;
    
    % 計算總需求車次
    totalDemand = sum(demand);
    
    % 計算理論最小需要卡車數量
    totalServiceTime = sum(demand .* serviceTime);
    avgTravelTime = mean(distance(1, 2:end)) + mean(mean(distance(2:end, 2:end)));
    avgTripTime = avgTravelTime + mean(serviceTime) + setupTime;
    theoreticalMinTrucks = ceil(totalDemand * avgTripTime / (maxWorkingHours * 60));
    
    % 計算時間窗限制因子
    timeWindowFactor = 0;
    for i = 1:numSites
        availableTime = timeWindows(i, 2) - timeWindows(i, 1);
        requiredTime = demand(i) * serviceTime(i);
        travelTime = distance(1, i+1) * 2 * demand(i);  % 假設每次從工廠出發
        
        if requiredTime + travelTime > availableTime
            timeWindowFactor = timeWindowFactor + (requiredTime + travelTime - availableTime);
        end
    end
    
    % 考慮優先級
   priorityFactor = 1;
    
    % 計算中斷時間影響
    interruptFactor = sum(maxInterrupt .* demand) / sum(demand);
    
    % 考慮GA參數的影響
    % 基於實際經驗，GA參數對結果的影響通常呈非線性關係
    gaEfficiency = 0.5 + 0.3 * (crossoverRate ^ 0.7) + 0.2 * (mutationRate ^ 0.3) * (1 - mutationRate ^ 2);
    gaScaleFactor = log(popSize) / log(200) * log(generations) / log(300);
    
    % 卡車利用率影響
    if truckCount < theoreticalMinTrucks
        % 卡車不足，大幅增加成本
        utilizationPenalty = (1 + (theoreticalMinTrucks - truckCount) / theoreticalMinTrucks * 5);
    else
        % 卡車富餘，略微增加成本
        utilizationPenalty = 1 + 0.02 * (truckCount - theoreticalMinTrucks) / theoreticalMinTrucks;
    end
    
    % 計算目標函數值（較低值表示更好的解決方案）
    % 合理化各個因素的權重
    baseValue = totalDemand * avgTripTime * utilizationPenalty;
    timeWindowCost = timeWindowFactor * 0.05 * penaltyLate;
    interruptCost = (maxInterrupt(1) / interruptFactor) * penaltyInterrupt * 0.02 * totalDemand;
    priorityCost = 0;
    
    % 最終目標函數值
    objective = (baseValue + timeWindowCost + interruptCost + priorityCost) / (gaEfficiency * gaScaleFactor);
    
    % 添加小幅度隨機性以模擬GA的隨機特性，但控制在較小範圍內
    randomVariation = 0.98 + 0.04 * rand();
    objective = objective * randomVariation;
end

% ===================== 視覺化函式 =====================
function visualizeDispatchPlan(params, siteInfo, distance, numSites)
    % 創建一個簡單的示例調度計劃來視覺化
    % 在實際應用中，這應該使用GA算法的實際最佳解
    
    % 提取參數
    truckCount = params.truckCount;
    demand = siteInfo.demand;
    timeWindows = siteInfo.timeWindows;
    serviceTime = siteInfo.serviceTime;
    
    % 創建視覺化界面
    figure('Position', [100, 100, 1200, 800], 'Name', '預拌混凝土車派遣示例');
    
    % 1. 創建工地需求圖
    subplot(2, 2, 1);
    bar(1:numSites, demand);
    title('各工地需求車次');
    xlabel('工地編號');
    ylabel('需求車次');
    grid on;
    
    % 2. 創建時間窗口圖
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
    xlim([7, 18]);  % 顯示7:00到18:00
    grid on;
    hold off;
    
    % 3. 簡單的甘特圖來展示卡車調度
    subplot(2, 2, [3, 4]);
    hold on;
    
    % 生成一個簡單的調度計劃用於展示
    % 這個部分在實際應用中應該使用GA算法的結果
    
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
            % 1. 從廠到工地的旅程
            t1 = truckStartTimes(truckIdx) - returnTime + travelToSite;
            plot([t1, arrivalTime], [truckIdx, truckIdx], '-', 'Color', colors(truckIdx,:), 'LineWidth', 2);
            
            % 2. 在工地的服務時間
            rectangle('Position', [arrivalTime, truckIdx-0.3, serviceTime(site), 0.6], ...
                'FaceColor', colors(truckIdx,:), 'EdgeColor', 'k');
            
            % 3. 從工地返回廠的旅程
            plot([completionTime, returnTime], [truckIdx, truckIdx], '--', 'Color', colors(truckIdx,:), 'LineWidth', 2);
            
            % 添加工地標籤
            text(arrivalTime + serviceTime(site)/2, truckIdx, ['S', num2str(site)], ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontWeight', 'bold', 'Color', 'white');
        end
    end
    
    % 設置圖表屬性
    title('預拌混凝土車派遣甘特圖 (示例)');
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
    
    % 添加圖例
    for t = 1:min(5, truckCount)  % 只顯示前5輛卡車的圖例以避免擁擠
        plot([0, 0], [0, 0], '-', 'Color', colors(t,:), 'LineWidth', 2, 'DisplayName', ['車輛 ', num2str(t)]);
    end
    legend('Location', 'NorthEast');
    
    hold off;
    
    % 添加參數信息
    annotation('textbox', [0.02, 0.02, 0.25, 0.1], 'String', ...
        sprintf('卡車數量: %d\n交配率: %.2f\n突變率: %.2f\n延遲懲罰: %.0f\n中斷懲罰: %.0f', ...
        params.truckCount, params.crossoverRate, params.mutationRate, ...
        params.penaltyLate, params.penaltyInterrupt), ...
        'EdgeColor', 'none', 'BackgroundColor', [0.95, 0.95, 0.95], 'FontSize', 9);
end