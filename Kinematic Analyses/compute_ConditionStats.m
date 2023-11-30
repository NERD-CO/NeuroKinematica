function conditionStats = compute_ConditionStats(data)
conditionStats.mean = mean(data, 'omitnan');
conditionStats.std = std(data, 'omitnan');
conditionStats.var = var(data, 'omitnan');
end
