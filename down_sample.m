function [new_cloud] = down_sample(origin, size)
  current_idx = 1; % unidrnd(origin.Count);
  result = zeros(1,min(origin.Count,size));
  d = inf(1, origin.Count);
  for i = 1 : size
    d2 = pdist2(origin.Location(current_idx, :), origin.Location)';
    d = min(d, d2');
    result(i) = current_idx;
    [~, current_idx] = max(d, [], 2);
  end
  new_cloud = pointCloud(origin.Location(result, :));
end
