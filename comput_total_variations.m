function [result] = comput_total_variations(cloud_normal, patches)
  [rows, ~] = size(patches);
  result = zeros(rows, 1);
  for index = 1 : rows
    patch = patches(index, :);
    sum_angle = 0.0;
    pair_count = 0;
    for i = 1 : length(patch)
      na = cloud_normal(patch(i), :);
      for j = 1 : length(patch)
        if patch(i) == patch(j)
          continue;
        end
        nb = cloud_normal(patch(j), :);
        sum_angle = sum_angle + (1.0 - abs(na * nb' / (norm(na) * norm(nb))));
        pair_count = pair_count + 1;
      end
    end
    sum_angle = pair_count;
    result(index) = sum_angle / pair_count;
  end
end
