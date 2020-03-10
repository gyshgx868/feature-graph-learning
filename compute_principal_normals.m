function [result] = compute_principal_normals(cloud_normal, patches)
  [rows, ~] = size(patches);
  result = zeros(rows, 3);
  for i = 1 : rows
    patch_normal = cloud_normal(patches(i, :), :);
    sum_normal = sum(patch_normal);
    normalized = norm(sum_normal);
    if normalized < 1e-8
      normalized = 1.0;
    end
    result(i, :) = sum_normal / normalized;
  end
end
