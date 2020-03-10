function [result] = build_patches(cloud, centers, size)
  result = zeros(length(centers), size);
  for i = 1 : length(centers)
    [indices, ~] = findNearestNeighbors(cloud, centers(i, :), size);
    result(i, :) = indices';
  end
end
