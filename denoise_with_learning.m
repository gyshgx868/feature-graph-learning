function [best_mse, best_snr, max_iterations] = denoise_with_learning( ...
  cloud, gt_cloud, model_name, ...
  kMaxIteration, kChangeTolerance, ...
  kDownSampleRate, kPatchSize, kPatchNeighborCount, ...
  alpha, ...
  use_matlab_down_sample, ...
  use_normal, ...
  use_cosine, ...
  use_color, ...
  C, ...
  save_iter_cloud, ...
  log_file, ...
  save_dir)
%% hyper parameters
if ~exist('kMaxIteration', 'var')
  kMaxIteration = 30;
end
if ~exist('alpha', 'var')
  alpha = 25.0 * (exp((1:kMaxIteration) / 20) - 1);
  alpha = repmat(alpha, 5, 1);
  alpha = reshape(alpha, 5 * kMaxIteration, 1);
end
if ~exist('kDownSampleRate', 'var')
  kDownSampleRate = 0.3;
end
if ~exist('kPatchSize', 'var')
  kPatchSize = 9;
end
if ~exist('kPatchNeighborCount', 'var')
  kPatchNeighborCount = 10;
end
if ~exist('kChangeTolerance', 'var')
  kChangeTolerance = 5e-6;
end
if ~exist('use_matlab_down_sample', 'var')
  use_matlab_down_sample = false;
end
if ~exist('use_normal', 'var')
  use_normal = true;
end
if ~exist('use_cosine', 'var')
  use_cosine = false;
end
if ~exist('use_color', 'var')
  use_color = ~isempty(cloud.Color);
end
if ~exist('C', 'var')
  C = 5.0;
end
if ~exist('save_iter_cloud', 'var')
  save_iter_cloud = true;
end
if ~exist('log_file', 'var')
  log_file = './log.txt';
end
if ~exist('save_dir', 'var')
  save_dir = './results';
end

%% prepare global variables
kFeatDims = 3 + 3;  % original coordinates + point featuers (default: coordinates)
if use_normal
  kFeatDims = kFeatDims + 3;
end
if use_color
  kFeatDims = kFeatDims + 3;
end
% save parameters
write_log(log_file, ['Denoising ', model_name, ':']);
write_log(log_file, 'Hyper parameter settings:');
write_log(log_file, '{');
write_log(log_file, ['  kMaxIteration: ', int2str(kMaxIteration)]);
write_log(log_file, ['  kDownSampleRate: ', num2str(kDownSampleRate)]);
write_log(log_file, ['  kPatchSize: ', int2str(kPatchSize)]);
write_log(log_file, ['  kPatchNeighborCount: ', int2str(kPatchNeighborCount)]);
write_log(log_file, ['  kChangeTolerance: ', num2str(kChangeTolerance)]);
write_log(log_file, ['  C: ', num2str(C, '%.2f ')]);
write_log(log_file, ['  use_normal: ', bool2str(use_normal)]);
write_log(log_file, ['  use_cosine: ', bool2str(use_cosine)]);
write_log(log_file, ['  use_color: ', bool2str(use_color)]);
write_log(log_file, ['  use_matlab_down_sample: ', bool2str(use_matlab_down_sample)]);
write_log(log_file, '}');

last_mse = meandistance(gt_cloud.Location, cloud.Location);
last_snr = compute_snr(cloud, last_mse);
write_log(log_file, ['Initial MSE: ', num2str(last_mse, 10)]);
write_log(log_file, ['Initial SNR: ', num2str(last_snr, 10)]);
best_mse = Inf;
best_snr = 0;
max_iterations = 0;

%% Begin iteration
for iter = 1 : kMaxIteration
  write_log(log_file, ['Iteration: ', num2str(iter)]);
  write_log(log_file, ['alpha: ', num2str(alpha(iter), '%.5f')]);
  write_log(log_file, 'Down sampling...');
  if use_matlab_down_sample
    down_sampled_cloud = pcdownsample(cloud, 'gridAverage', 1.4);
  else
    down_sampled_cloud = down_sample( ...
      cloud, round(cloud.Count * kDownSampleRate));
  end
  
  write_log(log_file, 'Building patches...');
  patches = build_patches(cloud, down_sampled_cloud.Location, kPatchSize);
  [patch_count, ~] = size(patches);
  
  centers = zeros(patch_count, 3);
  for i = 1 : patch_count
    patch_points = cloud.Location(patches(i, :), :);
    centers(i, :) = mean(patch_points);
  end
  cx = zeros(patch_count * kPatchSize, 1);
  cy = zeros(patch_count * kPatchSize, 1);
  cz = zeros(patch_count * kPatchSize, 1);
  for i = 1 : patch_count
    for j = 1 : kPatchSize
      cx((i - 1) * kPatchSize + j) = centers(i, 1);
      cy((i - 1) * kPatchSize + j) = centers(i, 2);
      cz((i - 1) * kPatchSize + j) = centers(i, 3);
    end
  end
  
  % build selection matrix
  selection = zeros(kPatchSize * patch_count, 2);
  index = 1;
  for i = 1 : patch_count
    patch = patches(i, :);
    for j = 1 : kPatchSize
      selection(index, 1) = (i - 1) * kPatchSize + j;
      selection(index, 2) = patch(j);
      index = index + 1;
    end
  end
  S = sparse( ...
    selection(:, 1), selection(:, 2), ...
    ones(kPatchSize * patch_count, 1), ...
    kPatchSize * patch_count, cloud.Count, kPatchSize * patch_count);
  
  % build patch connections
  patch_graph = zeros(down_sampled_cloud.Count, kPatchNeighborCount);
  for i = 1 : down_sampled_cloud.Count
    [indices, ~] = findNearestNeighbors( ...
      down_sampled_cloud, ...
      down_sampled_cloud.Location(i, :), ...
      kPatchNeighborCount + 1);
    patch_graph(i, :) = indices(2:end, 1);  % discard self-connection
  end

  pc_centers = zeros(cloud.Count, 3);
  for i = 1 : patch_count
    for p = 1 : kPatchSize
      j = patches(i, p);
      pc_centers(j, :) = centers(i, :);
    end
  end
  % prepare normals for patch graph
  write_log(log_file, 'Computing normals...');
  pc_normals = pcnormals(cloud, 100);
  p_normals = compute_principal_normals(pc_normals, patches);
  cosine = zeros(down_sampled_cloud.Count, kPatchNeighborCount);
  % compute the cosine of pricipal normals
  for i = 1 : down_sampled_cloud.Count
    na = p_normals(i, :);
    nb = p_normals(patch_graph(i, :), :);
    cosine(i, :) = abs(nb * na');
  end
  len_graph = patch_count * kPatchNeighborCount;
  cosine = reshape(cosine', len_graph, 1);
  cosine = reshape(repmat(cosine, 1, kPatchSize)', len_graph * kPatchSize, 1);
  % prepare coordinates
  pc_points = cloud.Location;
  % prepare color
  pc_colors = double(cloud.Color) / 255.0;
  % prepare point featuers
  pc_features = cloud.Location;

  write_log(log_file, 'Building graph...');
  
  laplacian_r = zeros(patch_count * kPatchNeighborCount, kPatchSize);
  laplacian_c = zeros(patch_count * kPatchNeighborCount, kPatchSize);
 
  lap_index = 1;
  % build patch graph
  lhs_features = zeros( ...
    patch_count * kPatchSize * kPatchNeighborCount, kFeatDims);
  rhs_features = zeros( ...
    patch_count * kPatchSize * kPatchNeighborCount, kFeatDims);
  last_from = 1;
  
  % inter-patch connections
  for i = 1 : down_sampled_cloud.Count
    p_indices = patch_graph(i, :);
    for j = 1 : length(p_indices)
      index = p_indices(j);
      if use_normal && use_color
        lhs_features(last_from:last_from + kPatchSize - 1, :) = ...
          [pc_points(patches(i, :), :), ...
           pc_features(patches(i, :), :), ...
           pc_normals(patches(i, :), :), ...
           pc_colors(patches(i, :), :)];
      elseif use_normal
        lhs_features(last_from:last_from + kPatchSize - 1, :) = ...
          [pc_points(patches(i, :), :), ...
           pc_features(patches(i, :), :), ...
           pc_normals(patches(i, :), :)];
      elseif use_color
        lhs_features(last_from:last_from + kPatchSize - 1, :) = ...
          [pc_points(patches(i, :), :), ...
           pc_features(patches(i, :), :), ...
           pc_colors(patches(i, :), :)];
      else
        lhs_features(last_from:last_from + kPatchSize - 1, :) = ...
          [pc_points(patches(i, :), :), ...
           pc_features(patches(i, :), :)];
      end

      r = zeros(1, kPatchSize);
      c = zeros(1, kPatchSize);
      c_points = zeros(kPatchSize, kFeatDims);
      patch_j_points = pc_points(patches(index, :), :);
      patch_j_features = pc_features(patches(index, :), :);
      if use_normal
        patch_j_normals = pc_normals(patches(index, :), :);
      end
      if use_color
        patch_j_colors = pc_colors(patches(index, :), :);
      end
      for p = 1 : kPatchSize
        dist = sum( ...
          abs( ...
            cloud.Location(patches(i, p), :) - centers(i, :) - ...
            (cloud.Location(patches(index, :), :) - centers(index, :)) ...
          ).^2, 2 ...
        ).^0.5;
        
        [~, min_index] = min(dist);
        r(1, p) = (i - 1) * kPatchSize + p;
        c(1, p) = (index - 1) * kPatchSize + min_index;
        if use_normal && use_color
          c_points(p, :) = [ ...
            patch_j_points(min_index, :), ...
            patch_j_features(min_index, :), ...
            patch_j_normals(min_index, :), ...
            patch_j_colors(min_index, :)];
        elseif use_normal
          c_points(p, :) = [ ...
            patch_j_points(min_index, :), ...
            patch_j_features(min_index, :), ...
            patch_j_normals(min_index, :)];
        elseif use_color
          c_points(p, :) = [ ...
            patch_j_points(min_index, :), ...
            patch_j_features(min_index, :), ...
            patch_j_colors(min_index, :)];
        else
          c_points(p, :) = [ ...
            patch_j_points(min_index, :), ...
            patch_j_features(min_index, :)];
        end
      end
      laplacian_r(lap_index, :) = r;
      laplacian_c(lap_index, :) = c;
      lap_index = lap_index + 1;

      rhs_features(last_from:last_from + kPatchSize - 1, :) = c_points;
      last_from = last_from + kPatchSize;
    end
  end

  write_log(log_file, 'Calculating weights...');
  [M, weights, iter_vals, min_iter] = block_gradient_descent( ...
    lhs_features, rhs_features, use_normal, use_cosine, ...
    use_color, C, log_file);

  laplacian_v = weights(:) .* cosine;
  laplacian_v(laplacian_v < 0.001) = 0;
  laplacian_r = laplacian_r';
  laplacian_c = laplacian_c';

  % build Laplacian matrix
  W = sparse( ...
    laplacian_r(:), laplacian_c(:), laplacian_v(:), ...
    kPatchSize * patch_count, kPatchSize * patch_count);
  mask1 = logical(tril(W));
  W = W - mask1' .* W;
  W = W + W';  % make symmetric
  D = diag(sum(W, 2) .^ -0.5);
  L = speye(kPatchSize * patch_count) - D * W * D;

  % solve optimization problem
  write_log(log_file, 'Solving...');
  X = double(cloud.Location(:, 1));
  Y = double(cloud.Location(:, 2));
  Z = double(cloud.Location(:, 3));
  I = speye(cloud.Count);
  rx = lsqr( ...
    alpha(iter) * I + S' * L * S, ...
    alpha(iter) * X + S' * L * cx, ...
    1e-6, ...
    50000);
  ry = lsqr( ...
    alpha(iter) * I + S' * L * S, ...
    alpha(iter) * Y + S' * L * cy, ...
    1e-6, ...
    50000);
  rz = lsqr( ...
    alpha(iter) * I + S' * L * S, ...
    alpha(iter) * Z + S' * L * cz, ...
    1e-6, ...
    50000);
  
  % save results
  to_save = pointCloud([rx ry rz]);
  to_save.Color = cloud.Color;
  mse = meandistance(to_save.Location, gt_cloud.Location);
  snr = compute_snr(to_save, mse);
  write_log(log_file, ['MSE: ', num2str(mse, 5), ', change = ', num2str(last_mse - mse)]);
  write_log(log_file, ['SNR: ', num2str(snr, 5), ', change = ', num2str(snr - last_snr)]);
  best_mse = min(best_mse, mse);
  best_snr = max(best_snr, snr);
  if last_mse - mse < kChangeTolerance
    write_log(log_file, 'Early stopped.');
    break;
  end
  if save_iter_cloud
    save_name = [model_name, '_iter_', num2str(iter), '.ply'];
    save_name = fullfile(save_dir, save_name);
    pcwrite(to_save, save_name);
    write_log(log_file, [save_name, ' saved.']);
  end
  max_iterations = iter;
  
  cloud = to_save;
  last_mse = mse;
  last_snr = snr;
end

if ~save_iter_cloud
  save_name = [model_name, '_denoised.ply'];
  save_name = fullfile(save_dir, save_name);
  pcwrite(cloud, save_name);
  write_log(log_file, [save_name, ' saved.']);
end
write_log(log_file, ['Best MSE: ', num2str(best_mse, 5)]);
write_log(log_file, ['Best SNR: ', num2str(best_snr, 5)]);

end
