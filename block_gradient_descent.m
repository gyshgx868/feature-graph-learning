function [M, weights, cur_val, min_iter] = block_gradient_descent( ...
  patch_s, patch_t, use_normal, use_cosine, use_color, C, log_file)

if ~exist('use_normal', 'var')
  use_normal = true;
end
if ~exist('use_cosine', 'var')
  use_cosine = false;
end
if ~exist('use_color', 'var')
  use_color = false;
end
if ~exist('C', 'var')
  C = 5.0;
end
if ~exist('log_file', 'var')
  log_file = './log.txt';
end

indices = 1:size(patch_s, 1);

% initialize variables
diff = patch_s(indices, 4:6) - patch_t(indices, 4:6);  % 9
if use_normal && use_cosine
  norm_s = sum(patch_s(indices, 7:9).^2, 2).^0.5;
  norm_t = sum(patch_t(indices, 7:9).^2, 2).^0.5;
  cos_diff = sum(patch_s(indices, 7:9) .* patch_t(indices, 7:9), 2);
  cos_diff = 1.0 - abs(cos_diff ./ norm_s ./ norm_t);
  diff = [diff, cos_diff];
elseif use_normal
  diff = [diff, patch_s(indices, 7:9) - patch_t(indices, 7:9)];
end
if use_color
  diff = [diff, patch_s(indices, 10:12) - patch_t(indices, 10:12)];
end
dist = sum(abs(patch_s(indices, 1:3) - patch_t(indices, 1:3)).^2, 2);

[node_dim, feat_dim] = size(diff);

max_algo_iter = 10;
max_diag_iter = 10;
max_block_iter = 10;
step = 1e-6;

M = eye(feat_dim) * C / feat_dim;
min_val = Inf;
min_M = M;
min_iter = 1;
cur_val = zeros(max_algo_iter, 1);
diag_val = zeros(max_diag_iter, 1);

for iter = 1 : max_algo_iter
  % optimize block
  for f = 1 : feat_dim
    M21 = [M(1:f-1, f); M(f+1:end, f)];
    M22 = [M(1:f-1, 1:f-1), M(1:f-1, f+1:end); ...
      M(f+1:end, 1:f-1), M(f+1:end, f+1:end)];
    g1 = diff(:, f);
    g2 = [diff(:, 1:f-1), diff(:, f+1:end)];
    theta_min = eigs(M22, 1, 'smallestabs');
    dist_hat = exp(-sum(g2*M22.*g2, 2)) .* dist;
    m11 = M(f, f);
    for block_iter = 1 : max_block_iter
      exp_obj = exp(-(g1.^2*m11 + 2.0*sum(g1*M21'.*g2,2)));
      grad = -2.0 * sum(g1.*g2.*exp_obj.*dist_hat);
      M21 = M21 - step * grad';
      if norm(M21, 2) > sqrt(theta_min * m11)
        M21 = M21 * sqrt(theta_min * m11) / norm(M21, 2);
      end
    end
    M(1:f-1, f) = M21(1:f-1);
    M(f+1:end, f) = M21(f:end);
    M(f, 1:f-1) = M21(1:f-1)';
    M(f, f+1:end) = M21(f:end)';
  end
  
  % optimize diagonal
  for diag_iter = 1 : max_diag_iter
    diag_M = diag(M);
    dist_tilde = exp(-sum(diff * (M - diag(diag_M)) .* diff, 2)) .* dist;
    row_M = sum(abs(M), 2) - abs(diag_M);

    grad_M = zeros(feat_dim, 1);
    temp = 0;
    for f = 1 : feat_dim
      temp = temp + diff(:, f).^2 * diag_M(f);
    end
    temp = exp(-temp);
    for f = 1 : feat_dim
      grad_M(f) = -sum(temp .* diff(:, f).^2 .* dist_tilde);
    end
    diag_M = diag_M - step * grad_M;
    temp_diag = min(max(diag_M, row_M), C);
    if sum(temp_diag) <= C
      diag_M = temp_diag;
    else
      diag_M = diag_M - (sum(diag_M) - C) / feat_dim;
      diag_M = min(max(diag_M, row_M), C);
    end
    M(1:feat_dim+1:end) = diag_M;
    
    weights = exp(-sum((diff * M) .* diff, 2));
    diag_val(diag_iter) = sum(weights .* dist);
  end

  weights = exp(-sum((diff * M) .* diff, 2));
  iter_val = sum(weights .* dist);
  cur_val(iter) = iter_val;
  if iter_val < min_val
    min_val = iter_val;
    min_iter = iter;
    min_M = M;
  end

  write_log(log_file, ['  ', int2str(iter), '/', int2str(max_algo_iter), ' done'], false);
end

M = min_M;
weights = exp(-sum((diff * M) .* diff, 2));

end
