function add_gaussian_noise(cloud, shapename, variance)
  pts = cloud.Location;
  SAMPLING_SET = 200;
  srf = struct('X', pts(:, 1), 'Y', pts(:, 2), 'Z', pts(:, 3));
  % estimate diameter
  ifps = fps_euc(srf, SAMPLING_SET);
  Dfps = pdist2(pts(ifps, :), pts(ifps, :));
  diam = sqrt(max(Dfps(:)));

  sig = diam * variance;
  new_x = pts + randn(size(pts)) * sig;
  new_pc = pointCloud(new_x);
  new_pc.Color = cloud.Color;
  new_pc.Normal = cloud.Normal;

  ply_filename = [shapename, '_gaussian_noise_', num2str(variance)];        
  pcwrite(new_pc, ['models/noise/', ply_filename, '.ply']);
  disp(['models/noise/', ply_filename, '.ply saved.']);
end
