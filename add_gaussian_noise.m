function add_gaussian_noise(cloud, shapename, variance)
  pts = cloud.Location;
  SAMPLING_SET = 200;
  % estimate diameter
  sampled = down_sample(cloud, SAMPLING_SET);
  Dfps = pdist2(pts, sampled.Location);
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
