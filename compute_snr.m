function [snr] = compute_snr(cloud, mse)
  s = mean(sum(cloud.Location .^ 2, 2) .^ 0.5);
  snr = 10.0 * log(s / mse);
end
