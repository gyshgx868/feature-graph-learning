function write_log(log_file, message, write_to_file)

if ~exist('write_to_file', 'var')
  write_to_file = true;
end
if write_to_file
  handle = fopen(log_file, 'at');
  fprintf(handle, '[%s] ', datestr(now, 'yyyy/mm/dd HH:MM:SS'));
end
fprintf('[%s] ', datestr(now, 'yyyy/mm/dd HH:MM:SS'));
if write_to_file
  fprintf(handle, '%s\n', message);
end
fprintf('%s\n', message);
if write_to_file
  fclose(handle);
end

end
