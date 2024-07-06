
function [output] = par_load(fname, var_str)
  
  load_data_struct = load(fname, var_str);
  output = load_data_struct.(var_str);
  
end