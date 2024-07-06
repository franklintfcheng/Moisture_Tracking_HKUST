%%
function ip = inParallel()
%/ Description: This function is to check if you are running on a parallel
%/              worker. Running parfor on a worker will throw an error.
%/
%/              See also https://www.mathworks.com/matlabcentral/answers/1982-am-i-inside-a-parfor

    job = getCurrentJob();
    ip = ~isempty(job);
end