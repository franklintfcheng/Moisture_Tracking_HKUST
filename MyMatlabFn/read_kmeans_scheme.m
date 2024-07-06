function [NoOfClust, Distance, MaxIter, Replicates, rmvalgt, min_clust_size, str_kmeans_info, pass_ratio_thres, n_mov] = read_kmeans_scheme(varargin)
    
    switch nargin
        case 1   %/ given 1 inputs
            kmeans_scheme = varargin{1};
            project_name  = 'CMIP6_MJO';         
        case 2   %/ given 2 inputs
            kmeans_scheme = varargin{1};
            project_name  = varargin{2};
    otherwise
        error('Unexpected inputs')
    end
    fprintf('*** Using Kmeans Scheme %d... ***\n', kmeans_scheme)

%%
    if isequal(project_name, 'CMIP6_MJO')
        if kmeans_scheme == 1
            NoOfClust      = 4;             %/ No of clusters
            Distance       = 'sqeuclidean'; %/ Squared Euclidean distance (default)
            MaxIter        = 10000;         %/ IMPORTANT kmeans param for the results to converge
            Replicates     = 100000;        %/ IMPORTANT kmeans param for finding the best solutions
            rmvalgt        = [];            %/ Remove values greater than xx piror to kmeans
            min_clust_size = 1;             %/ Default
            pass_ratio_thres = [];
            n_mov          = 3;             %/ n-point zonal moving mean on the raw data prior to kmeans

        elseif kmeans_scheme == 2
            NoOfClust      = 4;             %/ No of clusters
            Distance       = 'sqeuclidean'; %/ Squared Euclidean distance (default)
            MaxIter        = 10000;         %/ IMPORTANT kmeans param for the results to converge
            Replicates     = 100000;        %/ IMPORTANT kmeans param for finding the best solutions
            rmvalgt        = -5;            %/ Remove values greater than xx piror to kmeans
            min_clust_size = 1;             %/ Default
            pass_ratio_thres = [];
            n_mov          = 3;             %/ n-point zonal moving mean on the raw data prior to kmeans

        elseif kmeans_scheme == 3
            NoOfClust      = 4;             %/ No of clusters
            Distance       = 'sqeuclidean'; %/ Squared Euclidean distance (default)
            MaxIter        = 10000;         %/ IMPORTANT kmeans param for the results to converge
            Replicates     = 100000;        %/ IMPORTANT kmeans param for finding the best solutions
            rmvalgt        = -5;            %/ Remove values greater than xx piror to kmeans
            min_clust_size = 5;             %/ Make sure of minimum of 5 events in each kmeans cluster by omitting the ac hoc events
            pass_ratio_thres = [];
            n_mov          = 3;             %/ n-point zonal moving mean on the raw data prior to kmeans

        elseif kmeans_scheme == 4
            NoOfClust      = 4;             %/ No of clusters
            Distance       = 'correlation'; %/ correlation (after zero mean and unit standard deviation)
            MaxIter        = 10000;         %/ IMPORTANT kmeans param for the results to converge
            Replicates     = 10000000;      %/ IMPORTANT kmeans param for finding the best solutions
            rmvalgt        = -5;            %/ Remove values greater than xx piror to kmeans
            min_clust_size = 5;             %/ Make sure of minimum of 5 events in each kmeans cluster by omitting the ac hoc events
            pass_ratio_thres = [];
            n_mov          = 3;             %/ n-point zonal moving mean on the raw data prior to kmeans

        elseif kmeans_scheme == 5
            NoOfClust      = 4;             %/ No of clusters
            Distance       = 'sqeuclidean'; %/ correlation (after zero mean and unit standard deviation)
            MaxIter        = 10000;         %/ IMPORTANT kmeans param for the results to converge
            Replicates     = 10000000;        %/ IMPORTANT kmeans param for finding the best solutions
            rmvalgt        = [];            %/ Remove values greater than xx piror to kmeans
            min_clust_size = 5;             %/ Make sure of minimum of 5 events in each kmeans cluster by omitting the ac hoc events
            pass_ratio_thres = [];
            n_mov          = 3;             %/ n-point zonal moving mean on the raw data prior to kmeans

        elseif kmeans_scheme == 6
            NoOfClust      = 4;             %/ No of clusters
            Distance       = 'correlation'; %/ correlation (after zero mean and unit standard deviation)
            MaxIter        = 10000;         %/ IMPORTANT kmeans param for the results to converge
            Replicates     = 300;        %/ IMPORTANT kmeans param for finding the best solutions
            rmvalgt        = -5;            %/ Remove values greater than xx piror to kmeans
            min_clust_size = 5;             %/ Make sure of minimum of 5 events in each kmeans cluster by omitting the ac hoc events
            pass_ratio_thres = 80;
            n_mov          = 3;             %/ n-point zonal moving mean on the raw data prior to kmeans

        elseif kmeans_scheme == 7
            NoOfClust      = 4;             %/ No of clusters
            Distance       = 'correlation'; %/ Squared Euclidean distance (default)
            MaxIter        = 10000;         %/ IMPORTANT kmeans param for the results to converge
            Replicates     = 10000000;           %/ IMPORTANT kmeans param for finding the best solutions
            rmvalgt        = -5;            %/ Remove values greater than xx piror to kmeans
            min_clust_size = 5;             %/ Make sure of minimum of 5 events in each kmeans cluster by omitting the ac hoc events
            pass_ratio_thres = [];          %/ Minimum pass ratio (based on Silhouette score threshold of 0.06)
            n_mov          = 3;             %/ n-point zonal moving mean on the raw data prior to kmeans

        elseif kmeans_scheme == 8
            NoOfClust      = 4;             %/ No of clusters
            Distance       = 'correlation'; %/ Squared Euclidean distance (default)
            MaxIter        = 10000;         %/ IMPORTANT kmeans param for the results to converge
            Replicates     = 10000000;       %/ IMPORTANT kmeans param for finding the best solutions
            rmvalgt        = -5;            %/ Remove values greater than xx piror to kmeans
            min_clust_size = 5;             %/ Make sure of minimum of 5 events in each kmeans cluster by omitting the ac hoc events
            pass_ratio_thres = [];          %/ Minimum pass ratio (based on Silhouette score threshold of 0.06)
            n_mov          = 5;             %/ n-point zonal moving mean on the raw data prior to kmeans
        else
            error('Schemes not set for ''kmeans_scheme == %d''!', kmeans_scheme) 
        end
    else
        error('Schemes not set for project_name = %s!', project_name) 
    end

    % if ~isempty(rmvalgt)
    %     str_rm_val_gt = sprintf('_rmvalgt%.1f', rmvalgt);
    % else
    %     str_rm_val_gt = '';
    % end
    
    str_kmeans_info = sprintf('_KMSch%d', kmeans_scheme); %/ Simplified
end
