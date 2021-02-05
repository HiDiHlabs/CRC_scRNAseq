%% Non-negative matrix factorisation
% as implemented in Regev HNSCC paper
% after selecting for ~ 8,000 genes based on Regev approach 
% from oligodendroglioma paper

% set folder name
data.loc = ### INSERT DATA FOLDER HERE ###

% patient IDs
pid = [90121; 92047; 90130; 92492; 91526; 92507; 92460; 92489; 92074; 100399];

% set number of factors
k = 30;

opt = statset( 'MaxIter', 10, ...
                'TolFun', 1e-10, ...
                'TolX', 1e-10, ...
                'Display', 'iter');

% use mean-centered, log-transformed expression values
data = csvread(paste0(data.loc, 'NNMF_new10/nnmf_input_matrix_min0_new10.csv'),1,1);
rng(1)
[Factors_g,Factors_c,D] = nnmf(data,k, ...
                                    'algorithm', 'als', ...
                                    'replicates', 1);
filename_g = paste0(data.loc, 'NNMF/nnmf_new10_genes_k30_min0_defaultParams.csv');
filename_c = paste0(data.loc, 'NNMF/nnmf_new10_cells_k30_min0_defaultParams.csv');    
csvwrite(filename_g,Factors_g)
csvwrite(filename_c,transpose(Factors_c))   


for p=1:length(pid)
    data = csvread(strcat(paste0(data.loc, 'NNMF/nnmf_input_matrix_min0_Pat-'),num2str(pid(p)),'.csv'),1,1);
    data = data(1:100,1:100);
    rng(1)
    [Factors_g,Factors_c,D] = nnmf(data,k, ...
                                    'algorithm', 'mult', ...
                                    'replicates', 10, ...
                                    'options', opt);
    filename_g = strcat(paste0(data.loc, 'NNMF/nnmf_genes_k20_min0_Pat-'),num2str(pid(p)),'_defaultParams.csv');
    filename_c = strcat(paste0(data.loc, 'NNMF/nnmf_cells_k20_min0_Pat-'),num2str(pid(p)),'__defaultParams.csv');
    csvwrite(filename_g,Factors_g)
    csvwrite(filename_c,transpose(Factors_c))   
end



