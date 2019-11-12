%% Non-negative matrix factorisation
% as implemented in Regev HNSCC paper
% after selecting for ~ 8,000 genes based on Regev approach 
% from oligodendroglioma paper

% set folder name
folder = 'X';

% patient IDs
pid = [90121; 92047; 90130; 92492; 91526; 92507; 92460; 92489; 92074; 100399];

% set number of factors
k = 30;

opt = statset( 'MaxIter', 10, ...
                'TolFun', 1e-10, ...
                'TolX', 1e-10, ...
                'Display', 'iter');

% use mean-centered, log-transformed expression values
data = csvread(strcat('D:/Teresa/Colon/NNMF_new10/nnmf_input_matrix_min0_new10.csv'),1,1);
rng(1)
[Factors_g,Factors_c,D] = nnmf(data,k, ...
                                    'algorithm', 'als', ...
                                    'replicates', 1);
filename_g = 'D:/Teresa/Colon/NNMF_new10/nnmf_new10_genes_k30_min0_defaultParams.csv';
filename_c = 'D:/Teresa/Colon/NNMF_new10/nnmf_new10_cells_k30_min0_defaultParams.csv';    
csvwrite(filename_g,Factors_g)
csvwrite(filename_c,transpose(Factors_c))   


for p=1:length(pid)
    data = csvread(strcat('D:/Teresa/Colon/180226_NNMF/nnmf_input_matrix_min0_Pat-',num2str(pid(p)),'.csv'),1,1);
    data = data(1:100,1:100);
    rng(1)
    [Factors_g,Factors_c,D] = nnmf(data,k, ...
                                    'algorithm', 'mult', ...
                                    'replicates', 10, ...
                                    'options', opt);
    filename_g = strcat('D:/Teresa/Colon/180226_NNMF/nnmf_genes_k20_min0_Pat-',num2str(pid(p)),'_defaultParams.csv');
    filename_c = strcat('D:/Teresa/Colon/180226_NNMF/nnmf_cells_k20_min0_Pat-',num2str(pid(p)),'__defaultParams.csv');
    csvwrite(filename_g,Factors_g)
    csvwrite(filename_c,transpose(Factors_c))   
end



