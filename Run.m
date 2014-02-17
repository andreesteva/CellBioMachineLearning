%% Initialization - don't do this while CS229 folder still in directory
% Will lead to function overloading

% Add all of the external libraries in the directory
currentFolder = pwd;
addpath(genpath(currentFolder));

%% Data Pre-Processing - Data Conversion (fcs -> BioData)

% Create array of the listed fcs files in the directory './[root]'
root = '../CancerData/';
fcsFiles = {'H1_min5_s0.15_m10_debar1_NoDrug_PVO4.fcs',...
            'H2_min5_s0.10_m10_debar1_NoDrug_PVO4.fcs',...
            'H3_min3_s0.10_m10_debar1_NoDrug_PVO4.fcs',...
            'H4_min5_s0.20_m10_debar1_NoDrug_PVO4.fcs',...
            'H5_min3_s0.20_m10_debar1_NoDrug_PVO4.fcs'...            
            'H6_min5_s0.15_m10_debar1_NoDrug_PVO4.fcs'...  
            'H7_min5_s0.25_m10_debar1_NoDrug_PVO4.fcs'};
legendTitles = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7'};
numRandTrainExPerFile = 400; % 400 seems optimal for tsne 
useSurfaceProteinsOnly = 1;

bioData_array = [];
fcsFiles = strcat(root, fcsFiles);
for i = 1:length(fcsFiles)
    bioData_array = [bioData_array, BioData(fcsFiles{i})]; % #ok<AGROW>
end

%% Data Pre-Processing - Data Selection ({BioData} -> BioData -> [mxn])

% Create single BioData object out of desired data
DesiredCells = BioData.merge(bioData_array, numRandTrainExPerFile);

% Obtain data matrix and pre-process with arcsinh
if(useSurfaceProteinsOnly)
    data_stack = DesiredCells.getSurfaceProteinData();
else
    data_stack = DesiredCells.data;
end
data_stack = asinh(data_stack/5);

% Get data chunk indices - indices to the chunks of data that form the
% child object DesiredCells. Used for plotting
mergeGroups = DesiredCells.mergeGroups;

%% PCA %%

%Run Algorithm
display('Running PCA on Data')
tic;
[coeff,score_PCA,latent] = princomp(data_stack);
PCA_time = toc

% User alters these variables
figName = 'PCA';
CLR = [];
SYM = [];
SIZ = 20;

% Plot
gscatter(score_PCA(:,1), score_PCA(:,2), mergeGroups, CLR, SYM, SIZ);
title(figName);
makeFiguresPretty;

%% t-SNE (Max's Code) &&

% Want dimensionality reduction to 2
dim = 2;

% stopping criteria: number of iterations is no more than 100, runtime is
% no more than 30 seconds, and the relative tolerance in the embedding is 
% no less than 1e-3. Taken from Max's tsne example demo_swissroll.m
opts.maxit = 400; opts.runtime = 900; opts.tol = 1e-3;
opts.X0 = 1e-5*randn(size(data_stack, 1), dim);

% Run algorithm
display('Running t-SNE (Max Code) on Data');
tic;
[tsne_output, E, A, T] = alg_tsne(data_stack, dim, opts);
tsne_time = toc

% The plot parameters
figName = 'tSNE';
CLR = [];
SYM = [];
SIZ = 20;

% Plot
gscatter(tsne_output(:,1), tsne_output(:,2), mergeGroups, CLR, SYM, SIZ);
title(figName);
makeFiguresPretty;

%% s-SNE %%

% Want dimensionality reduction to 2
dim = 2;

% stopping criteria: number of iterations is no more than 100, runtime is
% no more than 30 seconds, and the relative tolerance in the embedding is 
% no less than 1e-3. Taken from Max's tsne example demo_swissroll.m
opts.maxit = 600; opts.runtime = 900; opts.tol = 1e-3;
opts.X0 = 1e-5*randn(size(data_stack, 1), dim);

% Run algorithm
tic;
[ssne_output, E, A, T] = alg_ssne(data_stack, dim, opts);
ssne_time = toc

% The plot parameters
figName = 'sSNE';
CLR = [];
SYM = [];
SIZ = 20;

% Plot
gscatter(ssne_output(:,1), ssne_output(:,2), mergeGroups, CLR, SYM, SIZ);
title(figName);
makeFiguresPretty;

%% EE %%

% Want dimensionality reduction to 2
dim = 2;

% stopping criteria: number of iterations is no more than 100, runtime is
% no more than 30 seconds, and the relative tolerance in the embedding is 
% no less than 1e-3. Taken from Max's tsne example demo_swissroll.m
opts.maxit = 100; opts.runtime = 900; opts.tol = 1e-3;
opts.X0 = 1e-5*randn(size(data_stack, 1), dim);

% Run algorithm
tic;
[ee_output, E, A, T] = alg_ee(data_stack, dim, opts);
ee_time = toc

% The plot parameters
figName = 'EE';
CLR = [];
SYM = [];
SIZ = 20;

% Plot
gscatter(ee_output(:,1), ee_output(:,2), mergeGroups, CLR, SYM, SIZ);
title(figName);
makeFiguresPretty;
