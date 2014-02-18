%% Initialization - don't do this while CS229 folder still in directory
% Will lead to function overloading

% Add all of the external libraries in the directory
currentFolder = pwd;
addpath(genpath(currentFolder));

%% Data Pre-Processing - Data Conversion (fcs -> BioData)

% Create array of the listed fcs files in the directory './[root]'
root = '../CancerData/';
% fcsFiles = {'H1_min5_s0.15_m10_debar1_NoDrug_PVO4.fcs',...
%             'H2_min5_s0.10_m10_debar1_NoDrug_PVO4.fcs',...
%             'H3_min3_s0.10_m10_debar1_NoDrug_PVO4.fcs',...
%             'H4_min5_s0.20_m10_debar1_NoDrug_PVO4.fcs',...
%             'H5_min3_s0.20_m10_debar1_NoDrug_PVO4.fcs'...            
%             'H6_min5_s0.15_m10_debar1_NoDrug_PVO4.fcs'...  
%             'H7_min5_s0.25_m10_debar1_NoDrug_PVO4.fcs'};
% legendTitles = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7'};
fcsFiles = {
'H3_min3_s0.10_m10_debar1_NoDrug_PVO4.fcs',...
'H4_min5_s0.20_m10_debar1_NoDrug_PVO4.fcs',...
'H5_min3_s0.20_m10_debar1_NoDrug_PVO4.fcs',...
            };
legendTitles = {'H3', 'H4','H5'};
figName = 'Healthy Cells - stim comparison';

numRandTrainExPerFile = 400; % 400 seems optimal for tsne 
markerType = 'all'; %Options - 'all', 'surface';
xAxisName = ['#' num2str(numRandTrainExPerFile) ' cells/file, markers=' num2str(markerType)];

bioData_array = [];
fcsFiles = strcat(root, fcsFiles);
for i = 1:length(fcsFiles)
    bioData_array = [bioData_array, BioData(fcsFiles{i})]; % #ok<AGROW>
end

%% Data Pre-Processing - Data Selection ({BioData} -> BioData -> [mxn])

% Create single BioData object out of desired data
DesiredCells = BioData.merge(bioData_array, numRandTrainExPerFile);

% Obtain data matrix and pre-process with arcsinh
switch markerType
    case 'all'
        [data_stack, markers] = DesiredCells.getAllMarkerData();
%         data_stack = DesiredCells.data;
%         markers = DesiredCells.column_headings(DesiredCells.protein_column_index:end);
    case 'surface'
        [data_stack, markers] = DesiredCells.getSurfaceMarkerData();
    otherwise 
        error('Unrecognized Marker Type');
end
% if(useSurfaceProteinsOnly)
%     data_stack = DesiredCells.getSurfaceProteinData();
% else
%     data_stack = DesiredCells.data;
% end
data_stack = asinh(data_stack/5);

% Get data chunk indices - indices to the chunks of data that form the
% child object DesiredCells. Used for plotting
mergeGroups = DesiredCells.mergeGroups;

%% Heat Plot for stims

idx = DesiredCells.mergeGroups;
N = length(unique(idx));
matrix = zeros(N, size(data_stack,2));
for i = 1:N
    rows = data_stack(idx == i,:);
    numrows = sum(idx == i);
    matrix(i,:) = sum(rows,1)/numrows;
end

figure;
HeatMap(matrix, 'ColumnLabels', markers, 'RowLabels', legendTitles, 'ColorMap', 'redgreencmap');
% bar3(matrix)

%% PCA %%

%Run Algorithm
display('Running PCA on Data')
tic;
[coeff,score_PCA,latent] = princomp(data_stack);
PCA_time = toc

% Plot
figure;
CLR = [];
SYM = [];
SIZ = 20;
gscatter(score_PCA(:,1), score_PCA(:,2), mergeGroups, CLR, SYM, SIZ);
title(['PCA: ' figName]);
legend(legendTitles);
xlabel(xAxisName);
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

% Plot
figure;
CLR = [];
SYM = [];
SIZ = 20;
gscatter(tsne_output(:,1), tsne_output(:,2), mergeGroups, CLR, SYM, SIZ);
title(['tSNE: ' figName]);
legend(legendTitles);
xlabel(xAxisName);
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

% Plot
figure;
CLR = [];
SYM = [];
SIZ = 20;
gscatter(ssne_output(:,1), ssne_output(:,2), mergeGroups, CLR, SYM, SIZ);
title(['sSNE: ' figName]);
legend(legendTitles);
xlabel(xAxisName);
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

%Plot
figure;
CLR = [];
SYM = [];
SIZ = 20;
gscatter(ee_output(:,1), ee_output(:,2), mergeGroups, CLR, SYM, SIZ);
title(['EE: ' figName]);
legend(legendTitles);
xlabel(xAxisName);
makeFiguresPretty;
