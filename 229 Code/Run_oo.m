% Rewrite the Run file using the CellData classes
% Purpose:
%   create a simple and intuitive interface for test driving new algorithms
%   The gory details of working with any of the fcs data should be
%   abstracted away

% Formalize this as the control script
clear all; clc;
% Add all of the external libraries in the directory
currentFolder = pwd;
addpath(genpath(currentFolder));

%% %%%%% DATA PRE-PROCESSING - 229 Project%%%%%

% Script Variables 
cell_data_array_all = []; % array of pointers to all CellData objects created

% Create array of all fcs files as CellData objects
D=dir;
fnames = {D.name};
for i = 1:length(fnames)
    if(CellData.isReadable(fnames{i}))
        cell_data_array_all = [cell_data_array_all, CellData(fnames{i})]; %#ok<AGROW>
    end
end

%% DATA PROCESSING - 229 Project

%%%%% VARIABLE INITIALIZATION %%%%%

% Cell Categories - definitions per the biologists
CancerCells = Set({'CancerCells'});
StemCells = Set({'HSC', 'MPP', 'CMP', 'GMP', 'MEP'});
BCells = Set({'Plasma cell', 'Pre-B I', 'Pre-B II', 'Immature B', 'Mature CD38lo B', 'Mature CD38mid B'});
TCells = Set({'Mature CD4+ T', 'Mature CD8+ T', 'Naive CD4+ T', 'Naive CD8+ T'});
NK = Set({'NK'});
pDC = Set({'Plasmacytoid DC'});
Monocytes = Set({'CD11b- Monocyte', 'CD11bhi Monocyte', 'CD11bmid Monocyte'});

% User Variables
% whichCellTypes = pDC & NK & Monocytes & TCells & BCells & StemCells;
% % whichCellTypes = TCells & Monocytes;
% whichCellTypes = Monocytes;
whichCellTypes = StemCells & BCells & TCells & NK & pDC & NK & Monocytes;
% whichCellTypes = NK;
numRandTrainExPerFile = 400; % 400 seems optimal for tsne 
whichStimLevels = Set({'Basal'}); % Either 'Basal' or 'PV04', can contain both
useSurfaceProteinsOnly = false;

%%%%% DATA PARSING %%%%%

% Keep the CellData objects whose cell_type is contained in the
% whichCellTypes variable
removeIndicies = [];
for i = 1:length(cell_data_array_all)
    ct = cell_data_array_all(i).cell_types; % will be a single string since no CellData objects have seen merger
    st = cell_data_array_all(i).cell_stimulation_levels; % also a string
    if(~whichCellTypes.contains(ct) || ~whichStimLevels.contains(st))
        removeIndicies = [removeIndicies, i]; %#ok<AGROW>
    end
end
cell_data_array = cell_data_array_all;
cell_data_array(removeIndicies) = [];

%% Data Processing - Research

% Create array of the listed fcs files in the directory './[root]'
root = '../../CancerData/';
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

cell_data_array = [];
fcsFiles = strcat(root, fcsFiles);
for i = 1:length(fcsFiles)
    if(CellData.isReadable(fcsFiles{i}))
        cell_data_array = [cell_data_array, CellData(fcsFiles{i})]; % #ok<AGROW>
    end
end

%% Condense Desired .fcs files into a single CellData object
% Create single CellData object out of desired data
DesiredCells = CellData.merge(cell_data_array, numRandTrainExPerFile);

% Obtain data matrix and pre-process with arcsinh
if(useSurfaceProteinsOnly)
    data_stack = DesiredCells.getSurfaceProteinData();
else
    data_stack = DesiredCells.getProteinData();
end
data_stack = asinh(data_stack/5);

% Get data chunk indices - indices to the chunks of data that form the
% child object DesiredCells. Used for plotting
chunk_indices = DesiredCells.data_subset_indicies_for_merged_data;

% Color of data
colors= DesiredCells.getRelativeColorsByAbsoluteLabels(); % no argument returns colors based on the original labeling of the data


%% PCA %%

%Run Algorithm
display('Running PCA on Data')
tic
[coeff,score_PCA,latent] = princomp(data_stack);
PCA_time = toc

% User alters these variables
processed_data = score_PCA;
figName = 'PCA';

% CLUSTERING METRIC BASED ON AVERAGING
[label_err new_labels centroids]  = CentroidClusteringMetric(processed_data(:,1:2), DesiredCells);
plot_title = [figName ': N/file=' num2str(numRandTrainExPerFile) ', err=' num2str(label_err)];

% Plot the processed data
plotTrueLabelsWithCentroids(figName, processed_data, DesiredCells, centroids, numRandTrainExPerFile, plot_title)
plot_title = ['Relabel ' plot_title];
plotNewLabelsWithCentroids(figName, processed_data, DesiredCells, centroids, numRandTrainExPerFile, plot_title, new_labels)

%Save error
pca_label_err = label_err;

% Makes figures more presentable for papers
makeFiguresPretty;

%% t-SNE (Max's Code) %%
% dimensionality reduction to dim = 2
% see the file alg_tsne for more details

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
tsne_time = toc;

% User alters these variables
processed_data = tsne_output;
figName = 'tSNE';
plot_title = ['tSNE: iter #' num2str(length(E)), ', e=' num2str(E(end)) ', t=' num2str(T(end)), 'N/file=' num2str(numRandTrainExPerFile)];

% CLUSTERING METRIC BASED ON AVERAGING
[label_err new_labels centroids]  = CentroidClusteringMetric(processed_data(:,1:2), DesiredCells);
tsne_label_err = label_err;
plot_title = [figName ': err=' num2str(label_err)];

% Plot the processed data
plotTrueLabelsWithCentroids(figName, processed_data, DesiredCells, centroids, numRandTrainExPerFile, plot_title)
plot_title = ['Relabel ' plot_title];
plotNewLabelsWithCentroids(figName, processed_data, DesiredCells, centroids, numRandTrainExPerFile, plot_title, new_labels)

% Makes figures more presentable for papers
makeFiguresPretty;

% Find k-means clusters from reduced data from t-SNE
% [tsne_centroid_indices,tsne_centroid_locations,tsne_cluster_point_separation] = kmeans(tsne_output,15);

%% Plot t-SNE
chunk_indices = DesiredCells.data_subset_indicies_for_merged_data;
figure; hold on;
for i =1:length(fcsFiles)
    lb = chunk_indices(i);
    ub = chunk_indices(i+1)-1;
    scatter(processed_data(lb:ub,1), processed_data(lb:ub,2), 20, [0.9 0.9 0.9]/i, 'fill'); 
    scatter(processed_data(lb:ub,1), processed_data(lb:ub,2), 20, rand(1,3), 'fill'); 
    group(lb:ub) = i;
end
% gscatter(processed_data(lb:ub,1), processed_data(lb:ub,2), group)
legend(legendTitles);

%% s-SNE %%
% dimensionality reduction to dim = 2
% see the file alg_ssne for more details

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
ssne_time = toc;

% User alters these variables
processed_data = ssne_output;
figName = 'sSNE: ';
plot_title = ['sSNE: iter #' num2str(length(E)), ', e=' num2str(E(end)) ', t=' num2str(T(end)), 'N/file=' num2str(numRandTrainExPerFile)];

% CLUSTERING METRIC BASED ON AVERAGING
[label_err new_labels centroids]  = CentroidClusteringMetric(processed_data(:,1:2), DesiredCells);
ssne_label_err = label_err;
plot_title = [figName ': err=' num2str(label_err)];

% Plot the processed data
plotTrueLabelsWithCentroids(figName, processed_data, DesiredCells, centroids, numRandTrainExPerFile, plot_title)
plot_title = ['Relabel ' plot_title];
plotNewLabelsWithCentroids(figName, processed_data, DesiredCells, centroids, numRandTrainExPerFile, plot_title, new_labels)

% Makes figures more presentable for papers
makeFiguresPretty;

% Find k-means clusters from reduced data from s-SNE
% [ssne_centroid_indices,ssne_centroid_locations,ssne_cluster_point_separation] = kmeans(ssne_output,15);

%% EE %%
% dimensionality reduction to dim = 2
% see the file alg_ee for more details

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
ee_time = toc;

% % Plot results
% figure('name','EE (Maxs Code)'); hold on;
% for i = 1:whichCellTypes.length()
%     lb = chunk_indices(i);
%     ub = chunk_indices(i+1)-1;        
%     scatter(ee_output(lb:ub,1),ee_output(lb:ub,2), 20, colors(i,:));
%     title(['EE: iter #' num2str(length(E)), ', e=' num2str(E(end)),...
%        ', t=' num2str(T(end)), ', N/file=' num2str(numRandTrainExPerFile)]);   
% end
% legend(whichCellTypes.list)

% User alters these variables
processed_data = ee_output;
figName = 'EE';
plot_title = ['EE: iter #' num2str(length(E)), ', e=' num2str(E(end)) ', t=' num2str(T(end)), 'N/file=' num2str(numRandTrainExPerFile)];

% CLUSTERING METRIC BASED ON AVERAGING
[label_err new_labels centroids]  = CentroidClusteringMetric(processed_data(:,1:2), DesiredCells);
ee_label_err = label_err;
plot_title = [figName ': err=' num2str(label_err)];

% Plot the processed data
plotTrueLabelsWithCentroids(figName, processed_data, DesiredCells, centroids, numRandTrainExPerFile, plot_title)
plot_title = ['Relabel ' plot_title];
plotNewLabelsWithCentroids(figName, processed_data, DesiredCells, centroids, numRandTrainExPerFile, plot_title, new_labels)

% Makes figures more presentable for papers
makeFiguresPretty;

%% Merge select figures into 1 with subplots

% Find figures
figHandles_all = findobj('Type','figure'); % Get all
figHandles = [1 2]; % indexes to figure numbers
nrows = 1;
ncols = 2;

mergeFiguresIntoSubplot(nrows, ncols, figHandles);

% % Get subplot positions
% nfindex = max(figHandles_all) + 1;
% nf = figure(nfindex); % Create new figure
% sppos = []
% for i = 1:length(figHandles)
%     ax = subplot(nrows, ncols,i);
%     sppos = [sppos; get(ax, 'pos')];
%     cla(ax);
% end
% cla(nf);
% 
% % Copy figures into subplots of new figure
% new_splots = {};
% for i = 1:length(figHandles)
%     new_splots{end +1} = copyobj(get(figHandles(i), 'children'), nfindex);
% end
% for i = 1:length(figHandles)
%     set(new_splots{i}, 'pos', sppos(i,:));
% end


%% %%%%%%%%%%%%% ALGORITHMS FROM DR TOOLBOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold off;
drawnow


% Find k-means clusters from reduced data from EE
% [ee_centroid_indices,ee_centroid_locations,ee_cluster_point_separation] = kmeans(ee_output,15);


%% Merge select figures into 1 with subplots

% % Find figures
% figHandles_all = findobj('Type','figure'); % Get all
% figHandles = [1 2]; % indexes to figure numbers
% nrows = 1;
% ncols = 2;
% 
% % Get subplot positions
% nfindex = max(figHandles_all) + 1;
% figure(nfindex); % Create new figure
% sppos = []
% for i = 1:length(figHandles)
%     sppos = [sppos; get(subplot(nrows, ncols,i), 'pos')];
% end
% 
% % Copy figures into subplots of new figure
% new_splots = {};
% for i = 1:length(figHandles)
%     new_splots{end +1} = copyobj(get(figHandles(i), 'children'), nfindex);
% end
% for i = 1:length(figHandles)
%     set(new_splots{i}, 'pos', sppos(i,:));
% end

%% %%%%%%%%%%%%% ALGORITHMS FROM DR TOOLBOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Naive Linear Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Multidimensional Scaling

% % Run non-classical MDS on the data
% display ('Running non-classical MDS')
% tic;
% dissimilarities = pdist(zscore(dataStack));
% [score_ncMDS,stress] = mdscale(dissimilarities,2,'criterion','metricstress');
% mds_time = toc; 
% 
% % Plot the Figure for non-classical MDS
% display ('Plotting non-classical MDS Result')
% figure('name','ncMDS');
% hold on;
% for i = 1:length(whichCellTypes)
%     lb = scoreIndices(i)+1;
%     ub = scoreIndices(i+1);
%     if(plotIn3D)
%         display('plotting 3D');
%         scatter3(score_ncMDS(lb:ub,1),score_ncMDS(lb:ub,2),score_ncMDS(lb:ub,3), 20, colors(i,:));
%     else
%         % Plot scatter for PCA
%         scatter(score_ncMDS(lb:ub,1),score_ncMDS(lb:ub,2), 20, colors(i,:));
%         xlabel('p1');
%         ylabel('p2');
%         title(strcat('non-classical MDS ',expr_title));
%     end
% %     scatter(score((i-1)*400+1:400*i,1),score(i:400*i + 1,2),'filled','b');
% %     scatter(score((i-1)*400+1:400*i,1),score((i-1)*400+1:400*i,2), c(i),'filled');
% %     scatter(score(lb:ub,1),score(lb:ub,2), c(i),'filled', 'markersize', 10);
% end
% legend(whichCellTypes)
% hold off;
% drawnow
%% 

% Run ICA on the data - good for separation

%%%%%%%%%%%%%% Non - Linear Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ISOMAP Algorithm

% Run Isomap Algorithm on Data
display('Running Isomap on Data');
tic;
[score_isomap, mapping_isomap] = isomap(data_stack);
isomap_time = toc;
%ISOMAP Runs the Isomap algorithm
%
%   [mappedX, mapping] = isomap(X, no_dims, k); 


% User alters these variables
processed_data = score_isomap;
figName = 'ISOMAP';

% % % CLUSTERING METRIC BASED ON AVERAGING - WONT WORK WITH ISOMAP
% % [label_err new_labels centroids]  = CentroidClusteringMetric(processed_data(:,1:2), DesiredCells);
% plot_title = [figName ': N/file=' num2str(numRandTrainExPerFile) ', err=' num2str(label_err)];
% plot_title = [figName];
% 
%     % Plot the processed data
% plotTrueLabelsWithCentroids(figName, processed_data, DesiredCells, [], numRandTrainExPerFile, plot_title)
% plot_title = ['Relabel ' plot_title];
% plotNewLabelsWithCentroids(figName, processed_data, DesiredCells, [], numRandTrainExPerFile, plot_title, new_labels)

figure('name','Isomap');
hold on;
scatter(score_isomap(:,1),score_isomap(:,2), 60, 'fill');
title(['Isomap: All Cell Types']);

% % Plot the Figure for Isomap
% % This code is just incorrect - isomap gives you fewer points then in the
% % original data space, so you don't scatter per cell type
% plotIn3D = false;
% display('Plotting Isomap Result')
% figure('name','Isomap');
% hold on;
% for i = 1:whichCellTypes.length()
%     lb = chunk_indices(i);
%     ub = chunk_indices(i+1)-1;
%     if ub>size(score_isomap,1)
%         ub=size(score_isomap,1)
%     end
%     if(plotIn3D)
%         display('plotting 3D');
%         if(size(score_isomap,2)<3)
%             display('Cannot plot Isomap results in 3D - need more data - Run Isomap with no_dims of >=3');
%         else
%             scatter3(score_isomap(lb:ub,1),score_isomap(lb:ub,2),score_isomap(lb:ub,3), 20, colors(i,:));
%         end
%     else
%         % Plot scatter for Isomap
%         scatter(score_isomap(lb:ub,1),score_isomap(lb:ub,2), 20, colors(i,:), 'fill');
%         title(['Isomap: N/file=' num2str(numRandTrainExPerFile)]);
%     end
% end
% legend(whichCellTypes.list)
% hold off;
% drawnow

% Find k-means clusters from reduced data from Isomap
[isomap_centroid_indices,isomap_centroid_locations,isomap_cluster_point_separation] = kmeans(score_isomap,15);

%% Locally Linear Embedding
% Run LLE on Data
display('Running LLE on the data')
tic;
score_LLE = lle(data_stack);
lle_time = toc;
%LLE Runs the locally linear embedding algorithm
%
%   mappedX = lle(X, no_dims, k, eig_impl)

% Plot the Figure for LLE
display('Plotting LLE Result');
figure('name','LLE');
hold on;
for i = 1:whichCellTypes.length()
    lb = chunk_indices(i);
    ub = chunk_indices(i+1)-1;
    if ub>size(score_LLE,1)
        ub=size(score_LLE,1)
    end
    if(plotIn3D)
        display('plotting 3D');
        if(size(score_LLE,2)<3)
            display('Cannot plot SNE results in 3D - need more data - Run tSNE with no_dims of >=3');
        else
            scatter3(score_LLE(lb:ub,1),score_LLE(lb:ub,2),score_LLE(lb:ub,3), 20, colors(i,:));
        end
    else
        % Plot scatter for tSNE
        scatter(score_LLE(lb:ub,1),score_LLE(lb:ub,2), 20, colors(i,:));
        title(['LLE: N/file=' num2str(numRandTrainExPerFile)]);
    end
end
legend(whichCellTypes.list)
hold off;
drawnow

% Find k-means clusters from reduced data from Isomap
% [LLE_centroid_indices,LLE_centroid_locations,LLE_cluster_point_separation] = kmeans(score_LLE,15);

%%%%%%%%%%% SNE & t-SNE ALGORITHMS (take a while to converge) %%%%%%%%%%%
%% 

% % Run SNE on Data 
% display('Running SNE on Data')
% tic;
% score_SNE = sne(data_stack);
% slow_sne_time = toc;
% %SNE Implementation of Stochastic Neighbor Embedding
% %
% %   mappedX = sne(X, no_dims, perplexity)

% % Plot the Figure for SNE
% display('Plotting SNE Result')
% figure('name','SNE');
% hold on;
% for i = 1:whichCellTypes.length()
%     lb = chunk_indices(i);
%     ub = chunk_indices(i+1)-1;
%     if(plotIn3D)
%         display('plotting 3D');
%         if(size(score_SNE,2)<3)
%             display('Cannot plot SNE results in 3D - need more data - Run tSNE with no_dims of >=3');
%         else
%             scatter3(score_SNE(lb:ub,1),score_SNE(lb:ub,2),score_SNE(lb:ub,3), 20, colors(i,:));
%         end
%     else
%         % Plot scatter for tSNE
%         scatter(score_SNE(lb:ub,1),score_SNE(lb:ub,2), 20, colors(i,:));
%         title(['SNE: N/file=' num2str(numRandTrainExPerFile)]);
%     end
% end
% legend(whichCellTypes.list)
% hold off;
% drawnow
%% 

% % Run tSNE on the data 
% display('Running t-SNE on data')
% tic;
% score_tSNE = tsne(data_stack);
% slow_tsne_time = toc;
% %TSNE Performs symmetric t-SNE on dataset X
% %
% %   mappedX = tsne(X, labels, no_dims, initial_dims, perplexity)
% %   mappedX = tsne(X, labels, initial_solution, perplexity)

% % Plot the Figure for t-SNE
% display('Plotting t-SNE Result')
% figure('name','t-SNE');
% hold on;
% for i = 1:whichCellTypes.length()
%     lb = chunk_indices(i);
%     ub = chunk_indices(i+1)-1;
%     if(plotIn3D)
%         display('plotting 3D');
%         if(size(score_tSNE,2)<3)
%             display('Cannot plot tSNE results in 3D - need more data - Run tSNE with no_dims of >=3');
%         else
%             scatter3(score_tSNE(lb:ub,1),score_tSNE(lb:ub,2),score_tSNE(lb:ub,3), 20, colors(i,:));
%         end
%     else
%         % Plot scatter for tSNE
%         scatter(score_tSNE(lb:ub,1),score_tSNE(lb:ub,2), 20, colors(i,:));
%         title(['t-SNE: N/file=' num2str(numRandTrainExPerFile)]);
%     end
% end
% legend(whichCellTypes.list)
% hold off;
% drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%% Classification Algorithms %%%%%%%%%%%%%%%%%%%%%%%%

%% Liner Support Vector Machine

%%%%%%%%% MULTICLASS CLASSIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

actual_labels = DesiredCells.data_celltype_indices;
% Randomly permute rows before processing
permutation_mat = randperm(size(actual_labels,2));
actual_labels = actual_labels';
actual_labels = actual_labels(permutation_mat,1);
data_stack = data_stack(permutation_mat,:);
TrainingSet = data_stack(1:5000,:);
GroupTrain = actual_labels(1:5000);
TestSet = data_stack(5001:end,:);

u=unique(GroupTrain);
N=length(u);
if(N>2)
    disp('multi class problem');
    itr=1;
    classes=0;
    results = zeros(size(TestSet,1),1);
    while(itr<=length(u))
        %this while loop is the multiclass SVM Trick
        c1=(GroupTrain==u(itr));
        newClass=c1;
        options.MaxIter = 100000;
        svmStruct = svmtrain(TrainingSet,newClass,'showplot',true,'Options',options);
        classes = svmclassify(svmStruct,TestSet,'showplot',true)
        % Populate the results 
        for i=1:length(classes)
            if(classes(i)==1)
                if results(i)~=0
                    results(i) = u(itr);
                end
            end
        end
        itr=itr+1
    end
%clc;
itr=itr-1;
disp(itr)
end


%results = multisvm(TrainingSet, GroupTrain, TestSet); 

%%%%%% BINARY CLASSIFICATION - TRIAL UNTIL WE GET CANCER DATA %%%%%%%%%%%

% % Prepare the class labels
% % First prepare the labels for the healthy cells (label as zero)
% healthy_cell_labels = zeros(size(data_stack,1),1);
% % (TBD) Prepare the labels for the cancer cells
% % cancer_cell_labels
% % (TBD) Combine the data and the labels for healthy + cancer cells & randomize
% 
% % Run the algorithm on the data
% SVMStruct = svmtrain(data_stack,healthy_cell_labels);
