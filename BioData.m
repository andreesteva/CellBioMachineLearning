% BioData
% this class reads a .fcs file and stores its data in a matrix
% Purpose: abstract away the logistical handling of the fcs data and be
%   able to merge multiple sets of data as needed
% 
% BioData is a simpler and more higher-level version of CellData, intended
% to work for a variety of data sets (including Karen's latest Cancer Data
% set as of January 2014
%
% FCS file data format
%   contains a large matrix where rows correspond to individual cells and
%   the columns correspond to some variable. This variable is usually a
%   protein, in which case the matrix entry would be the counts of that
%   protein in the processed cell. Other variables include time, etc.
%
%   contains a header struct, which labels the columns of the matrix as
%   well as contains data like the name of the cell type

classdef BioData < handle

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PROPERTIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    properties(Constant = true)
        fcs_reader = @fca_readfcs;
        surface_protein_tag = 'CD';  
        marker_colIdx_start = 3;
        marker_colIdx_end = 42;

    end
    
    properties
        data = []; % [m x n] matrix read from the .fcs file   
        column_headings = {}; % string of the names of the columns
        header = {}; % data struct as read in by the reader
    end
    
    properties(SetAccess = private)
        mergeGroups = [] % For merged objects, these are indices (for the rows ) to the parent from whence that row came. Primarily for use with gscatter
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Static = true)
        
        % MERGE %
        % Creates a child object from all the parents in pointer_array
        % pointer_array is an array of pointers to BioData objects
        %
        % throws an error if the column_headings are not identical (case
        % sensitive)
        % the child has the concatenated data of the parents (in no
        %   particular order) along the 1st dimension of the matrix (cell
        %   rows)
        % the child has en empty header
        function obj = merge(pointer_array, numRandDataPts)
            if((length(pointer_array) ~= 1) && ~isequal(pointer_array.column_headings))
                error('pointers being merged must have identical column headers');
            end
            
            if(~exist('numRandDataPts', 'var'));
                numRandDataPts = inf;
            end
            
            % Create new object and add variable values
            obj = BioData();
            obj.column_headings = pointer_array(1).column_headings;             
            for i = 1:length(pointer_array)
                % Add Data
                ad = pointer_array(i).data;
                r = randperm(size(ad,1));
                r = r(1:min(numRandDataPts, size(ad,1)));
                ad = ad(r,:);
                obj.data = [obj.data; ad];
                
                % Update mergeGroup
                m = ones(size(ad, 1), 1) * i;
                obj.mergeGroups = [obj.mergeGroups; m];                
            end
        end
    end
    
    methods
        % Class Constructor - BioData
        function this = BioData(fcs_filename)
            %Empty constructor (MATLAB does not support overloading fncs            
            if(nargin == 0)
                return;
            end
            
            % Parse data and header
            [this.data, this.header] = BioData.fcs_reader(fcs_filename);
            this.column_headings = {this.header.par.name2};
        end
        
        % getSurfaceProteinData %
        % returns a matrix of values corresponding to the surface protein
        % data only. These are proteins with [surface_protein_tag] in their
        % name
        function [matrix, markers] = getSurfaceMarkerData(this)
            proteinNames = this.column_headings;
            tf = strfind(proteinNames, BioData.surface_protein_tag);
            idxSP = ~cellfun('isempty', tf);
            matrix = this.data(:, idxSP);  
            markers = proteinNames(idxSP);
        end
        
        function [matrix, markers] = getAllMarkerData(this)
            markers = this.column_headings(BioData.marker_colIdx_start:BioData.marker_colIdx_end);
            matrix = this.data(:, BioData.marker_colIdx_start:BioData.marker_colIdx_end);
        end
        
        
    end
    
    
end
