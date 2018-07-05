%% Denis Schapiro - Bodenmiller Lab - University of Zurich

% Load SVCA output and visualize connection betweeen cell-cell
% interactions, cell density and cell amount using the histoCAT dataset:
% https://www.nature.com/articles/nmeth.4391

%% (1) Load all necessary files
% Load all SVCA output files
Input_Folder = '/Users/denis/Desktop/Damien_Revision/signatures_per_image/';

% Load histoCAT paper data:
% (http://www.bodenmillerlab.org/histoCAT_downloads/Data_52_BreastCancerSamples.zip)
load('/Users/denis/Desktop/histoCAT_sessionResubmission.mat');

%% (2) Prepare for loading the SVCA signatures
% Get all image names in the order of the sessionData
Image_Names_all_full = gates(:,1);

% 52 images in the data set --> remove PhenoGraph names (HARDCODED!!!)
Image_number = 52;

% Clean variables before the for-loop
Split_Image_Names = {};
Cell_Cell_All = [];
Remove_position = [];

%% (3) Load SVCA signatures using a for-loop over all images
for i=1:Image_number
    
    % Load the individual SVCA signatures
    % Split string to just keep the image names
    Split_Image_Names_Temp = strsplit(Image_Names_all_full{i,:},'_');
    Split_Image_Names(i,1) = Split_Image_Names_Temp(1,1);
    Split_Image_Names_Temp = [];
    
    % Check if file exist
    if exist(fullfile(Input_Folder,['svca_signatures_',Split_Image_Names{i,1}]),...
            'file') == 2
        % if yes
        % Create table
        T_Temp= readtable(fullfile(Input_Folder,['svca_signatures_',Split_Image_Names{i,1}]),...
            'Delimiter',' ','ReadVariableNames',true,'ReadRowNames',true);
        
        % Flip the table to improve handling
        T_Temp_Array = table2array(T_Temp);
        YourNewTable = array2table(T_Temp_Array.');
        YourNewTable.Properties.RowNames = T_Temp.Properties.VariableNames;
        YourNewTable.Properties.VariableNames = T_Temp.Properties.RowNames;
        
        % Save the "env" variable representing cell-cell interactions
        % (HARDCODED!!!)
        Cell_Cell_All = [Cell_Cell_All, YourNewTable.env];
        
        % Save that the position in "Split_Image_Names" should be saved
        Remove_position(i,:) = i;
        
    else
        % Save that the position in "Split_Image_Names" should be removed
        Remove_position(i,:) = 0;
        
    end
end

%% (4) Remove not used images
% Remove zeros (if not all images used in comparison to the histoCAT dataset)
Remove_position(Remove_position==0)=[];

% Get all names
Split_Names_Cut = Split_Image_Names(Remove_position,1);

% Create table with only cell-cell interactions
Cell_Cell_Table = array2table(Cell_Cell_All,'RowNames',YourNewTable.Properties.RowNames,...
    'VariableNames',Split_Names_Cut');

%% (5)Create average + STD for all markers for cell-cell interactions across all images
% Get average+std columns
Results_mean = varfun(@mean,Cell_Cell_Table);
Results_std = varfun(@std,Cell_Cell_Table);

%% (6)Create average + STD for top 5 markers for cell-cell interactions across all images
% Prepare variables
Save_sorted = {};
Save_mean_top5 = [];
Save_std_top5 = [];

% Get top 5 markers (ToDo)
for m=1:size(Cell_Cell_Table,2)
    Sorted_row_all = sortrows(Cell_Cell_Table,m,'descend');
    Save_sorted{m} = Sorted_row_all(:,m);
    Save_mean_top5(1,m) = mean(table2array(Sorted_row_all(1:5,m)));
    Save_std_top5(1,m) = std(table2array(Sorted_row_all(1:5,m)));
end

%% (7) Get average + STD for the amount of cells, percent touching and amount neighbors 
% Create empty channels
Neighbor_results_mean = [];
Neighbor_results_std = [];
Touching_results_mean = [];
Touching_results_std = [];
Amount_cells = [];

% Get average+std amount neighbors and percent touching
for k=Remove_position'
    
    % Where are the channels
    % Find all channels with the correspoding name
    idx_neighbors = strfind(gates{k,3},'Number_Neighbors');
    idx_touching = strfind(gates{k,3},'Percent_Touching');
    % Select the position (only one! HARDCODED)
    Position_neighbors = find(~cellfun(@isempty,idx_neighbors));
    Position_touching = find(~cellfun(@isempty,idx_touching));
    
    % Calculate mean and std for amount neighbors
    Neighbor_results_mean(1,k) = mean(sessionData(gates{k,2},Position_neighbors));
    Neighbor_results_std(1,k) = std(sessionData(gates{k,2},Position_neighbors));
    % Calculate mean and std for amount neighbors
    Touching_results_mean(1,k) = mean(sessionData(gates{k,2},Position_touching));
    Touching_results_std(1,k) = std(sessionData(gates{k,2},Position_touching));
    
    Amount_cells(1,k)= size(gates{k,2},2);
end

% Remove zeros (since some images could be missing (SVCA vs. histoCAT)
Neighbor_results_mean(Neighbor_results_mean==0)=[];
Neighbor_results_std(Neighbor_results_std==0)=[];
Touching_results_mean(Touching_results_mean==0)=[];
Touching_results_std(Touching_results_std==0)=[];
Amount_cells(Amount_cells==0)=[];

%% (8) Visualize results - Play around!
% Visualize and calculate
neighorsVScellcellall= [Neighbor_results_mean;table2array(Results_mean)];
scatter(neighorsVScellcellall(1,:),neighorsVScellcellall(2,:));

figure()
touchingVScellcellall= [Touching_results_mean;table2array(Results_mean)];
scatter(touchingVScellcellall(1,:),touchingVScellcellall(2,:));

figure()
amountVScellcellall= [Amount_cells;table2array(Results_mean)];
scatter(amountVScellcellall(1,:),amountVScellcellall(2,:));

figure()
neighorsVScellcelltop5 = [Neighbor_results_mean;Save_mean_top5];
scatter(neighorsVScellcelltop5(1,:),neighorsVScellcelltop5(2,:));
hold on

p = polyfit(neighorsVScellcelltop5(1,:),neighorsVScellcelltop5(2,:),1);
f = polyval(p,neighorsVScellcelltop5(1,:));
plot(neighorsVScellcelltop5(1,:),neighorsVScellcelltop5(2,:),'o',...
    neighorsVScellcelltop5(1,:),f,'-')
