% Run through all clusters
Get_all_cluster_cut = [];
for p = 1:29
    %% Denis Schapiro - Bodenmiller Lab - University of Zurich
    
    % Load SVCA output and visualize connection betweeen cell-cell
    % interactions, cell density and cell amount using the histoCAT dataset:
    % https://www.nature.com/articles/nmeth.4391
    
    %% (0) Clear work space
    %clear all
    %clc
    
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
    
    % Cluster_number for images:
    cluster_number = p;
    
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
    Position_PG = [];
    
    % Get average+std amount neighbors and percent touching
    for k=Remove_position'
        
        % Where are the channels
        % Find all channels with the correspoding name
        idx_neighbors = strfind(gates{k,3},'Number_Neighbors');
        idx_touching = strfind(gates{k,3},'Percent_Touching');
        idx_area = strfind(gates{k,3},'Area');
        idx_PG =  strfind(gates{k,3},'Phenograph');
        
        % Select the position (only one! HARDCODED)
        Position_neighbors = find(~cellfun(@isempty,idx_neighbors));
        Position_touching = find(~cellfun(@isempty,idx_touching));
        Position_area = find(~cellfun(@isempty,idx_area));
        
        % Get the PG cluster's
        Position_PG = find(~cellfun(@isempty,idx_PG));
        Get_all_cluster = sessionData(gates{k,2},Position_PG);
        Get_cluster_7(1,k) = numel(find(Get_all_cluster==cluster_number));
        
        % Calculate mean and std for amount neighbors
        Neighbor_results_mean(1,k) = mean(sessionData(gates{k,2},Position_neighbors));
        Neighbor_results_std(1,k) = std(sessionData(gates{k,2},Position_neighbors));
        
        % Calculate mean and std for percent touching
        Touching_results_mean(1,k) = mean(sessionData(gates{k,2},Position_touching));
        Touching_results_std(1,k) = std(sessionData(gates{k,2},Position_touching));
        
        % Calculate mean and std for area
        Area_results_mean(1,k) = mean(sessionData(gates{k,2},Position_area));
        Area_results_std(1,k) = std(sessionData(gates{k,2},Position_area));
        
        Amount_cells(1,k)= size(gates{k,2},2);
    end
    
    % Remove zeros (since some images could be missing (SVCA vs. histoCAT)
    Neighbor_results_mean(Neighbor_results_mean==0)=[];
    Neighbor_results_std(Neighbor_results_std==0)=[];
    Touching_results_mean(Touching_results_mean==0)=[];
    Touching_results_std(Touching_results_std==0)=[];
    Area_results_mean(Area_results_mean==0)=[];
    Area_results_std(Area_results_std==0)=[];
    
    
    Amount_cells(Amount_cells==0)=[];
    
    %% (8) Visualize results - Play around!
    %     % Visualize and calculate
    %     neighorsVScellcellall= [Neighbor_results_mean;table2array(Results_mean)];
    %     scatter(neighorsVScellcellall(1,:),neighorsVScellcellall(2,:),[],);
    %     title('Neighbors vs. SVCA Cell-Cell signature')
    %     ylabel('Cell-Cell interactions explained - average over all markers')
    %     xlabel('Average amount of neighbors')
    %
    % figure()
    % touchingVScellcellall= [Touching_results_mean;table2array(Results_mean)];
    % scatter(touchingVScellcellall(1,:),touchingVScellcellall(2,:));
    % title('Touching vs. SVCA Cell-Cell signature')
    % ylabel('Cell-Cell interactions explained - average over all markers')
    % xlabel('Average amount of percent touching')
    %
    % figure()
    % amountVScellcellall= [Amount_cells;table2array(Results_mean)];
    % scatter(amountVScellcellall(1,:),amountVScellcellall(2,:));
    % title('Amount cells vs. SVCA Cell-Cell signature')
    % ylabel('Cell-Cell interactions explained - average over all markers')
    % xlabel('Amount of cells')
    %
    % figure()
    % areaVScellcellall= [Area_results_mean;table2array(Results_mean)];
    % scatter(areaVScellcellall(1,:),areaVScellcellall(2,:));
    % title('Area vs. SVCA Cell-Cell signature')
    % ylabel('Cell-Cell interactions explained - average over all markers')
    % xlabel('Average Area')
    %
    % figure()
    % neighorsamountVScellcellall = [(Amount_cells./Neighbor_results_mean);table2array(Results_mean)];
    % scatter(neighorsamountVScellcellall(1,:),neighorsamountVScellcellall(2,:));
    % title('Cell for each present neighbor vs. SVCA Cell-Cell signature')
    % ylabel('Cell-Cell interactions explained - average over all markers')
    % xlabel('Cells for each present neighbor')
    %
    % figure()
    % neighorsVScellcelltop5 = [Neighbor_results_mean;Save_mean_top5];
    % scatter(neighorsVScellcelltop5(1,:),neighorsVScellcelltop5(2,:));
    % title('Neighbors vs. TOP5 SVCA Cell-Cell signature')
    % ylabel('Cell-Cell interactions explained - average over TOP5 markers')
    % xlabel('Average amount of neighbors')
    %
    % figure()
    % touchingVScellcelltop5 = [Touching_results_mean;Save_mean_top5];
    % scatter(touchingVScellcelltop5(1,:),touchingVScellcelltop5(2,:));
    % title('Touching vs. TOP5 SVCA Cell-Cell signature')
    % ylabel('Cell-Cell interactions explained - average over TOP5 markers')
    % xlabel('Average amount of percent touching')
    %
    % figure()
    % amountVScellcelltop5 = [Amount_cells;Save_mean_top5];
    % scatter(amountVScellcelltop5(1,:),amountVScellcelltop5(2,:));
    % title('Amount cells vs. TOP5 SVCA Cell-Cell signature')
    % ylabel('Cell-Cell interactions explained - average over TOP5 markers')
    % xlabel('Amount of cells')
    %
    % figure()
    % areaVScellcelltop5= [Area_results_mean;Save_mean_top5];
    % scatter(areaVScellcelltop5(1,:),areaVScellcelltop5(2,:));
    % title('Area vs. TOP5 SVCA Cell-Cell signature')
    % ylabel('Cell-Cell interactions explained - average over TOP5 markers')
    % xlabel('Average Area')
    %
    % figure()
    % neighorsamountVScellcellall = [(Amount_cells./Neighbor_results_mean);Save_mean_top5];
    % scatter(neighorsamountVScellcellall(1,:),neighorsamountVScellcellall(2,:));
    % title('Cell for each present neighbor vs. SVCA Cell-Cell signature')
    % ylabel('Cell-Cell interactions explained - average over TOP5 markers')
    % xlabel('Cells for each present neighbor')
    %
    % hold on
    %
    % p = polyfit(neighorsVScellcelltop5(1,:),neighorsVScellcelltop5(2,:),1);
    % f = polyval(p,neighorsVScellcelltop5(1,:));
    % plot(neighorsVScellcelltop5(1,:),neighorsVScellcelltop5(2,:),'o',...
    %     neighorsVScellcelltop5(1,:),f,'-')
    
    %% (9) Check for cluster 7 in PhenoGraph and its distribiution compared to the SVCA signature
    
    % Remove not used images
    Get_cluster_7_cut=Get_cluster_7(:,Remove_position);
    Get_all_cluster_cut = [Get_all_cluster_cut;Get_cluster_7_cut];
    
    %     % Creat figure for amount
    %     h = figure()
    %     PG7VScellcellall = [Get_cluster_7_cut;Save_mean_top5];
    %     scatter(PG7VScellcellall(1,:),PG7VScellcellall(2,:));
    %     title(strcat('Cell cluster ', num2str(cluster_number), ' vs. SVCA Cell-Cell signature'));
    %     ylabel('Cell-Cell interactions explained - average over TOP5 markers')
    %     xlabel(strcat('Amount cell cluster',num2str(cluster_number)))
    %     % Calculate correlations
    %     [Amount_all_R,Amount_all_P] = corrcoef(PG7VScellcellall(1,:),PG7VScellcellall(2,:));
    %     % Plot R value and P value on figure
    %     text(0.8, 0.9, strcat('R value: ',num2str(Amount_all_R(1,2))), ...
    %         'Units', 'normalized', ...
    %         'HorizontalAlignment', 'left', ...
    %         'VerticalAlignment', 'top');
    %
    %     text(0.8, 0.8, strcat('P value: ',num2str(Amount_all_P(1,2))), ...
    %         'Units', 'normalized', ...
    %         'HorizontalAlignment', 'left', ...
    %         'VerticalAlignment', 'top');
    %
    %     % Save figure
    %     savefig(h,strcat('PGCluster',num2str(cluster_number),'.fig'));
    %     print(strcat('PGCluster',num2str(cluster_number)),'-dtiffn');
    %
    %     % Create figure for percentage
    %     g = figure()
    %     PG7_Per_VScellcellall = [Get_cluster_7_cut./Amount_cells;Save_mean_top5];
    %     scatter(PG7_Per_VScellcellall(1,:),PG7_Per_VScellcellall(2,:));
    %     title(strcat('Percentage cell cluster ', num2str(cluster_number), ' vs. SVCA Cell-Cell signature'));
    %     ylabel('Cell-Cell interactions explained - average over TOP5 markers')
    %     xlabel(strcat('Percentage cell cluster',num2str(cluster_number)))
    %     % Calculate correlations
    %     [Percentage_all_R,Percentage_all_P] = corrcoef(PG7_Per_VScellcellall(1,:),PG7_Per_VScellcellall(2,:));
    %     % Plot R value and P value on figure
    %     text(0.8, 0.9, strcat('R value: ',num2str(Percentage_all_R(1,2))), ...
    %         'Units', 'normalized', ...
    %         'HorizontalAlignment', 'left', ...
    %         'VerticalAlignment', 'top');
    %
    %     text(0.8, 0.8, strcat('P value: ',num2str(Percentage_all_P(1,2))), ...
    %         'Units', 'normalized', ...
    %         'HorizontalAlignment', 'left', ...
    %         'VerticalAlignment', 'top');
    %     % Save figure
    %     savefig(g,strcat('Percentage_PGCluster',num2str(cluster_number),'.fig'));
    %     print(strcat('Percentage_PGCluster',num2str(cluster_number)),'-dtiffn');
    %
    %         %% (8) Visualize results - Play around!
    %     % Visualize and calculate
    %     neighorsVScellcellall= [Neighbor_results_mean;table2array(Results_mean)];
    %     scatter(neighorsVScellcellall(1,:),neighorsVScellcellall(2,:),[],PG7_Per_VScellcellall(1,:));
    %         % Calculate correlations
    %     [R,P] = corrcoef(neighorsVScellcellall(1,:),neighorsVScellcellall(2,:));
    %     % Plot R value and P value on figure
    %     text(0.6, 0.9, strcat('R value: ',num2str(R(1,2))), ...
    %         'Units', 'normalized', ...
    %         'HorizontalAlignment', 'left', ...
    %         'VerticalAlignment', 'top');
    %     text(0.6, 0.8, strcat('P value: ',num2str(P(1,2))), ...
    %         'Units', 'normalized', ...
    %         'HorizontalAlignment', 'left', ...
    %         'VerticalAlignment', 'top');
    %     title('Neighbors vs. SVCA Cell-Cell signature')
    %     ylabel('Cell-Cell interactions explained - average over all markers')
    %     xlabel('Average amount of neighbors')
end

%% (10) Save as CSV for Damien
% Save all readouts (Careful - Order HARDCODED!!!)
All_Readouts = [table2array(Results_mean);table2array(Results_std);...
    Neighbor_results_mean;Neighbor_results_std;...
    Touching_results_mean;Touching_results_std;...
    Amount_cells;Area_results_mean(Remove_position);...
    Area_results_std(Remove_position);...
    Save_mean_top5;Save_std_top5;Get_all_cluster_cut];

% Create string for row
Row_names = {'Results_mean','Results_std','Neighbor_results_mean',...
    'Neighbor_results_std','Touching_results_mean','Touching_results_std',...
    'Amount_cells','Area_results_mean','Area_results_std','Save_mean_top5',...
    'Save_std_top5',gates{53:end,1}};

% Create table
CSV_Output = array2table(All_Readouts,'RowNames',Row_names,'VariableNames',gates(Remove_position,1))

% Create CSV file
writetable(CSV_Output,'CSV_Output.csv','Delimiter',',','WriteRowNames',true);

%% (11) Load fcs files from cluster of interest
% Cluster of interest is 8 and 14
% Load via fcs files
[dat_8,hdr_8] = fca_readfcs...
    ('/Users/denis/Desktop/Damien_Revision/custom_gates_0/custom_gates_0/Phenograph5671971535Cluster_8.fcs');

[dat_14,hdr_14] = fca_readfcs...
    ('/Users/denis/Desktop/Damien_Revision/custom_gates_0/custom_gates_0/Phenograph5671971535Cluster_14.fcs');

% Load full session
load('/Users/denis/Desktop/Damien_Revision/custom_gates_0/Full_session.mat');

% Sort rows
sorted_sessionData = sortrows(sessionData,[1,2]);
sorted_dat_8 = sortrows(dat_8,[1,2]);
sorted_dat_14 = sortrows(dat_14,[1,2]);

% Find same image ids in sessionData
[C,position_in_session_8,ib] = intersect...
    (sorted_sessionData(:,1:2),sorted_dat_8(:,1:2),'rows');

[C,position_in_session_14,ib] = intersect...
    (sorted_sessionData(:,1:2),sorted_dat_14(:,1:2),'rows');

% Extract information
sessionData_PG8 = sessionData(position_in_session_8,:);
sessionData_PG14 = sessionData(position_in_session_14,:);

% Get biggest gate - HARDCODED!!!
Table_PG8 = array2table(sessionData_PG8,'VariableNames',gates{33,3});
Table_PG14 = array2table(sessionData_PG14,'VariableNames',gates{33,3});
Table_all = array2table(sessionData,'VariableNames',gates{33,3});

% Create CSV files
writetable(Table_PG8,'Table_PG8.csv','Delimiter',',','WriteRowNames',true);
writetable(Table_PG14,'Table_PG14.csv','Delimiter',',','WriteRowNames',true);
writetable(Table_all,'Table_all.csv','Delimiter',',','WriteRowNames',true);

x1 = Table_PG8.Cell_CD68Nd146Di;
x2 = Table_PG14.Cell_CD68Nd146Di;
x3 = Table_all.Cell_CD68Nd146Di;

x = [x1; x2; x3];
g = [zeros(size(x1,1), 1); ones(size(x2,1), 1); 2*ones(size(x3,1), 1)];
boxplot(x, g)
