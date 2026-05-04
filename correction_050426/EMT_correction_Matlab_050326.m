% EMT correction in Matlab
% 050326

%% Section 1: Data import and wrangling 
% =========================================

clearvars

% raw mixed content (text + numbers)
raw = readcell('EpiToMesen.TGFB.nPoP_trial1_ProtByCellMatrix_NSThreshDART_medIntCrNorm_unimputed.txt');

% numeric data (NA -> NaN)
num = readmatrix('EpiToMesen.TGFB.nPoP_trial1_ProtByCellMatrix_NSThreshDART_medIntCrNorm_unimputed.txt');

% Assign as you requested:
Cell_Ids = raw(1, 1:end);       % second row, all columns
ProtIDs_L = raw(2:end, 1);       % first column, from row 2 downward
dat_L = num(:, 2:end);

IDs = readcell('cellIDToTimepoint.txt');

% Convert relevant columns to string arrays
col2 = string( IDs(:,2) );        % timepoint column as string
ids_col1 = string( IDs(:,1) );    % ID column as string
cell_ids_str = string( Cell_Ids );% column names as string

% build groups
ind.d0 = ids_col1( col2 == "d0" );
ind.d3 = ids_col1( col2 == "d3" );
ind.d9 = ids_col1( col2 == "d9" );

% assign numeric labels
ins = zeros( size(cell_ids_str) );
ins( ismember(cell_ids_str, ind.d0) ) = 1;
ins( ismember(cell_ids_str, ind.d3) ) = 2;
ins( ismember(cell_ids_str, ind.d9) ) = 3;

%% Filter proteins to a subset with less missing data

in1 = sum(~isnan(dat_L),2) > 30;
in2 = sum( ~isnan(dat_L(:, ins==1 )), 2) > 4 & ...
      sum( ~isnan(dat_L(:, ins==2 )), 2) > 4 & ...
      sum( ~isnan(dat_L(:, ins==3 )), 2) > 4;
in3 = (in1 + in2)>1 ;

fprintf( 'Selected Proteins: %d\n\n', sum( in3 ) );
ProtIDs = ProtIDs_L(in3,:);
dat_nan = dat_L(in3,:);

%% Without IMPUTAION commute pairwise protein correlations for days 0, 3, and 9

dat_Lin = dat_L(in3,:);
r1 = pearNaN ( dat_Lin(:, ins == 1 )'  );
r2 = pearNaN ( dat_Lin(:, ins == 2 )'  );
r3 = pearNaN ( dat_Lin(:, ins == 3 )'  );

%% oOad SwissProt Annotation 

load ('Human_ONLYsp_Swiss_2_GeneName.mat');
[Lia,Locb] = ismember( ProtIDs,  HumanONLYspSwiss2GeneName(:, 1 ) );

%% Section 2: Identify and plot correlation CLUSTERS
% ====================================================

rd = r3 - r1;

% Select proteins with complete pairwise correlations
Mean = mean ( isnan( rd ) );
in_NAN =  find( Mean == 0 ); 
Rd = rd( in_NAN, in_NAN );


% Select correlations without missing data and large difference 
SumR = sum( Rd.^2 );
in = SumR > prctile( SumR, 50 ); 
fprintf( 'Number of selected proteins: %d \n\n', sum(in) );
ProtIDs2 = ProtIDs( in_NAN(in) );
In = in_NAN(in);

% Cluster the matrix of correlation differences 
%Clu = clu( Rd(in,in), 3); rb(1);

%% Cluster matrix of correlation differences using K-means (Figure 2)

rng(1, 'twister');

X = [r1(In,In), r2(In,In), r3(In,In)];
[idx, C] = kmeans(X, 3, 'Replicates', 1);

Idx = [find(idx==1); find(idx==2); find(idx==3)];

Kmeans = 1;
if Kmeans == 1
    Inx = In(Idx);
else
    Inx = In(Clu);
end

%% Heatmaps

figure;

Lim1 = 0.7;
Lim2 = 0.5;

ncol = 200;
blue_white_red = [
    linspace(0,1,ncol/2)', linspace(0,1,ncol/2)', ones(ncol/2,1);
    ones(ncol/2,1), linspace(1,0,ncol/2)', linspace(1,0,ncol/2)'
    ];

colormap(blue_white_red);

subplot(2,3,1)
imagesc(r1(Inx, Inx));
set(gca, 'CLim', [-Lim1 Lim1], 'Xtick', [], 'Ytick', []);
h(1) = title('Day 0');

subplot(2,3,2)
imagesc(r2(Inx, Inx));
set(gca, 'CLim', [-Lim1 Lim1], 'Xtick', [], 'Ytick', []);
h(2) = title('Day 3');

subplot(2,3,3)
imagesc(r3(Inx, Inx));
set(gca, 'CLim', [-Lim1 Lim1], 'Xtick', [], 'Ytick', []);
h(3) = title('Day 9');
hb(1) = colorbar;
set(hb(1), 'Ytick', [-Lim1 0 Lim1]);

subplot(2,3,4)
imagesc(r2(Inx, Inx) - r1(Inx, Inx));
set(gca, 'CLim', [-Lim2 Lim2], 'Xtick', [], 'Ytick', []);
h(4) = title('Day 3 - Day 0');

subplot(2,3,5)
imagesc(r3(Inx, Inx) - r2(Inx, Inx));
set(gca, 'CLim', [-Lim2 Lim2], 'Xtick', [], 'Ytick', []);
h(5) = title('Day 9 - Day 3');

subplot(2,3,6)
imagesc(r3(Inx, Inx) - r1(Inx, Inx));
set(gca, 'CLim', [-Lim2 Lim2], 'Xtick', [], 'Ytick', []);
h(6) = title('Day 9 - Day 0');
hb(2) = colorbar;
set(hb(2), 'Ytick', [-Lim2 0 Lim2]);

set(h, 'FontSize', 18);
set(hb, 'FontSize', 12);

%%
pdf('EMT_Correct/EMT_Global_correlation-Patterns_Kmeans-EMT_6', [12 7], 1 );

%% Figure 3a

close all, clear h

ci = In(idx == 1);

x = squareform(r1(ci,ci) - diag(diag(r1(ci,ci))), 'tovector');
y = squareform(r3(ci,ci) - diag(diag(r3(ci,ci))), 'tovector');

plot(x, y, 'ko'); hold on
plot([-0.6 0.9], [-0.6 0.9], 'r', 'LineWidth', 3)

grid minor
xlim([-0.6 0.9]);
ylim([-0.6 0.9]);

set(gca, 'FontSize', 28);
set(gca, 'Ytick', -0.6:0.3:0.9);
set(gca, 'Xtick', -0.6:0.3:0.9);

h(1) = ylabel('Day 9');
h(2) = xlabel('Day 0');
h(3) = text(-0.5, 0.75, 'c2 correlations');
h(4) = text(0.3, -0.4, sprintf('$r = %1.2f$', corr(x(:), y(:))));
set(h(4), 'interpreter', 'latex');
set(h, 'FontSize', 40);

%%
pdf('EMT_Correct/Corr-scaling', [9 6], 1 );

%% Display mean correlations for each cluster across time 
clear h, close all
cMean =  zeros( 3 ); 
for i = 1:3
    ci = In(idx==i);
    
    cMean(1,i) = mean( squareform( r1(ci,ci)-diag(diag(r1(ci,ci))), 'tovector') ); 
    cMean(2,i) = mean( squareform( r2(ci,ci)-diag(diag(r2(ci,ci))), 'tovector') ); 
    cMean(3,i) = mean( squareform( r3(ci,ci)-diag(diag(r3(ci,ci))), 'tovector') ); 
end

plot( 1:3, cMean, 'o-', 'MarkerSize', 16, 'Linewidth', 5 );
xlim( [0.5 3.9]);
ylim( [0.00 0.30]);
grid minor
set(gca, 'Xtick', 1:3, 'XtickLabel', {'Day 0', 'Day 3', 'Day 9'}, 'FontSize', 32);
set(gca, 'Ytick', 0:0.1:0.3);
h(1) = ylabel( 'Mean correlation');
h(2) = legend( 'c2  ', 'c1   ', 'c3');
h(2) = legend('c2', 'c1', 'c3', 'Location', 'northeast');  % or 'southeast'
%set(h(2), 'Orientation', 'horizontal' );
set(h, 'FontSize', 40 );
%%
pdf('EMT_Correct/Mean-cluster-correlation', [9 6], 1 );



%% Display mean protein abundance for each cluster across time 
clear h, close all
aMean =  zeros( 3 ); 
for i = 1:3
    ci = In(idx==i);
    
    aMean(1,i) = mean(nanmean( dat_Lin(ci, ins == 1)  )); 
    aMean(2,i) = mean(nanmean( dat_Lin(ci, ins == 2)  )); 
    aMean(3,i) = mean(nanmean( dat_Lin(ci, ins == 3)  )); 
end

plot( 1:3, aMean, 'o-', 'MarkerSize', 16, 'Linewidth', 5 );
xlim( [0.5 3.9]);
ylim( [-0.2 0.2]);
grid minor
set(gca, 'Xtick', 1:3, 'XtickLabel', {'Day 0', 'Day 3', 'Day 9'}, 'FontSize', 32);
set(gca, 'Ytick', -0.2: 0.1 :0.2 );
h(1) = ylabel( 'Mean abundance');
h(2) = legend( 'c2  ', 'c1   ', 'c3');
%set(h(2), 'Orientation', 'horizontal' );
set(h, 'FontSize', 40 );
%%
pdf('EMT_Correct/Mean-cluster-abundance', [9 6], 1 );


%% Display correlation vs mean protein abundance 
clear h, close all
r_mean_corr_1 = pear2i( aMean', cMean' )

[r_mean_corr_2, p_val] = corr( aMean(:), cMean(:) )

plot( cMean, aMean, 'o', 'MarkerSize', 16, 'Linewidth', 5 );

xlim( [0.0 0.35]);
ylim( [-0.2 0.15]);
grid minor
set(gca, 'FontSize', 32);
set(gca, 'Ytick', -0.2: 0.1 :0.2 );
h(1) = ylabel( 'Mean abundance');
h(2) = xlabel( 'Mean correlation');
h(3) = legend( 'c2  ', 'c1   ', 'c3');
h(4) = text(0.35, -0.18, sprintf('$r = %1.2f$', r_mean_corr_2));
set(h(4), 'interpreter', 'latex', 'HorizontalAlignment', 'right');
set(h(4), 'interpreter', 'latex' );
set(h, 'FontSize', 40 );
%%
pdf('EMT_Correct/Mean-clustercorr-vs--abundance', [9 6], 1 );

%% Analyze protein clusters (Table) 
fprf( ProtIDs2, 3, 'EMT_Correct/Prot_CorrClusters_Background.txt');
fprf( ProtIDs2(idx==1), 3, 'EMT_Correct/Prot_CorrCluster-1.txt');
fprf( ProtIDs2(idx==2), 3, 'EMT_Correct/Prot_CorrCluster-2.txt');
fprf( ProtIDs2(idx==3), 3, 'EMT_Correct/Prot_CorrCluster-3.txt');
