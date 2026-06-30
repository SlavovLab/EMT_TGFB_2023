% Scripts modified from EMT_2
% Forked on Sept 9th
% Goal: Load DART-ID updated data and compute correlations without imputaion

%% Section 1: Data import and wrangling 
% =========================================

Path = [path 'Collaborators/EmiliLab/EMT.TGFB_scProt_1/001-SCPipeline_Output/'];

[dat_L, txt] = xlsread( [Path 'EpiToMesen.TGFB.nPoP_trial1_1PercDartFDRTMTBulkDIA.WallE_unimputed.xlsx']);
Cell_Ids = txt(1,1:end);
ProtIDs_L = txt(2:end,1);

[~, IDs] = xlsread( [Path 'cellIDToTimepoint.xlsx']);
ind.d0 = IDs( ismember( IDs(:,2), 'd0' ), 1 );
ind.d3 = IDs( ismember( IDs(:,2), 'd3' ), 1 );
ind.d9 = IDs( ismember( IDs(:,2), 'd9' ), 1 );


ins = zeros( size(Cell_Ids) );
ins ( ismember( Cell_Ids, ind.d0 ) ) = 1;
ins ( ismember( Cell_Ids, ind.d3 ) ) = 2;
ins ( ismember( Cell_Ids, ind.d9 ) ) = 3;

%% Filter proteins to a subset with less missing data

in1 = sum(~isnan(dat_L),2) > 30;
in2 = sum( ~isnan(dat_L(:, ins==1 )), 2) > 4 & ...
      sum( ~isnan(dat_L(:, ins==2 )), 2) > 4 & ...
      sum( ~isnan(dat_L(:, ins==3 )), 2) > 4;
in3 = (in1 + in2)>1 ;

fprintf( 'Selected Proteins: %d\n\n', sum( in3 ) );
ProtIDs = ProtIDs_L(in3,:);
dat_nan = dat_L(in3,:);


%% Without IMPUTAION comute pairwise protein correlations for days 0, 3, and 9
dat_Lin = dat_L(in3,:);
r1 = pearNaN ( dat_Lin(:, ins == 1 )'  );
r2 = pearNaN ( dat_Lin(:, ins == 2 )'  );
r3 = pearNaN ( dat_Lin(:, ins == 3 )'  );


%% oOad SwissProt Annotation 
load ( [path 'bin/Human_ONLYsp_Swiss_2_GeneName'] );
[Lia,Locb] = ismember( ProtIDs,  HumanONLYspSwiss2GeneName(:, 1 ) );
%%



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
Clu = clu( Rd(in,in), 3); rb(1);
%%
Kmeans = 1;

if Kmeans == 1
    In = In(Idx); 
else
    In = In(Clu);
end


Lim = 1;
subplot(2,3,1)
imagesc( r1(In,In)); 
set( gca, 'CLim', [-Lim Lim] );
set( gca, 'Xtick', [], 'Ytick', [] );
h(1) = title( 'Day 0' );


subplot(2,3,2)
imagesc( r2(In,In)); 
set( gca, 'CLim', [-Lim Lim] );
set( gca, 'Xtick', [], 'Ytick', [] );
h(2) = title( 'Day 3' );

subplot(2,3,3)
imagesc( r3(In,In)); hb(1) = rb(Lim); 
set( gca, 'Xtick', [], 'Ytick', [] );
h(3) = title( 'Day 9' );


% Differences 
Lim = 0.7;
subplot(2,3,4)
imagesc( r2(In,In) - r1(In,In) );
set( gca, 'Xtick', [], 'Ytick', [] );
h(4) = title( 'Day 3 - Day 0' );

subplot(2,3,5)
imagesc( r2(In,In) - r1(In,In) ); 
set( gca, 'Xtick', [], 'Ytick', [] );
h(5) = title( 'Day 9 - Day 3' );

subplot(2,3,6)
imagesc( r3(In,In) - r1(In,In)  ); hb(2) = rb(Lim); 
set( gca, 'Xtick', [], 'Ytick', [] );
h(6) = title( 'Day 9 - Day 0' );


set(h,  'FontSize', 18 );
set(hb, 'FontSize', 12 );

set( hb(1), 'Ytick', [-1 0 1] );
set( hb(2), 'Ytick', [-0.7 0 0.7] );

%%
pdf( [path 'NS/PDFs/EMT_Global_correlation-Patterns_Kmeans'], [12 7], 1 );




%% Analyze protein clusters

%[Lia,Locb] = ismember( ProtIDs2,  HumanONLYspSwiss2GeneName(:, 1 ) );

fprf( ProtIDs2(Clu), 3, [path 'NS/txt/Prot_CorrClusters.txt']  );
fprf( ProtIDs2(flip(Clu,2)), 3, [path 'NS/txt/Prot_CorrClusters_flip.txt']  );




%% Cluster the matrix of correlation differences using K-means

Rd2 = Rd(in,in);
[idx,C] = kmeans( [Rd2 r1(In,In) r3(In,In)], 3); %

Idx = [find(idx==1); find(idx==2);  find(idx==3) ];

imagesc( Rd2(Idx,Idx) );  rb(1);

%%
fprf( ProtIDs2, 3, [path 'NS/txt/b/Prot_CorrClusters_Background.txt']  );
fprf( ProtIDs2(idx==1), 3, [path 'NS/txt/b/Prot_CorrCluster-1.txt']  );
fprf( ProtIDs2(idx==2), 3, [path 'NS/txt/b//Prot_CorrCluster-2.txt']  );
fprf( ProtIDs2(idx==3), 3, [path 'NS/txt/b/Prot_CorrCluster-3.txt']  );


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
ylim( [0.05 0.42]);
grid minor
set(gca, 'Xtick', 1:3, 'XtickLabel', {'Day 0', 'Day 3', 'Day 9'}, 'FontSize', 32);
set(gca, 'Ytick', 0.1: 0.1 :0.4);
h(1) = ylabel( 'Mean correlation');
h(2) = legend( 'c2  ', 'c1   ', 'c3');
%set(h(2), 'Orientation', 'horizontal' );
set(h, 'FontSize', 40 );
%%
pdf( [path 'NS/PDFs/Mean-cluster-correlation'], [9 6], 1 );



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
pdf( [path 'NS/PDFs/Mean-cluster-abundance'], [9 6], 1 );


%% Display correlation vs mean protein abundance 
clear h, close all
r_mean_corr_1 = pear2i( aMean', cMean' )

[r_mean_corr_2, p_val] = corr( aMean(:), cMean(:) )

plot( cMean, aMean, 'o', 'MarkerSize', 16, 'Linewidth', 5 );

xlim( [0.05 0.42]);
ylim( [-0.2 0.2]);
grid minor
set(gca, 'FontSize', 32);
set(gca, 'Ytick', -0.2: 0.1 :0.2 );
h(1) = ylabel( 'Mean abundance');
h(2) = xlabel( 'Mean correlation');
h(3) = legend( 'c2  ', 'c1   ', 'c3');
h(4) = text( 0.1, 0.15, sprintf( '$r = %1.2f$', r_mean_corr_2 ) );
set(h(4), 'interpreter', 'latex' );
set(h, 'FontSize', 40 );
%%
pdf( [path 'NS/PDFs/Mean-clustercorr-vs--abundance'], [9 6], 1 );
















