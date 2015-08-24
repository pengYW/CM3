% **************************************************************** %
% ***********            Data Input:                    ********** %
% **************************************************************** %
% 1.) logpvec1,   vector of -log10 p values of the focal trait
% 2.) zvec1,      vector of z score values of the focal trait
% 3.) hvec,       vector of heterozygosity of SNPs
% 4.) LDmat1,     Binary Matlab sparse matrix for pairwise LD r2>=0.8
% 5.) LDmat2,     Binary Matlab sparse matrix for pairwise LD r2>=0.1
% 6.) M,          matrix of enrichment factors
% 7.) zmat,       matrix of z score values for each sub-study(by column)
% 8.) Neffvec,   vector of effective sample sizes for each sub-study

%%% load('example_data.mat');
% Input Parameters:
% **************************************************************** %
% ***********            Input Parameters:              ********** %
% **************************************************************** %
opts_struct = struct(...
    'nrep', 50,...       % number of random draws in disc/repl experiment
    'ntop', 100000,...    % number of top SNPs in empirical repl figure
    'pthresh', 0.05,...   % replication threshold
    'pruneflag', true,... % wheter performing random pruning
    'logpthresh', 3,...   % threhold of -log10 p for RES
    'nbins', 5,...       % N bins in stratifying SNPs use RES
    'dflist', [0.3 0.4],... % fractions in construct disc/repl empirical dataset
    'plotflag', true,...  % whether plot in each step to debuging
    'low_RES', 1,...      % minimum percentils for interopolating RES
    'high_RES',99.99,...  % maximum percentils for interopolating RES
    'n_RES',1000,...      % number of points for interopolating RES
    'low_QQ', 0,...       % minmum -log10 p for contruct conditional Q-Q
    'high_QQ', 300,...    % maximum -log10 p for contruct conditional Q-Q
    'n_QQ', 10000,...     % number of points betwen min/max for cond Q-Q
    'low_z', -10.0,...    % minimum value for interpolating z scores
    'high_z', 10.0,...    % maximum value for interpolating z scores
    'n_z', 201,...        % number of interpolating points for z scores
    'n_lookup', 101,...   % number of interpolating points in lookup table
    'maxiter', 10000,...  % number of iteration in fitting to empirical data
    'TolFun', 1e-2,...    % funtional tolerance for fitting to empirical data
    'TolX', 1e-3);        % termination tolerance for fitting to empirical data
% *** predict enrichment scores ** %
logpvec1_pruned = fast_prune(logpvec1, LDmat1);
[yvec_pred, yvec] = predicted_RES(logpvec1_pruned, M, opts_struct);
% *** stratify SNPs  ************* %
[iimat, qqmat, hv] = stratify_RES(logpvec1, yvec, yvec_pred, opts_struct);

% *** construct empirical discovery/replication dataset  ************* %
Neff = sum(Neffvec); nsamp = size(zmat, 2); studyindlist_disc = [1:nsamp];
zvals_disc = linspace(opts_struct.low_z, opts_struct.high_z, opts_struct.n_z);
datastructmat = construct_disc_repl_data(zmat, Neffvec, iimat,...
    LDmat1, zvals_disc, opts_struct, false, M);

% *** estimate parameter values  ************* %
% initial guess
pi1=0.0512; sigd=0.0248; sige=0.0078; sig0=1.0078; x0 = [pi1 sigd sige sig0];
xfits = fit_mixture(datastructmat, hvec, zvals_disc, x0, opts_struct);
% *** estimate posteriori z-score, z-score squared, replication probablity ** %
plotstruct = construct_empi_fitted_data(xfits, datastructmat, zvals_disc,...
        hvec, opts_struct);

% * predict z-score, z-score squared, replication probablity for future sample %
% % repvec_pred, predicted replication probablity for future sample
% % zvec_pred, predicted z score for future sample
% % z2vec_pred, predicted z-score^2 for future sample
[repvec_pred,zvec_pred,z2vec_pred,~] = z2z_lookuptables(...
        xfits(1:opts_struct.nbins), zvals_disc, hvec, Neff,...
        Neff, zvec1, yvec_pred, iimat(:,1:opts_struct.nbins),...
        linspace(1,opts_struct.nbins, opts_struct.n_lookup),opts_struct);

% *** make plots *********%
legname = cell(opts_struct.nbins+1, 1);
for i = 1:length(legname), legname{i} = sprintf('Bin%d', i); end
legname{end} = 'All SNPs';
plot_zz(plotstruct,zvals_disc, opts_struct.dflist, legname, true)
plot_z2z2(plotstruct,zvals_disc, opts_struct.dflist, legname, true)
plot_repl(plotstruct,zvals_disc,  opts_struct.dflist,  legname, true)

% *** empirical replication rate by resampling experiment *********%
for dfrac=opts_struct.dflist
    [repratemat, repratemat_cm] = construct_empi_replrate_data(zmat,Neffvec,...
          hvec, xfits, iimat, yvec_pred, zvals_disc, dfrac, LDmat2,...
          opts_struct);
    plot([mean(repratemat, 1)' mean(repratemat_cm, 1)'], 1:opts_struct.ntop);
    legend({'sorted by p' 'sorted by repl rate'});
    xlim([0.5 1]); ylim([0 310]); set(gca, 'Xdir', 'reverse');
end
