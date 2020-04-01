function [ params ] = setup_project_ht_WUSTL
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

params = struct();
%WUSTL dataset
env = environment();
if strcmp(env, 'ciirc')
    params.data.dir = '/home/lucivpav/InLoc_dataset';
    params.data.netvlad.dir = '/home/lucivpav/NetVLAD';
elseif strcmp(env, 'cmp')
    params.data.dir = '/datagrid/personal/lucivpav/InLoc_dataset';
    params.data.netvlad.dir = '/datagrid/personal/lucivpav/NetVLAD';
elseif strcmp(env, 'laptop')
    params.data.dir = '/Volumes/GoogleDrive/Můj disk/ARTwin/InLoc_dataset';
    params.data.netvlad.dir = '/Volumes/GoogleDrive/Můj disk/ARTwin/InLocCIIRC_dataset/NetVLAD';
end
params.data.netvlad.pretrained = fullfile(params.data.netvlad.dir, 'vd16_pitts30k_conv5_3_vlad_preL2_intra_white.mat');
%database
params.data.db.dir = 'database';
params.data.db.subsets_name = {'DUC1', 'DUC2', 'CSE3', 'CSE4', 'CSE5'};
params.data.db.subsets_header = {'DUC_', 'DUC_', 'cse_', 'cse_', 'cse_'};
%%scan
params.data.db.scan.dir = fullfile(params.data.db.dir, 'scans');
params.data.db.scan.header = strcat(params.data.db.subsets_header, 'scan_');
params.data.db.scan.imgformat = '.ptx.png';
params.data.db.scan.matformat = '.ptx.mat';
%%cutouts
params.data.db.cutout.dir = fullfile(params.data.db.dir, 'cutouts');
params.data.db.cutout.header = strcat(params.data.db.subsets_header, 'cutout_');
params.data.db.cutout.imgformat = '.jpg';
params.data.db.cutout.matformat = '.mat';
%%alignments
params.data.db.trans.dir = fullfile(params.data.db.dir, 'alignments');
params.data.db.trans.header = strcat(params.data.db.subsets_header, 'trans_');
%query
params.data.q.dir = 'query/iphone7';
params.data.q.imgformat = '.JPG';
params.data.q.fl = 4032*28/36;
params.data.queryRefPoses.dir = fullfile(params.data.dir, 'queryRefPoses');

%input
params.input.dir = fullfile(params.data.dir, 'inputs');
params.input.dblist_matname = fullfile(params.input.dir, 'cutout_imgnames_all.mat');%string cell containing cutout image names
params.input.qlist_matname = fullfile(params.input.dir, 'query_imgnames_all.mat');%string cell containing query image names
params.input.score_matname = fullfile(params.input.dir, 'scores.mat');%retrieval score matrix
params.input.feature.dir = fullfile(params.input.dir, 'features');
params.input.feature.db_matformat = '.features.dense.mat';
params.input.feature.q_matformat = '.features.dense.mat';
params.input.feature.db_sps_matformat = '.features.sparse.mat';
params.input.feature.q_sps_matformat = '.features.sparse.mat';


%output
params.output.dir = fullfile(params.data.dir, 'outputs');
params.output.gv_dense.dir = fullfile(params.output.dir, 'gv_dense');%dense matching results (directory)
params.output.gv_dense.matformat = '.gv_dense.mat';%dense matching results (file extention)
params.output.gv_sparse.dir = fullfile(params.output.dir, 'gv_sparse');%sparse matching results (directory)
params.output.gv_sparse.matformat = '.gv_sparse.mat';%sparse matching results (file extention)

params.output.pnp_dense_inlier.dir = fullfile(params.output.dir, 'PnP_dense_inlier');%PnP results (directory)
params.output.pnp_dense.matformat = '.pnp_dense_inlier.mat';%PnP results (file extention)
params.output.pnp_sparse_inlier.dir = fullfile(params.output.dir, 'PnP_sparse_inlier');%PnP results (directory)
params.output.pnp_sparse_inlier.matformat = '.pnp_sparse_inlier.mat';%PnP results (file extention)

params.output.pnp_sparse_origin.dir = fullfile(params.output.dir, 'PnP_sparse_origin');%PnP results (directory)
params.output.pnp_sparse_origin.matformat = '.pnp_sparse_origin.mat';%PnP results (file extention)

params.output.synth.dir = fullfile(params.output.dir, 'synthesized');%View synthesis results (directory)
params.output.synth.matformat = '.synth.mat';%View synthesis results (file extention)


%groundtruth - TODO: where is the data?
params.gt.dir = 'Refposes';
params.gt.matname = 'DUC_refposes_all.mat';

% evaluation
params.evaluation.dir = fullfile(params.data.dir, 'evaluation');
params.evaluation.query_vs_synth.dir = fullfile(params.evaluation.dir, 'queryVsSynth');
params.evaluation.errors.path = fullfile(params.evaluation.dir, 'errors.csv');
params.evaluation.summary.path = fullfile(params.evaluation.dir, 'summary.txt');

end

