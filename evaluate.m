%% initialization
load(params.input.qlist_matname, 'query_imgnames_all');
densePV_matname = fullfile(params.output.dir, 'densePV_top10_shortlist.mat');
load(densePV_matname, 'ImgList');

if exist(params.evaluation.dir, 'dir') ~= 7
    mkdir(params.evaluation.dir)
end

%% visual evaluation
if exist(params.evaluation.query_vs_synth.dir, 'dir') ~= 7
    mkdir(params.evaluation.query_vs_synth.dir);
end

for i=1:size(query_imgnames_all,2)
    queryName = query_imgnames_all{i};
    queryImage = imread(fullfile(params.data.dir, params.data.q.dir, queryName));
    
    fun = @(x) strcmp(ImgList(x).queryname,queryName);
    tf = arrayfun(fun, 1:numel(ImgList));
    ImgListRecord = ImgList(find(tf));
    cutoutPath = ImgListRecord.topNname{1};
    synthName = strsplit(cutoutPath, '/');
    synthName = synthName{end};
    synthName = synthName(1:end-size(params.data.db.cutout.imgformat,2));
    synthPath = fullfile(params.output.synth.dir, queryName, [synthName, params.output.synth.matformat]);
    load(synthPath, 'RGBpersp');
    numRows = size(queryImage,1);
    numCols = size(queryImage,2);
    
    synthImage = RGBpersp;
    if isempty(synthImage)
        synthImage = zeros(numRows, numCols, 3, 'uint8');
    else
        synthImage = imresize(synthImage, [numRows numCols]);
    end

    imshowpair(queryImage, synthImage, 'montage');
    queryName = strrep(queryName, 'JPG', 'jpg');
    saveas(gcf, fullfile(params.evaluation.query_vs_synth.dir, queryName));
end

%% quantitative results
nQueries = size(query_imgnames_all,2);
errors = struct();
retrievedPoses = struct();

rFixes = permn([0.0, 90.0, -90.0, 180.0], 3);
initialDirections = [1.0, 0.0, 0.0; -1.0, 0.0, 0.0;
                        0.0, 1.0, 0.0; 0.0, -1.0, 0.0;
                        0.0, 0.0, 1.0; 0.0, 0.0, -1.0];
rFixOrders = ['XYZ'; 'ZYX'];
Rot1s = [1;2;3;4];
nCandidates = size(rFixes,1)*2*size(initialDirections,1)*2*size(rFixOrders,1)*2*size(Rot1s,1);
fprintf('Considering %d candidates\n', nCandidates);
lowestAngularError = inf;
j1_=1;j2_=6;j3_=1;j4_=2; % approximate checkpoint
restoreCheckpoint = true;
for j1=1:size(rFixes,1)
for j2=1:size(rFixes,1) 
for j3=1:size(initialDirections,1)
for j4=1:size(initialDirections,1)
for j5=1:size(rFixOrders,1)
for j6=1:size(rFixOrders,1)
for j7=1:size(Rot1s,1)
    if restoreCheckpoint
        restoreCheckpoint = false;
        j1=j1_;j2=j2_;j3=j3_;j4=j4_;
    end
queries = [1,2,4,5,11];
% for i=1:nQueries
errorSum = 0.0;
for ii=1:size(queries,2)
    i = queries(1,ii);
    queryName = query_imgnames_all{i};
    queryImage = imread(fullfile(params.data.dir, params.data.q.dir, queryName));
    
    fun = @(x) strcmp(ImgList(x).queryname,queryName);
    tf = arrayfun(fun, 1:numel(ImgList));
    ImgListRecord = ImgList(tf);
    
    cutoutPath = ImgListRecord.topNname{1};
    cutoutPath = strsplit(cutoutPath, '/');
    spaceName = cutoutPath{1};
    sweepId = cutoutPath{2};
    if contains(spaceName, 'DUC')
        spaceType = 'DUC';
    else
        spaceType = 'cse';
    end
    transPath = fullfile(params.data.dir, params.data.db.trans.dir, spaceName, 'transformations', ...
                            [spaceType, '_trans_', sweepId, '.txt']);
    P1 = load_WUSTL_transformation(transPath);
    R1 = P1(1:3,1:3);
    T1 = P1(1:3,4);
    
    P2 = ImgListRecord.P{1};
    R2 = P2(1:3,1:3); % in fact this is K*R - are you sure? viz demo

    rFix1 = rFixes(j1,:);
    rFix2 = rFixes(j2,:);
    initialDirection1 = initialDirections(j3,:)';
    initialDirection2 = initialDirections(j4,:)';
    Rfix1 = rotationMatrix(deg2rad(rFix1), rFixOrders(j5,:));
    Rfix2 = rotationMatrix(deg2rad(rFix2), rFixOrders(j6,:));
    Rot1 = Rot1s(j7);
    if Rot1 == 1
        Rot1 = R1;
    elseif Rot1 == 2
        Rot1 = R2;
    elseif Rot1 == 3
        Rot1 = R1*R2;
    elseif Rot1 == 4
        Rot1 = R2*R1;
    end
    
    if any(isnan(P2(:)))
        T = nan(3,1);
        orientation = nan(3,1);
    else
        T2 = P2(1:3,4);
        P = P2*P1;
        T = -inv(P2(1:3,1:3))*P2(1:3,4);
        
        orientation = Rfix1 * (Rot1 * initialDirection1);
    end
    
    queryId = strsplit(queryName, '.');
    queryId = queryId{1};
    posePath = fullfile(params.data.queryRefPoses.dir, [queryId, '.pose.txt']);
    if exist(posePath, 'file') ~= 2
        errors(i).queryId = queryId;
        errors(i).translation = -1;
        errors(i).orientation = -1;
        continue;
    end
    referenceP = load_CIIRC_transformation(posePath);
    referenceP = referenceP(1:3,1:4);
    referenceR = referenceP(1:3,1:3);

    referenceT = referenceP(1:3,4);
    referenceOrientation = Rfix2 * (referenceR * initialDirection2);
    %fprintf('orient: %g\n', orientation)
    %fprintf('refOrient: %g\n', referenceOrientation)
    %fprintf('T: %g\n', T);
    %fprintf('refT: %g\n', referenceT);
    
    errors(i).queryId = queryId;
    errors(i).translation = norm(T - referenceT);
    errors(i).orientation = atan2d(norm(cross(orientation,referenceOrientation)),dot(orientation,referenceOrientation));
%     fprintf('query %s errors: t: %g, orient: %g\n', queryId, errors(i).translation, errors(i).orientation);
    errorSum = errorSum + errors(i).orientation;
end
    avgError = errorSum / size(queries,2);
    if avgError < lowestAngularError
        lowestAngularError = avgError;
        fprintf('found a better solution with config: %d,%d,%d,%d,%d,%d,%d\n', j1,j2,j3,j4,j5,j6,j7);
        for ii=1:size(queries,2)
            i = queries(1,ii);
            fprintf('query %s errors: t: %g, orient: %g\n', queryId, errors(i).translation, errors(i).orientation);
        end
    end
end
end
end
end
end
end
end

% errors
errorsTable = struct2table(errors);
errors = table2struct(sortrows(errorsTable, 'queryId'));
errorsFile = fopen(params.evaluation.errors.path, 'w');
fprintf(errorsFile, 'id,translation,orientation\n');
for i=1:nQueries
    fprintf(errorsFile, '%s,%0.2f,%0.2f\n', errors(i).queryId, errors(i).translation, errors(i).orientation);
end
fclose(errorsFile);

%% summary
summaryFile = fopen(params.evaluation.summary.path, 'w');
thresholds = [[0.25 10], [0.5 10], [1 10]];
scores = zeros(1, size(thresholds,2)/2);
fprintf(summaryFile, 'Conditions: ');
nQueriesWithKnownReferencePoses = 0;
for i=1:2:size(thresholds,2)
    if i > 1
        fprintf(summaryFile, ' / ');
    end
    fprintf(summaryFile, '(%g [m], %g [deg])', thresholds(i), thresholds(i+1));
    
    count = 0;
    for j=1:size(errors,1)
        if errors(j).translation == -1
            continue
        end
        nQueriesWithKnownReferencePoses = nQueriesWithKnownReferencePoses + 1;
        if errors(j).translation < thresholds(i) && errors(j).orientation < thresholds(i+1)
            count = count + 1;
        end
    end

    scores((i-1)/2+1) = count / nQueriesWithKnownReferencePoses * 100;
end
fprintf(summaryFile, '\n');
for i=1:size(scores,2)
    if i > 1
        fprintf(summaryFile, ' / ');
    end
    fprintf(summaryFile, '%g [%%]', scores(i));
end
fprintf(summaryFile, '\n');
fclose(summaryFile);