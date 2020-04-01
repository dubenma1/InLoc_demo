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
for i=1:nQueries
    %i = 278; % 41 % TODO: remove
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
    if any(isnan(P2(:)))
        T = nan(3,1);
        orientation = nan(3,1);
    else
        R2 = P2(1:3,1:3); % in fact this is K*R - are you sure? viz demo
        T2 = P2(1:3,4);

        P = P2*P1;
        T = P2(:,4);
        T = -inv(P(1:3,1:3))*P(1:3,4);
        %T = P1(1:3,4) + P2(1:3,4);
        %T = P(1:3,4);
        initialDirection = [-1.0; 0.0; 0.0];
        rFix = [0.0, 0.0, 0.0];
        Rfix = rotationMatrix(deg2rad(rFix), 'XYZ');
        orientation = Rfix * (R2 * initialDirection);
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

    referenceT = -inv(referenceR)*referenceP(:,4);
    initialDirection = [-1.0; 0.0; 0.0];
    rFix = [-90.0, 0.0, 0.0];
    Rfix = rotationMatrix(deg2rad(rFix), 'XYZ');
    referenceOrientation = Rfix * (referenceR * initialDirection);
    %fprintf('orient: %g\n', orientation)
    %fprintf('refOrient: %g\n', referenceOrientation)
    %fprintf('T: %g\n', T);
    %fprintf('refT: %g\n', referenceT);
    
    errors(i).queryId = queryId;
    errors(i).translation = norm(T - referenceT);
    errors(i).orientation = atan2d(norm(cross(orientation,referenceOrientation)),dot(orientation,referenceOrientation));
    fprintf('errors: t: %g, orient: %g\n', errors(i).translation, errors(i).orientation);
    xxx = 42;
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
for i=1:2:size(thresholds,2)
    if i > 1
        fprintf(summaryFile, ' / ');
    end
    fprintf(summaryFile, '(%g [m], %g [deg])', thresholds(i), thresholds(i+1));
    
    count = 0;
    for j=1:size(errors,1)
        if errors(j).translation < thresholds(i) && errors(j).orientation < thresholds(i+1)
            count = count + 1;
        end
    end

    scores((i-1)/2+1) = count / size(errors,1) * 100;
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