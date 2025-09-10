%% EOD Fish Tracking Across Minutes (cm & Hz)
% - Reads 1-min WAVs from a folder
% - Filters (300–2000 Hz), Hilbert envelope, event detection
% - Localizes each EOD event on a 2D grid via Gaussian spatial weighting
% - Clusters events per minute (DBSCAN) => fish per file
% - Builds amplitude fingerprints for clusters (8-ch pattern)
% - Tracks fish identity across consecutive minutes with a Kalman filter
%   using [x,y,vx,vy,f] and amplitude-aware association
%
% Requirements: Signal Processing Toolbox (filtfilt, filter, hilbert),
%               Statistics and Machine Learning Toolbox (dbscan).
% Optional:     munkres.m on path (else greedy assignment is used).

% ------------ USER: set default root folder here -------------
addpath('D:\KIT3'); addpath('D:\KIT');
rootDefault = 'D:\datosgrid2022\10 a 11_marzo 2023\';

clearvars -except rootDefault; close all; clc;

%% Select folder and enumerate .wav files
myDir = uigetdir(rootDefault, 'Select WAV folder containing 1-min WAV files');
if isequal(myDir,0), error('No folder selected.'); end
files = dir(fullfile(myDir, '*.wav'));
if isempty(files), error('No .wav files found in %s', myDir); end

%% Geometry: 8 sensors laid out on a 3x4 grid (cm)
% Sensor coordinates (x,y) in cm (edit if your array changes)
XYmeas = [
    0   0;   % ch1
    0  80;   % ch2
    40   0;   % ch3
    40  40;   % ch4
    40  80;   % ch5
    40 120;   % ch6
    80  40;   % ch7
    80 120];  % ch8
numChExpected = size(XYmeas,1);  % should be 8

% Dense interpolation grid for localization (cm)
xlimGrid = [0 80]; ylimGrid = [0 120]; step = 1;  % cm
[xg, yg] = meshgrid(xlimGrid(1):step:xlimGrid(2), ylimGrid(1):step:ylimGrid(2));
XYgrid = [xg(:) yg(:)];

% Spatial kernel (Gaussian) precomputed once (per geometry)
sigmaSpatial = 30;  % cm; try 10–30
Dists = pdist2(XYmeas, XYgrid);
Wproto = exp(-Dists.^2 / (2*sigmaSpatial^2));  % (8 x nGrid)

%% Detection parameters (tune if needed)
minEventsForFish = 200;   % min consolidated events in a file to consider fish present
minPkHeight      = 0.015; % envelope threshold (Hz/amp-unit dependent)
mpd_ms           = 5;    % MinPeakDistance in ms (per channel)
merge_ms         = 2;    % cross-channel merge window in ms

% DBSCAN clustering of event positions (per file)
epsPhys  = 20;  % cm neighborhood radius
minPts   = 30;  % min events per cluster

% Kalman/association parameters are inside kf_update_tracks_amp()

%% Containers for results across files
trackLog = [];  % rows: [fileIdx, ID, x, y, f]
fileSumm = struct('file',{},'Fs',{},'nEvents',{},'nClusters',{},'centroids',{},'freqs',{},'trackIDs',{});

%% Main loop over files
for i = 1:numel(files)
    fname = files(i).name;
    fpath = fullfile(myDir, fname);
    fprintf('\n--- Processing %d/%d: %s ---\n', i, numel(files), fname);

    %% Read audio
    try
        [y, Fs] = audioread(fpath);    % y: [N x C]
    catch ME
        fprintf('Archivo %d %s defectuoso (read): %s\n', i, fname, ME.message);
        continue;
    end
    if isempty(y)
        fprintf('Archivo %d %s defectuoso (empty).\n', i, fname);
        continue;
    end
    if size(y,2) ~= numChExpected
        fprintf('Archivo %d %s channel mismatch (got %d, expected %d). Skipping.\n',...
            i, fname, size(y,2), numChExpected);
        continue;
    end

    tic;

    %% Filter (bandpass 300–2000 Hz) and compute Hilbert envelope
    Wn = [300 2000] / (Fs/2);
    [b, a] = butter(3, Wn, 'bandpass');
    dat = filter(b, a, y);
    env = abs(hilbert(dat));            % envelope per channel
    env = env - median(env,1);          % per-channel baseline
    env = max(env, 0);                  % nonnegative
    env = env.^2; % help to put signal far from noise, after the EOD detection we use nthroot(max(snap,0), 3 or 4) 

    %% Per-channel peak detection
    mpd = max(1, round((mpd_ms/1000) * Fs));
    EODti = [];
    for c = 1:numChExpected
        [~, locs] = findpeaks(env(:,c), 'MinPeakHeight', minPkHeight, 'MinPeakDistance', mpd);
        EODti = [EODti; locs]; %#ok<AGROW>
    end
    if isempty(EODti)
        fprintf('Archivo %d %s no hay peces (sin picos).\n', i, fname);
        continue;
    end

    %% Merge across channels (dedupe near-simultaneous detections)
    EODtis = sort(EODti);
    mergeSamp = max(1, round((merge_ms/1000) * Fs));
    keep = [true; diff(EODtis) >= mergeSamp];
    EODt = EODtis(keep);

    if numel(EODt) < minEventsForFish
        fprintf('Archivo %d %s no hay peces (pocos eventos: %d).\n', i, fname, numel(EODt));
        continue;
    end

    %% Amplitude snapshots at event times
    snap = env(EODt, :);   % (#events x 8)

    %% Event localization via spatial interpolation
    vals = nthroot(max(snap,0), 3);    % compress dynamic range. This “softens” big amplitudes so a single hot channel doesn.t dominate the fingerprint. It3s like a gentle, log-like compression but defined at 0 (no -Inf).
    %vals = snap;
    rs = sum(vals,2); rs(rs==0) = eps; vals = vals ./ rs;  % per-event normalization
    weighted = vals * Wproto;          % (#events x nGrid)
    [~, kmax] = max(weighted, [], 2);
    X = XYgrid(kmax,1);  % cm
    Y = XYgrid(kmax,2);  % cm

        figure; plot(X,Y,'or') %to check the X and Y are ok
    %% Cluster events into fish (per file)
    try
        labels = dbscan([Y X], epsPhys, minPts);  % note [y x] to match earlier plots
    catch ME
        error('dbscan() required: %s', ME.message);
    end
    good = labels > 0;
    if ~any(good)
        fprintf('Archivo %d %s no hay peces (sin clusters).\n', i, fname);
        continue;
    end

    uniq = unique(labels(good));
    M = numel(uniq);

    %% Per-cluster frequency (Hz), centroid (x,y), and amplitude fingerprint
    centroids = nan(M,2);
    freqs     = nan(M,1);

    S = compute_amp_fingerprints(snap, labels, 20);  % fingerprints per cluster label
    ampSigs = nan(M, numChExpected);

    for m = 1:M
        k = uniq(m);
        ev = find(labels == k);
        % Frequency as median inverse ISI for events in this cluster
        tt = sort(EODt(ev));
        if numel(tt) >= 2
            isis = diff(tt) / Fs;
            freqs(m) = 1 ./ median(isis, 'omitnan');
        end
        centroids(m,:) = [median(X(ev)), median(Y(ev))];
        hit = find([S.k] == k, 1);
        if ~isempty(hit), ampSigs(m,:) = S(hit).sig; end
    end

    % Sort clusters by frequency (descending)
    [freqs, order] = sort(freqs, 'descend', 'MissingPlacement','last');
    centroids = centroids(order,:);
    ampSigs   = ampSigs(order,:);

    %% Update tracker (keeps IDs across minutes)
    tracks = kf_update_tracks_amp(centroids, freqs, ampSigs, 60);  % dt = 60 s per file

    % Log current estimates
    if ~isempty(tracks)
        cur = arrayfun(@(t) [t.id, t.x(1), t.x(2), t.x(5)], tracks, 'uni', 0);
        cur = vertcat(cur{:});
        for r = 1:size(cur,1)
            trackLog = [trackLog; i, cur(r,:)]; %#ok<AGROW>
        end
        Ttbl = array2table(cur, 'VariableNames', {'ID','x_cm','y_cm','f_Hz'});
        disp(Ttbl);
    end

    % Save file summary
    fileSumm(end+1) = struct('file', fname, 'Fs', Fs, 'nEvents', numel(EODt), ...
        'nClusters', M, 'centroids', centroids, 'freqs', freqs, ...
        'trackIDs', ~isempty(tracks)); %#ok<SAGROW>

    fprintf('Archivo %d %s finalizado en %.2f s\n', i, fname, toc);
end

%% Save outputs
try
    save(fullfile(myDir, 'tracking_results.mat'), 'trackLog', 'fileSumm');
    fprintf('\nSaved tracking_results.mat in %s\n', myDir);
catch ME
    warning('tracker:SaveFailed', 'Could not save results: %s', ME.message);
end

%% Optional quick plot: tracks across files
if ~isempty(trackLog)
    figure('Name','Fish tracks (x,y by ID)'); clf; hold on; grid on; axis equal;
    ids = unique(trackLog(:,2));
    for u = reshape(ids,1,[])
        rows = trackLog(:,2)==u;
        plot(trackLog(rows,3), trackLog(rows,4), '-o', 'DisplayName', sprintf('ID %d', u));
    end
    xlim(xlimGrid); ylim(ylimGrid); xlabel('x (cm)'); ylabel('y (cm)'); legend show;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sigs = compute_amp_fingerprints(AUXC, idx, minEventsPerCluster)
% Compute an L2-normalized 8-channel amplitude fingerprint per cluster.
% AUXC: (#events x #channels) amplitudes at consolidated EOD times
% idx : (#events x 1) cluster labels for events (<=0 means noise)
% minEventsPerCluster : minimum events to build a fingerprint

if nargin<3, minEventsPerCluster = 20; end
C = size(AUXC,2);
labels = setdiff(unique(idx(:)), 0);
sigs = struct('k',{},'sig',{},'n',{});
for u = reshape(labels,1,[])
    ev = find(idx==u);
    if numel(ev) < minEventsPerCluster, continue; end
    V = max(AUXC(ev,:), 0);
    V = nthroot(V, 4);                % compress (try 3–5)
    rs = sum(V,2); rs(rs==0) = eps;
    V = V ./ rs;                      % per-event normalization (sum=1)
    fp = median(V,1,'omitnan');
    nrm = norm(fp); if nrm==0, fp(1)=1; nrm=1; end
    fp = fp / nrm;                    % L2-normalize
    sigs(end+1) = struct('k', u, 'sig', fp, 'n', numel(ev)); %#ok<AGROW>
end
end

function tracks = kf_update_tracks_amp(centroidsNEW, jo, ampSigs, dt)
% Kalman tracker with amplitude-aware association.
% Inputs per minute:
%   centroidsNEW: Nx2 [x y] cm
%   jo          : Nx1 frequency (Hz)
%   ampSigs     : NxC L2-normalized fingerprints (NaN row if unavailable)
%   dt          : seconds since previous file (use 60)

if nargin<4, dt = 60; end

% Tuned params for this tank (cm & Hz)
measSig_pos = 6;        % cm RMS measurement noise
measSig_f   = 0.15;     % Hz RMS measurement noise
wander1min  = 40;       % cm RMS position wander per minute
q_acc       = 3*(wander1min^2) / (dt^3);  % cm^2/s^3
freqDrift1m = 0.3;      % Hz RMS per minute
q_f         = (freqDrift1m)^2;            % Hz^2 per step

% Chi-square gate (~99.7% for 3 dims)
try
    gate2 = chi2inv(0.997, 3);
catch
    gate2 = 14.16;
end
hardGate_pos = 60;      % cm hard gate
hardGate_f   = 2.0;     % Hz hard gate

% Amplitude similarity (cosine distance)
sigGate  = 0.35;  % reject if cos-dist > 0.35
sigmaSig = 0.20;  % scale for cost term
w_sig    = 6.0;   % weight of amp term vs Mahalanobis^2
betaSig  = 0.2;   % EMA update rate for track.sig

% Persistent store
persistent T nextID
if isempty(T), T = struct('id',{},'x',{},'P',{},'miss',{},'age',{},'lastDet',{},'sig',{}); nextID = 1; end

% Build detections
N = size(centroidsNEW,1);
dets = struct('z',{},'R',{},'sig',{});
for k = 1:N
    z = [centroidsNEW(k,1); centroidsNEW(k,2); jo(k)];
    R = diag([measSig_pos^2, measSig_pos^2, measSig_f^2]);
    s = ampSigs(k,:);
    if any(isnan(s)), s = []; else, s = s(:).'; end
    dets(end+1) = struct('z', z, 'R', R, 'sig', s); %#ok<AGROW>
end

% Kalman matrices
F = [1 0 dt 0  0;
    0 1 0  dt 0;
    0 0 1  0  0;
    0 0 0  1  0;
    0 0 0  0  1];
H = [1 0 0 0 0;
    0 1 0 0 0;
    0 0 0 0 1];
Qcv = [dt^3/3 dt^2/2; dt^2/2 dt] * q_acc;  % for x and y
Q = blkdiag(Qcv, Qcv, q_f);


% Predict
for i = 1:numel(T)
    T(i).x = F*T(i).x;
    T(i).P = F*T(i).P*F' + Q;

    % ---- Clamp to tank bounds after predict ----
    xmin=0; xmax=80; ymin=0; ymax=120;
    % clip position
    if T(i).x(1) < xmin, T(i).x(1) = xmin; T(i).x(3) = 0; end  % zero vx at wall
    if T(i).x(1) > xmax, T(i).x(1) = xmax; T(i).x(3) = 0; end
    if T(i).x(2) < ymin, T(i).x(2) = ymin; T(i).x(4) = 0; end  % zero vy at wall
    if T(i).x(2) > ymax, T(i).x(2) = ymax; T(i).x(4) = 0; end

    % (optional) extra damping of velocity when the track is missing
    if isfield(T,'miss') && T(i).miss > 0
        damp = 0.9;                 % 0.9–0.98; smaller = stronger damping
        T(i).x(3:4) = damp*T(i).x(3:4);
    end
end

% Association cost (Mahalanobis^2 + amplitude term)
Cmat = inf(numel(T), numel(dets));
for i = 1:numel(T)
    xi = T(i).x; Pi = T(i).P;
    for j = 1:numel(dets)
        z = dets(j).z; R = dets(j).R;
        % Hard gates
        dp = hypot(z(1)-xi(1), z(2)-xi(2));
        df = abs(z(3)-xi(5));
        if dp > hardGate_pos || df > hardGate_f, continue; end
        % Base Mahalanobis^2
        S  = H*Pi*H' + R;
        v  = z - H*xi;
        d2 = v' / S * v;
        if d2 >= gate2, continue; end
        % Amplitude cosine distance penalty
        addCost = 0;
        if isfield(T,'sig') && ~isempty(T(i).sig) && ~isempty(dets(j).sig)
            s1 = T(i).sig(:); s2 = dets(j).sig(:);
            n1 = norm(s1); if n1==0, n1=1; end
            n2 = norm(s2); if n2==0, n2=1; end
            cosSim  = (s1.'*s2)/(n1*n2);
            cosDist = max(0, 1 - cosSim);
            if cosDist > sigGate, continue; end
            addCost = w_sig * (cosDist / sigmaSig)^2;
        end
        Cmat(i,j) = d2 + addCost;
    end
end

% Assign (Hungarian if available; else greedy)
assign = zeros(0,2);
if ~isempty(Cmat) && any(isfinite(Cmat(:)))
    try
        [asgnIdx, ~] = munkres(Cmat);  % if available
        for i = 1:numel(asgnIdx)
            j = asgnIdx(i);
            if j>0 && isfinite(Cmat(i,j)), assign(end+1,:) = [i,j]; end %#ok<AGROW>
        end
    catch
        % Greedy fallback
        usedT = false(1,size(Cmat,1)); usedD = false(1,size(Cmat,2));
        while true
            [minrow, rIdx] = min(Cmat, [], 2);
            [mincost, ii]  = min(minrow);
            if ~isfinite(mincost), break; end
            jj = rIdx(ii);
            if usedT(ii) || usedD(jj), Cmat(ii,jj) = inf; continue; end
            assign(end+1,:) = [ii,jj]; %#ok<AGROW>
            usedT(ii) = true; usedD(jj) = true;
            Cmat(ii,:) = inf; Cmat(:,jj) = inf;
        end
    end
end
assignedTracks = assign(:,1);
assignedDets   = assign(:,2);
unTr = setdiff(1:numel(T), assignedTracks);
unDt = setdiff(1:numel(dets), assignedDets);

% Update matched
for k = 1:size(assign,1)
    i = assign(k,1); j = assign(k,2);
    x = T(i).x; P = T(i).P; z = dets(j).z; R = dets(j).R;

    S = H*P*H' + R; K = P*H'/S;
    x = x + K*(z - H*x);
    P = (eye(5) - K*H)*P;

    % ---- Clamp to tank bounds after update ----
    xmin=0; xmax=80; ymin=0; ymax=120;
    x(1) = min(max(x(1), xmin), xmax);
    x(2) = min(max(x(2), ymin), ymax);

    T(i).x = x; T(i).P = P; T(i).miss = 0; T(i).age = T(i).age + 1; T(i).lastDet = z;
    ...
end

% Missed tracks
for i = unTr
    T(i).miss = T(i).miss + 1; T(i).age = T(i).age + 1;
end

% Initialize new tracks
for j = unDt
    z = dets(j).z;
    x0 = [z(1); z(2); 0; 0; z(3)];
    P0 = diag([ (measSig_pos*2)^2, (measSig_pos*2)^2, ...
        (measSig_pos)^2,   (measSig_pos)^2, ...
        max(measSig_f, freqDrift1m)^2 ]);
    sig0 = [];
    if ~isempty(dets(j).sig), sig0 = dets(j).sig; end
    T(end+1) = struct('id', nextID, 'x', x0, 'P', P0, 'miss', 0, 'age', 1, 'lastDet', z, 'sig', sig0); %#ok<AGROW>
    nextID = nextID + 1;
end

% Prune stale tracks
T = T([T.miss] <= 3);

tracks = T;
end
