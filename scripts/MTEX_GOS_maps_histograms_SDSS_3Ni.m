%% MTEX_GOS_maps_histograms_SDSS_3Ni.m
% Grain orientation spread (GOS) maps and number-fraction histograms
% for LPBF SDSS 2507 + 3%Ni EBSD .ang files.
%
% This script is configured for the folder shown in your MATLAB Online
% screenshot. It uses the current folder as rootDir and processes all .ang
% files present there:
%
%   AS_500x500.ang
%   SR400_500x500.ang
%   SR450_500x500.ang
%   SR500_500x500.ang
%   SR550_500x500.ang
%   SA1100_500x500.ang
%
% Outputs:
%   1) GOS grain maps
%   2) Number-fraction GOS histograms
%   3) Phase-resolved GOS histograms
%   4) Combined all-sample GOS histogram
%   5) Grain-level GOS CSV table
%   6) Summary-statistics CSV table
%
% Notes:
%   - GOS is calculated per reconstructed grain as grains.GOS.
%   - Histogram normalization is number fraction, not area fraction.
%   - Very small grains are removed before the final GOS calculation.

clear; close all; clc;

%% ------------------------------------------------------------------------
%  INITIALIZE MTEX
% -------------------------------------------------------------------------

rootDir = pwd;

% This starts MTEX from the local mtex-6.1.0 folder if MTEX is not already
% available. The validation check is intentionally based on crystalSymmetry,
% not calcGrains, because calcGrains can be method-dispatched in MTEX.
initializeMTEX(rootDir);

%% ------------------------------------------------------------------------
%  USER SETTINGS
% -------------------------------------------------------------------------

% Output directory.
outDir = fullfile(rootDir, 'MTEX_GOS_outputs');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% Grain reconstruction settings.
pars.boundaryAngle_deg      = 5;      % grain boundary threshold in degrees
pars.boundaryAngle          = pars.boundaryAngle_deg * degree;
pars.minGrainPixels         = 5;      % remove grains smaller than this pixel count
pars.smoothIterations       = 1;      % only smooths displayed grain boundaries

% Optional confidence-index filtering.
% Leave empty [] unless you have a defensible CI threshold.
% Example: pars.ciMin = 0.1;
pars.ciMin                  = [];

% Histogram settings.
pars.gosBinWidth_deg        = 0.25;   % GOS bin width in degrees
pars.histMax_deg            = [];     % [] = automatic maximum; or set e.g. 10

% Map color range.
% [] = automatic per map.
% For direct visual comparison among samples, use a fixed range, e.g. [0 5].
pars.mapColorRange_deg      = [];

% Figure/export settings.
% false is safer for batch processing in MATLAB Online.
pars.visibleFigures         = false;
pars.exportResolution       = 600;    % PNG resolution in dpi

% EDAX .ang reference-frame setting.
% For EDAX .ang files, setting 2 is commonly appropriate, but verify map
% orientation against your vendor IPF maps.
pars.edaxReferenceSetting   = 'setting 2';

% Phase definitions.
% Based on the uploaded .ang headers:
%   Phase 1 = Austenite / FCC
%   Phase 2 = Ferrite / BCC
%
% Lattice parameters are approximate and mainly needed for symmetry handling.
CS = { ...
    'notIndexed', ...
    crystalSymmetry('m-3m', [3.59 3.59 3.59], ...
        'mineral', 'Austenite', 'color', [0.85 0.30 0.30]), ...
    crystalSymmetry('m-3m', [2.87 2.87 2.87], ...
        'mineral', 'Ferrite', 'color', [0.30 0.45 0.85]) ...
    };

phaseNames = {'Austenite', 'Ferrite'};
phaseIds   = [1, 2];

%% ------------------------------------------------------------------------
%  FIND .ANG FILES
% -------------------------------------------------------------------------

angFiles = dir(fullfile(rootDir, '*.ang'));

if isempty(angFiles)
    error('No .ang files found in the current folder: %s', rootDir);
end

% Sort files by name for reproducible processing order.
[~, sortIdx] = sort({angFiles.name});
angFiles = angFiles(sortIdx);

fprintf('\nFound %d .ang files in:\n%s\n\n', numel(angFiles), rootDir);

for i = 1:numel(angFiles)
    fprintf('  %d) %s\n', i, angFiles(i).name);
end
fprintf('\n');

%% ------------------------------------------------------------------------
%  STORAGE TABLES
% -------------------------------------------------------------------------

allGrainRows = table();
summaryRows  = table();

combinedGOS = struct( ...
    'sample', {}, ...
    'all', {}, ...
    'Austenite', {}, ...
    'Ferrite', {} ...
    );

%% ------------------------------------------------------------------------
%  MAIN PROCESSING LOOP
% -------------------------------------------------------------------------

for k = 1:numel(angFiles)

    fname = fullfile(angFiles(k).folder, angFiles(k).name);
    [~, sampleNameRaw, ~] = fileparts(fname);

    % Remove common suffix and make a valid MATLAB-safe name.
    sampleName = erase(sampleNameRaw, '_500x500');
    sampleName = matlab.lang.makeValidName(sampleName);

    fprintf('------------------------------------------------------------\n');
    fprintf('[%d/%d] Processing sample: %s\n', k, numel(angFiles), sampleName);
    fprintf('File: %s\n', fname);

    %% Load EBSD data.
    ebsd = loadSDSSang(fname, CS, pars);

    %% Optional CI filtering.
    ebsd = applyCIFilter(ebsd, pars.ciMin, sampleName);

    %% Keep indexed data only.
    ebsd = ebsd('indexed');

    if isempty(ebsd)
        warning('No indexed pixels retained for %s. Skipping.', sampleName);
        continue;
    end

    fprintf('Indexed EBSD points retained: %d\n', length(ebsd));

    %% First grain reconstruction.
    [grains0, ebsd.grainId, ebsd.mis2mean] = calcGrains( ...
        ebsd, ...
        'angle', pars.boundaryAngle ...
        );

    fprintf('Initial reconstructed grains: %d\n', length(grains0));

    %% Remove very small grains.
    if pars.minGrainPixels > 1

        smallGrains = grains0(grains0.grainSize < pars.minGrainPixels);

        fprintf('Small grains removed, grainSize < %d pixels: %d\n', ...
            pars.minGrainPixels, length(smallGrains));

        if ~isempty(smallGrains)
            ebsd(smallGrains) = [];
        end
    end

    %% Final grain reconstruction after cleanup.
    [grains, ebsd.grainId, ebsd.mis2mean] = calcGrains( ...
        ebsd, ...
        'angle', pars.boundaryAngle ...
        );

    if isempty(grains)
        warning('No grains reconstructed for %s after cleanup. Skipping.', sampleName);
        continue;
    end

    fprintf('Final reconstructed grains: %d\n', length(grains));

    %% Calculate GOS in degrees.
    gosAll_deg = grains.GOS ./ degree;

    fprintf('Mean GOS:   %.3f deg\n', mean(gosAll_deg, 'omitnan'));
    fprintf('Median GOS: %.3f deg\n', median(gosAll_deg, 'omitnan'));
    fprintf('Max GOS:    %.3f deg\n', max(gosAll_deg));

    %% Grain-level table.
    grainTable = makeGrainTable(sampleName, grains, gosAll_deg, phaseIds, phaseNames);

    if isempty(allGrainRows)
        allGrainRows = grainTable;
    else
        allGrainRows = [allGrainRows; grainTable]; %#ok<AGROW>
    end

    %% Summary table: all grains.
    summaryTable = summarizeGOS(sampleName, 'All indexed grains', gosAll_deg);

    if isempty(summaryRows)
        summaryRows = summaryTable;
    else
        summaryRows = [summaryRows; summaryTable]; %#ok<AGROW>
    end

    %% Summary table: phase-resolved grains.
    for p = 1:numel(phaseNames)

        gp = getPhaseGrains(grains, phaseNames{p}, phaseIds(p));

        if isempty(gp)
            warning('No %s grains found for %s.', phaseNames{p}, sampleName);
            phaseSummary = summarizeGOS(sampleName, phaseNames{p}, []);
        else
            gosPhase_deg = gp.GOS ./ degree;
            phaseSummary = summarizeGOS(sampleName, phaseNames{p}, gosPhase_deg);
        end

        summaryRows = [summaryRows; phaseSummary]; %#ok<AGROW>
    end

    %% Store for combined histogram.
    combinedGOS(end+1).sample = sampleName; %#ok<SAGROW>
    combinedGOS(end).all = gosAll_deg;

    for p = 1:numel(phaseNames)

        gp = getPhaseGrains(grains, phaseNames{p}, phaseIds(p));

        if isempty(gp)
            combinedGOS(end).(phaseNames{p}) = [];
        else
            combinedGOS(end).(phaseNames{p}) = gp.GOS ./ degree;
        end
    end

    %% Export figures.
    exportGOSMap(grains, gosAll_deg, sampleName, outDir, pars);

    exportSampleHistograms( ...
        grains, ...
        sampleName, ...
        outDir, ...
        pars, ...
        phaseNames, ...
        phaseIds ...
        );

    %% Save MATLAB workspace for this sample.
    save(fullfile(outDir, [sampleName '_GOS_workspace.mat']), ...
        'fname', ...
        'sampleName', ...
        'ebsd', ...
        'grains', ...
        'gosAll_deg', ...
        'pars', ...
        '-v7.3' ...
        );

    fprintf('Finished sample: %s\n\n', sampleName);
end

%% ------------------------------------------------------------------------
%  EXPORT TABLES
% -------------------------------------------------------------------------

if ~isempty(allGrainRows)

    grainCsv = fullfile(outDir, 'grain_level_GOS_values.csv');
    writetable(allGrainRows, grainCsv);

    fprintf('Saved grain-level GOS table:\n%s\n', grainCsv);
else
    warning('No grain-level GOS table was generated.');
end

if ~isempty(summaryRows)

    summaryCsv = fullfile(outDir, 'GOS_summary_statistics.csv');
    writetable(summaryRows, summaryCsv);

    fprintf('Saved GOS summary table:\n%s\n', summaryCsv);
else
    warning('No GOS summary table was generated.');
end

%% ------------------------------------------------------------------------
%  COMBINED HISTOGRAMS
% -------------------------------------------------------------------------

if ~isempty(combinedGOS)
    exportCombinedHistogram(combinedGOS, outDir, pars);
end

fprintf('\n============================================================\n');
fprintf('GOS analysis complete.\n');
fprintf('Outputs saved in:\n%s\n', outDir);
fprintf('============================================================\n');

%% ========================================================================
%  LOCAL FUNCTIONS
% ========================================================================

function initializeMTEX(rootDir)

    % Check whether MTEX is already available.
    if exist('crystalSymmetry', 'file') == 2 || exist('crystalSymmetry', 'class') == 8
        fprintf('MTEX appears to be already available on the MATLAB path.\n');
        return;
    end

    % Try local MTEX folder.
    localMtexDir = fullfile(rootDir, 'mtex-6.1.0');
    localStartup = fullfile(localMtexDir, 'startup_mtex.m');

    if exist(localStartup, 'file') == 2

        fprintf('Starting MTEX from:\n%s\n', localMtexDir);

        currentFolder = pwd;
        cleanupObj = onCleanup(@() cd(currentFolder));

        cd(localMtexDir);
        startup_mtex;

    elseif exist('startup_mtex', 'file') == 2

        fprintf('Starting MTEX using startup_mtex found on MATLAB path.\n');
        startup_mtex;

    else

        error([ ...
            'MTEX was not found. ', ...
            'Run startup_mtex manually or place the mtex-6.1.0 folder in the current directory.' ...
            ]);
    end

    % Safer validation:
    % crystalSymmetry is sufficient to confirm that the MTEX path has been
    % initialized. Do not test calcGrains here because in MTEX it may be
    % resolved as a class method rather than a plain file.
    if ~(exist('crystalSymmetry', 'file') == 2 || exist('crystalSymmetry', 'class') == 8)
        error('MTEX startup was attempted, but crystalSymmetry is still unavailable.');
    end

    fprintf('MTEX startup check passed.\n');

end

function ebsd = loadSDSSang(fname, CS, pars)

    % Load EDAX .ang data with explicit phase definitions.
    %
    % Primary route:
    %   EDAX reference-frame setting 2.
    %
    % Fallbacks:
    %   1) explicit CS with convertEuler2SpatialReferenceFrame
    %   2) header/autodetect import

    try

        ebsd = EBSD.load( ...
            fname, ...
            CS, ...
            'interface', 'ang', ...
            'convertEuler2SpatialReferenceFrame', pars.edaxReferenceSetting ...
            );

    catch ME1

        warning(['Primary EDAX import failed for:\n%s\n', ...
                 'Reason:\n%s\n', ...
                 'Trying fallback import with convertEuler2SpatialReferenceFrame only.'], ...
                 fname, ME1.message);

        try

            ebsd = EBSD.load( ...
                fname, ...
                CS, ...
                'interface', 'ang', ...
                'convertEuler2SpatialReferenceFrame' ...
                );

        catch ME2

            warning(['Explicit-CS fallback import failed for:\n%s\n', ...
                     'Reason:\n%s\n', ...
                     'Trying header/autodetect import.'], ...
                     fname, ME2.message);

            ebsd = EBSD.load( ...
                fname, ...
                'interface', 'ang', ...
                'convertEuler2SpatialReferenceFrame' ...
                );
        end
    end

end

function ebsd = applyCIFilter(ebsd, ciMin, sampleName)

    if isempty(ciMin)
        return;
    end

    ci = [];

    if isfield(ebsd.prop, 'ci')
        ci = ebsd.prop.ci;
    elseif isfield(ebsd.prop, 'CI')
        ci = ebsd.prop.CI;
    elseif isprop(ebsd, 'ci')
        ci = ebsd.ci;
    end

    if isempty(ci)
        warning('CI filter requested for %s, but no CI field was found. No CI filter applied.', sampleName);
        return;
    end

    nBefore = length(ebsd);
    ebsd(ci < ciMin) = [];
    nAfter = length(ebsd);

    fprintf('CI filter applied to %s: CI >= %.3f\n', sampleName, ciMin);
    fprintf('Pixels before CI filter: %d\n', nBefore);
    fprintf('Pixels after CI filter:  %d\n', nAfter);

end

function gp = getPhaseGrains(grains, phaseName, phaseId)

    % Prefer mineral-name selection. If this fails, fall back to numeric
    % phase selection.

    gp = [];

    try
        gp = grains(phaseName);
        return;
    catch
    end

    try
        gp = grains(grains.phase == phaseId);
    catch
        gp = [];
    end

end

function T = makeGrainTable(sampleName, grains, gos_deg, phaseIds, phaseNames)

    n = length(grains);

    grainId = nan(n, 1);
    phaseId = nan(n, 1);
    phaseName = strings(n, 1);
    grainSizePixels = nan(n, 1);
    area_um2 = nan(n, 1);
    equivalentDiameter_um = nan(n, 1);

    % Grain ID.
    try
        grainId = double(grains.id(:));
    catch
        grainId = (1:n).';
    end

    % Phase ID and phase name.
    try
        phaseId = double(grains.phase(:));

        for i = 1:n
            idx = find(phaseIds == phaseId(i), 1, 'first');

            if isempty(idx)
                phaseName(i) = "Unknown";
            else
                phaseName(i) = string(phaseNames{idx});
            end
        end

    catch
        phaseName(:) = "Unknown";
    end

    % Grain size in pixels.
    try
        grainSizePixels = double(grains.grainSize(:));
    catch
    end

    % Grain area.
    try
        area_um2 = double(grains.area(:));
    catch
    end

    % Equivalent circular diameter.
    try
        equivalentDiameter_um = 2 .* sqrt(area_um2 ./ pi);
    catch
    end

    T = table( ...
        repmat(string(sampleName), n, 1), ...
        grainId, ...
        phaseId, ...
        phaseName, ...
        grainSizePixels, ...
        area_um2, ...
        equivalentDiameter_um, ...
        double(gos_deg(:)), ...
        'VariableNames', { ...
            'sample', ...
            'grain_id', ...
            'phase_id', ...
            'phase_name', ...
            'grain_size_pixels', ...
            'area_um2', ...
            'equivalent_diameter_um', ...
            'GOS_deg' ...
        } ...
        );

end

function T = summarizeGOS(sampleName, phaseName, gos_deg)

    gos_deg = double(gos_deg(:));
    gos_deg = gos_deg(isfinite(gos_deg));

    if isempty(gos_deg)

        T = table( ...
            string(sampleName), ...
            string(phaseName), ...
            0, ...
            NaN, ...
            NaN, ...
            NaN, ...
            NaN, ...
            NaN, ...
            NaN, ...
            NaN, ...
            'VariableNames', { ...
                'sample', ...
                'phase', ...
                'n_grains', ...
                'mean_GOS_deg', ...
                'std_GOS_deg', ...
                'median_GOS_deg', ...
                'p75_GOS_deg', ...
                'p90_GOS_deg', ...
                'p95_GOS_deg', ...
                'max_GOS_deg' ...
            } ...
            );

        return;
    end

    T = table( ...
        string(sampleName), ...
        string(phaseName), ...
        numel(gos_deg), ...
        mean(gos_deg, 'omitnan'), ...
        std(gos_deg, 'omitnan'), ...
        median(gos_deg, 'omitnan'), ...
        localPercentile(gos_deg, 75), ...
        localPercentile(gos_deg, 90), ...
        localPercentile(gos_deg, 95), ...
        max(gos_deg), ...
        'VariableNames', { ...
            'sample', ...
            'phase', ...
            'n_grains', ...
            'mean_GOS_deg', ...
            'std_GOS_deg', ...
            'median_GOS_deg', ...
            'p75_GOS_deg', ...
            'p90_GOS_deg', ...
            'p95_GOS_deg', ...
            'max_GOS_deg' ...
        } ...
        );

end

function q = localPercentile(x, p)

    x = double(x(:));
    x = x(isfinite(x));
    x = sort(x);

    if isempty(x)
        q = NaN;
        return;
    end

    if numel(x) == 1
        q = x;
        return;
    end

    pos = 1 + (numel(x) - 1) * p / 100;
    lo = floor(pos);
    hi = ceil(pos);

    if lo == hi
        q = x(lo);
    else
        q = x(lo) + (pos - lo) * (x(hi) - x(lo));
    end

end

function edges = makeEdges(gos_deg, pars)

    gos_deg = double(gos_deg(:));
    gos_deg = gos_deg(isfinite(gos_deg));

    if isempty(gos_deg)
        edges = 0:pars.gosBinWidth_deg:1;
        return;
    end

    if isempty(pars.histMax_deg)
        xmax = ceil(max(gos_deg) / pars.gosBinWidth_deg) * pars.gosBinWidth_deg;
        xmax = max(xmax, pars.gosBinWidth_deg);
    else
        xmax = pars.histMax_deg;
    end

    edges = 0:pars.gosBinWidth_deg:xmax;

    if numel(edges) < 2
        edges = [0, pars.gosBinWidth_deg];
    end

end

function exportGOSMap(grains, gos_deg, sampleName, outDir, pars)

    fig = makeFigure(pars, [100 100 950 780]);

    grainsToPlot = grains;

    if pars.smoothIterations > 0
        try
            grainsToPlot = smooth(grains, pars.smoothIterations);
        catch
            grainsToPlot = grains;
        end
    end

    plot(grainsToPlot, gos_deg, 'micronbar', 'on');
    hold on;

    try
        plot(grainsToPlot.boundary, 'lineWidth', 0.5);
    catch
    end

    hold off;

    mtexColorbar('title', 'GOS (degree)');

    title(sprintf('%s - GOS map', strrep(sampleName, '_', '\_')));

    if ~isempty(pars.mapColorRange_deg)
        try
            setColorRange(pars.mapColorRange_deg);
        catch
            try
                caxis(pars.mapColorRange_deg);
            catch
            end
        end
    end

    outFile = fullfile(outDir, [sampleName '_GOS_map.png']);
    saveFigure(fig, outFile, pars.exportResolution);

    fprintf('Saved GOS map: %s\n', outFile);

    if ~pars.visibleFigures
        close(fig);
    end

end

function exportSampleHistograms(grains, sampleName, outDir, pars, phaseNames, phaseIds)

    gosAll = grains.GOS ./ degree;
    edges = makeEdges(gosAll, pars);

    fig = makeFigure(pars, [100 100 900 950]);

    tiledlayout(3, 1, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact' ...
        );

    % All indexed grains.
    nexttile;

    histogram( ...
        gosAll, ...
        'BinEdges', edges, ...
        'Normalization', 'probability' ...
        );

    xlabel('GOS (degree)');
    ylabel('Number fraction');
    title(sprintf('%s - all indexed grains', strrep(sampleName, '_', '\_')));
    grid on;
    box on;

    % Phase-resolved histograms.
    for p = 1:numel(phaseNames)

        gp = getPhaseGrains(grains, phaseNames{p}, phaseIds(p));

        nexttile;

        if isempty(gp)

            text( ...
                0.5, ...
                0.5, ...
                sprintf('No %s grains found', phaseNames{p}), ...
                'HorizontalAlignment', 'center', ...
                'Units', 'normalized' ...
                );

            axis off;

        else

            histogram( ...
                gp.GOS ./ degree, ...
                'BinEdges', edges, ...
                'Normalization', 'probability' ...
                );

            xlabel('GOS (degree)');
            ylabel('Number fraction');
            title(sprintf('%s - %s', strrep(sampleName, '_', '\_'), phaseNames{p}));
            grid on;
            box on;

        end
    end

    outFile = fullfile(outDir, [sampleName '_GOS_number_fraction_histograms.png']);
    saveFigure(fig, outFile, pars.exportResolution);

    fprintf('Saved GOS histogram: %s\n', outFile);

    if ~pars.visibleFigures
        close(fig);
    end

end

function exportCombinedHistogram(combinedGOS, outDir, pars)

    allVals = [];

    for i = 1:numel(combinedGOS)
        allVals = [allVals; combinedGOS(i).all(:)]; %#ok<AGROW>
    end

    edges = makeEdges(allVals, pars);

    fig = makeFigure(pars, [100 100 1000 720]);

    hold on;

    for i = 1:numel(combinedGOS)

        histogram( ...
            combinedGOS(i).all, ...
            'BinEdges', edges, ...
            'Normalization', 'probability', ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 1.5, ...
            'DisplayName', combinedGOS(i).sample ...
            );

    end

    hold off;

    xlabel('GOS (degree)');
    ylabel('Number fraction');
    title('GOS number-fraction histograms - all indexed grains');
    legend('Location', 'best', 'Interpreter', 'none');
    grid on;
    box on;

    outFile = fullfile(outDir, 'combined_GOS_number_fraction_histogram_all_samples.png');
    saveFigure(fig, outFile, pars.exportResolution);

    fprintf('Saved combined GOS histogram: %s\n', outFile);

    if ~pars.visibleFigures
        close(fig);
    end

end

function fig = makeFigure(pars, position)

    if pars.visibleFigures
        fig = figure('Color', 'w', 'Position', position);
    else
        fig = figure('Color', 'w', 'Position', position, 'Visible', 'off');
    end

end

function saveFigure(fig, outFile, dpi)

    try
        exportgraphics(fig, outFile, 'Resolution', dpi);
    catch
        [folderName, baseName, ~] = fileparts(outFile);
        print(fig, fullfile(folderName, baseName), '-dpng', sprintf('-r%d', dpi));
    end

end
