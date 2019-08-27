function demo_4_psychometric_function_interleave_staircases
% Alex S Baldwin, McGill Vision Research, August 2019
% Demo of fitting a psychometric function to data obtained from the object-
% oriented psychophysical staircase. Requires staircase.m class definition 
% file to be in same folder or on Matlab path. Also requires the Palamedes
% toolbox for function fitting, available from www.palamedestoolbox.org.
% This minimal example shows use of staircase in a simulated NAFC task. The
% intensity units are given in dB, as might be used for contrast detection.
% N.B. it is likely that fits to individual simulation runs will be
% unstable (Palamedes will warn of bootstrap fits not converging), a method
% to avoid this problem will be presented in the next demo
% From: https://github.com/alexsbaldwin/MatlabStaircase

close all

if ~exist('PAL_version','file')
    warning('Palamedes toolbox is either not installed or not on path. Please install Palamedes from www.palamedestoolbox.org')
    error('Failed to locate PAL_version function')
end

PF = @PAL_Logistic;          % using logistic psychometric function
numBootstraps         = 500; % number of bootstrap samples (e.g. 500)
isBootstrapParametric = 1;   % can have parametric or non-parametric bootstrap

% setting the simulated noise level of the model subject for this demo
SIM.simNoiseStdDev  = sqrt(2); % performance-limiting noise in simulation
SIM.nSimulations    = 4;       % number of simulated experiment runs
SIM.trialPauseSec   = 0.01;    % quick pause to allow graphs to be plotted
SIM.numAlternatives = 2;       % e.g. 2 for a 2-alternative forced-choice task
SIM.guessRate       = 1/SIM.numAlternatives;

% Simulate running four repetitions of an experiment condition
for iSim = SIM.nSimulations:-1:1
    sc = run_experiment(SIM);
    allSC(iSim) = sc; % save repetition data for analysis below
end

% In a real experiment, the data should be saved and then loaded from file
% to perform the analysis given below in a separate script. The two functions
% are combined in a single file here for demonstration purposes.

% Perform separate analysis on the data from each simulated run
figure(1)
figpos = [200 200 600 600];
set(gcf, 'Units', 'pixels','PaperUnits', 'points', 'Position', ...
    figpos, 'PaperPosition', figpos, 'Color', [1 1 1]);
hold on

for iSim = 1:SIM.nSimulations
    
    subplot(2,2,iSim)
    hold on
    xlabel('log(stimulus level)')
    ylabel('P(correct)')
    
    sc     = allSC(iSim); % get staircase object for run that is being plotted
    logLev = sc.levels;   % log-transformed stimulus levels
    nT     = sc.nTrials;  % number of trials per stimulus level
    nC     = sc.nCorrect; % number of correct responses per level

    axis([min(logLev),max(logLev),-0.05,1.05])

    % horizontal line at the guess-rate
    plot([min(logLev),max(logLev)],[SIM.guessRate,SIM.guessRate], ...
         'color', [0.5,0.5,0.5], 'linewidth', 2, 'linestyle', '--')

    pCorrect = nC ./ nT;
    % calculate 95% confidence intervals for each stimulus level
    [~,pci]  = binofit(nC, nT);
    eValUpp  = -(pCorrect' - pci(:,2));
    eValLow  = pCorrect' - pci(:,1);

    for i = 1:length(logLev) % for each stimulus level
        if nT(i) > 0
            % scale marker size by number of trials at that stimulus level
            markerSize = 4+ceil(4*nT(i)/max(nT));
            errorbar(logLev(i), pCorrect(i), eValLow(i), eValUpp(i),      ...
                     'marker', 'o', 'color', [0,0,0], 'linestyle', 'none',   ...
                     'markeredgecolor', [0,0,0], 'markerfacecolor', [0,0,0], ...
                     'markersize', markerSize);
        end
    end    

    % searchGrid gives parameter space for initial brute-force search
    % threshold/location parameter: search within range of tested values
    searchGrid.alpha  = linspace(min(logLev),max(logLev),100);
    % psychometric slope parameter: search within a reasonable range
    searchGrid.beta   = logspace(0.1,100,100);
    % guess-rate: fixed at the reciprocal of the number of alternatives
    searchGrid.gamma  = SIM.guessRate;
    % lapse-rate: fixed at 0.5%
    searchGrid.lambda = 0.005;
    % fitting threshold and slope, fixing guess-rate and lapse-rate
    paramsFree = [1, 1, 0, 0];

    % fit psychometric function to data
    paramsValues = PAL_PFML_Fit(logLev, nC, nT, searchGrid, paramsFree, PF);

    % generate and plot model-prediction curve from fitted function
    pfX = linspace(min(logLev),max(logLev),1000);
    pfY = PF(paramsValues,pfX);
    plot(pfX, pfY, 'color', [0.0, 0.5, 1.0], 'linewidth', 2)
    
    % bootstrap psychometric function to data
    if isBootstrapParametric
        [~, paramsSim, ~, ~] = PAL_PFML_BootstrapParametric(logLev, nT, paramsValues, ...
                                                            paramsFree, numBootstraps, PF);
    else
        [~, paramsSim, ~, ~] = PAL_PFML_BootstrapNonParametric(logLev, nC, ...
                                                nT, [], paramsFree, numBootstraps, PF, ...
                                                'searchGrid', searchGrid);
    end

    text(min(logLev)+1,1.0,sprintf('Simulation: %0.0f', iSim))
    
    % estimated threshold from staircase reversals with 95% confidence
    % interval approximated assuming normally-distributed error
    revThresh      = sc.curReversalThresh;
    revThreshUppCI = sc.curReversalThresh + 1.96*sc.curReversalError;
    revThreshLowCI = sc.curReversalThresh - 1.96*sc.curReversalError;
    textStr1 = sprintf('Reversal Threshold = %0.1f, 95%% CI [%0.1f to %0.1f]', ...
                      revThresh, revThreshLowCI, revThreshUppCI);
    text(min(logLev)+1,0.15,textStr1)
    
    % threshold estimate from psychometric function fit with 95% confidence
    % intervals obtained from bootstrapping (parametric or non-parametric)
    fitThresh      = paramsValues(1);
    fitThreshUppCI = prctile(paramsSim(:,1), 97.5);
    fitThreshLowCI = prctile(paramsSim(:,1), 02.5);
    textStr2 = sprintf('Logistic Threshold = %0.1f, 95%% CI [%0.1f to %0.1f]', ...
                        fitThresh, fitThreshLowCI, fitThreshUppCI);
	text(min(logLev)+1,0.05,textStr2)
    
    pause(SIM.trialPauseSec);

end

% Perform analysis on data combined over every run
figure(2)
figpos = [200 200 400 400];
set(gcf, 'Units', 'pixels','PaperUnits', 'points', 'Position', ...
    figpos, 'PaperPosition', figpos, 'Color', [1 1 1]);
hold on
xlabel('log(stimulus level)')
ylabel('P(correct)')

% combine data across all of the simulations by summing
for iSim = 1:SIM.nSimulations
    if iSim == 1
        nT = allSC(iSim).nTrials;
        nC = allSC(iSim).nCorrect;
    else
        nT = nT + allSC(iSim).nTrials;
        nC = nC + allSC(iSim).nCorrect;
    end
end

% get average reversal threshold across simulations with standard error
allRevThresh    = [allSC(:).curReversalThresh];
revThreshMean   = mean(allRevThresh);
revThreshStdErr = std(allRevThresh)./sqrt(length(allRevThresh));

axis([min(logLev),max(logLev),-0.05,1.05])

% horizontal line at the guess-rate
plot([min(logLev),max(logLev)],[SIM.guessRate,SIM.guessRate], ...
     'color', [0.5,0.5,0.5], 'linewidth', 2, 'linestyle', '--')

pCorrect = nC ./ nT;
[~,pci]  = binofit(nC, nT);
eValUpp  = -(pCorrect' - pci(:,2));
eValLow  = pCorrect' - pci(:,1);

for i = 1:length(logLev) % for each stimulus level
    if nT(i) > 0
        % scale marker size by number of trials at that stimulus level
        markerSize = 4+ceil(4*nT(i)/max(nT));
        errorbar(logLev(i), pCorrect(i), eValLow(i), eValUpp(i),      ...
                 'marker', 'o', 'color', [0,0,0], 'linestyle', 'none',   ...
                 'markeredgecolor', [0,0,0], 'markerfacecolor', [0,0,0], ...
                 'markersize', markerSize);
    end
end

% searchGrid gives parameter space for initial brute-force search
% threshold/location parameter: search within range of tested values
searchGrid.alpha  = linspace(min(logLev),max(logLev),100);
% psychometric slope parameter: search within a reasonable range
searchGrid.beta   = logspace(0.1,100,100);
% guess-rate: fixed at the reciprocal of the number of alternatives
searchGrid.gamma  = SIM.guessRate;
% lapse-rate: fixed at 0.5%
searchGrid.lambda = 0.005;
% fitting threshold and slope, fixing guess-rate and lapse-rate
paramsFree = [1, 1, 0, 0];

% fit psychometric function to data
paramsValues = PAL_PFML_Fit(logLev, nC, nT, searchGrid, paramsFree, PF);

% generate and plot model-prediction curve from fitted function
pfX = linspace(min(logLev),max(logLev),1000);
pfY = PF(paramsValues,pfX);
plot(pfX, pfY, 'color', [0.0, 0.5, 1.0], 'linewidth', 2)

% bootstrap psychometric function to data
if isBootstrapParametric
    [~, paramsSim, ~, ~] = PAL_PFML_BootstrapParametric(logLev, nT, paramsValues, ...
                                                        paramsFree, numBootstraps, PF);
else
    [~, paramsSim, ~, ~] = PAL_PFML_BootstrapNonParametric(logLev, nC, ...
                                            nT, [], paramsFree, numBootstraps, PF, ...
                                            'searchGrid', searchGrid);
end
  
% estimated threshold from staircase reversals with 95% confidence
% interval approximated assuming normally-distributed error
revThresh      = revThreshMean;
revThreshUppCI = revThreshMean + 1.96*revThreshStdErr;
revThreshLowCI = revThreshMean - 1.96*revThreshStdErr;
textStr1 = sprintf('Reversal Threshold = %0.1f, 95%% CI [%0.1f to %0.1f]', ...
                   revThresh, revThreshLowCI, revThreshUppCI);
text(min(logLev)+1,0.15,textStr1)

% threshold estimate from psychometric function fit with 95% confidence
% intervals obtained from bootstrapping (parametric or non-parametric)
fitThresh      = paramsValues(1);
fitThreshUppCI = prctile(paramsSim(:,1), 97.5);
fitThreshLowCI = prctile(paramsSim(:,1), 02.5);
textStr2 = sprintf('Logistic Threshold = %0.1f, 95%% CI [%0.1f to %0.1f]', ...
                    fitThresh, fitThreshLowCI, fitThreshUppCI);
text(min(logLev)+1,0.05,textStr2)

return

function sc = run_experiment(SIM)

    logStimLevels   = -24:3:24; % stimulus levels (here in dB logarithmic units)
    initStepSize    = 6;        % staircase step size before first reversal
    stepSize        = 3;        % staircase step size after first reversal
    nRightToDescend = 3;        % number of "correct" responses before descending
    nWrongToAscend  = 1;        % number of "incorrect" responses before ascending
    maxNumTrials    = 100;      % maximum number of trials for staircase
    maxNumReversals = 12;       % maximum number of reversals for staircase
    startLevel      = 12;       % starting level of staircase (here in dB units)
    verbose         = 1;        % set to 1 for verbose staircase output

    % class constructor function returns staircase object
    sc = staircase(logStimLevels, initStepSize, stepSize,         ...
                   nRightToDescend, nWrongToAscend, maxNumTrials, ...
                   maxNumReversals, startLevel, verbose);

    while ~sc.isFinished           % keep looping until staircase concludes

        logStimLevel = sc.curLevel;          % testing level from staircase
        linStimLevel = 10^(logStimLevel/20); % convert to linear units

        % in any real study the stimulus presentation and response collection
        % code would go here, for example purposes we instead take responses
        % from a simulated subject with response behaviour as defined by SIM
        isSimCorrect = do_sim(SIM, linStimLevel); % 0 = incorrect, 1 = correct

        % sending the Boolean (0 or 1) isCorrect to the staircase
        sc.doResp(isSimCorrect); % calls "doResp" method of staircase

        % current rough estimate of threshold from staircase reversals
        % n.b. currently uses ALL reversals after the first reversal
        fprintf('Current threshold from reversals: %0.2f +/- %0.2f\n\n', ...
                sc.curReversalThresh, sc.curReversalError)

    end

return

function isSimCorrect = do_sim(SIM, linStimLevel)
    % Alex S Baldwin, McGill Vision Research, July 2019
    % Simulate subject behaviour for NAFC task

    % simulated "noisy internal responses" from the n intervals
    respT = randn * SIM.simNoiseStdDev + linStimLevel;      % target interval
    numNullChannels = SIM.numAlternatives - 1;
    respN = randn(numNullChannels,1) .* SIM.simNoiseStdDev; % null interval
    
    % simulated subject responds correctly if noisy response to target is 
    % greater than all noisy responses to null intervals
    isSimCorrect = respT > max(respN); 

return