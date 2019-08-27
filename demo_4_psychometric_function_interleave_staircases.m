function demo_4_psychometric_function_interleave_staircases
% Alex S Baldwin, McGill Vision Research, August 2019
% Demo of fitting a psychometric function to data obtained from the object-
% oriented psychophysical staircase. In this example two staircases target
% different percent-correct points on the psychometric function. This does
% a better job of constraining the psychometric function fit. Requires 
% staircase.m class definition file to be in same folder or on Matlab path.
% Also requires the Palamedes toolbox for function fitting, available from 
% www.palamedestoolbox.org. This minimal example shows use of staircase in
% a simulated NAFC task. The intensity units are given in dB, as might be 
% used for contrast detection.
% From: https://github.com/alexsbaldwin/MatlabStaircase

close all

if ~exist('PAL_version','file')
    warning('Palamedes toolbox is either not installed or not on path. Please install Palamedes from www.palamedestoolbox.org')
    error('Failed to locate PAL_version function')
end

PF = @PAL_Logistic;           % using logistic psychometric function
numBootstraps         = 500; % number of bootstrap samples (e.g. 500)
isBootstrapParametric = 1;    % can have parametric or non-parametric bootstrap

% setting the simulated noise level of the model subject for this demo
SIM.simNoiseStdDev  = sqrt(2); % performance-limiting noise in simulation
SIM.nSimulations    = 4;       % number of simulated experiment runs
SIM.trialPauseSec   = 0.01;    % quick pause to allow graphs to be plotted
SIM.numAlternatives = 2;       % e.g. 2 for a 2-alternative forced-choice task
SIM.guessRate       = 1/SIM.numAlternatives;
% Array of staircase "rules" for the number of correct responses needed to
% make the level of the staircase decrease to a more difficult value.
% Different numbers for this rule cause the staircase to target different
% points on the psychometric function. Interleaving different rules results
% in sampling more broadly across the psychometric function, giving a
% better constraint on the psychometric function fit. e.g. [2, 3] gives
% two staircases, one is a two-down one-up, the other a three-down one-up.
SIM.scRightRules    = [2, 3]; 
SIM.numStaircases   = length(SIM.scRightRules);

% Simulate running an experiment condition
staircaseArray = run_experiment(SIM);

% In a real experiment, the data should be saved and then loaded from file
% to perform the analysis given below in a separate script. The two functions
% are combined in a single file here for demonstration purposes.

% Perform analysis on data combined over staircases
figure(2)
figpos = [200 200 400 400];
set(gcf, 'Units', 'pixels','PaperUnits', 'points', 'Position', ...
    figpos, 'PaperPosition', figpos, 'Color', [1 1 1]);
hold on
xlabel('log(stimulus level)')
ylabel('P(correct)')

scRevThreshold      = zeros(SIM.numStaircases,1);
scRevThresholdError = zeros(SIM.numStaircases,1);
% combine data across all of the staircases by summing
for iSC = 1:SIM.numStaircases
    if iSC == 1
        logLev = staircaseArray(iSC).levels;
        nT     = staircaseArray(iSC).nTrials;
        nC     = staircaseArray(iSC).nCorrect;
    else
        nT = nT + staircaseArray(iSC).nTrials;
        nC = nC + staircaseArray(iSC).nCorrect;
    end
    scRevThreshold(iSC)      = staircaseArray(iSC).curReversalThresh;
    scRevThresholdError(iSC) = staircaseArray(iSC).curReversalError;
end

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
for iSC = 1:SIM.numStaircases
    revThresh      = staircaseArray(iSC).curReversalThresh;
    revThreshUppCI = revThresh + 1.96*staircaseArray(iSC).curReversalError;
    revThreshLowCI = revThresh - 1.96*staircaseArray(iSC).curReversalError;
    textStr1 = sprintf('Staircase %0.0f: Reversal Threshold = %0.1f, 95%% CI [%0.1f to %0.1f]', ...
                       iSC, revThresh, revThreshLowCI, revThreshUppCI);
    text(min(logLev)+1,0.05+0.1*iSC,textStr1)
end

% threshold estimate from psychometric function fit with 95% confidence
% intervals obtained from bootstrapping (parametric or non-parametric)
fitThresh      = paramsValues(1);
fitThreshUppCI = prctile(paramsSim(:,1), 97.5);
fitThreshLowCI = prctile(paramsSim(:,1), 02.5);
textStr2 = sprintf('Logistic Threshold = %0.1f, 95%% CI [%0.1f to %0.1f]', ...
                    fitThresh, fitThreshLowCI, fitThreshUppCI);
text(min(logLev)+1,0.05,textStr2)

return

function staircaseArray = run_experiment(SIM)

    logStimLevels   = -24:3:24; % stimulus levels (here in dB logarithmic units)
    initStepSize    = 6;        % staircase step size before first reversal
    stepSize        = 3;        % staircase step size after first reversal
    nWrongToAscend  = 1;        % number of "incorrect" responses before ascending
    maxNumTrials    = 120;      % maximum number of trials for staircase
    maxNumReversals = 14;       % maximum number of reversals for staircase
    startLevel      = 15;       % starting level of staircase (here in dB units)
    verbose         = 1;        % set to 1 for verbose staircase output

    for iSC = SIM.numStaircases:-1:1
        % number of "correct" responses before descending
        nRightToDescend = SIM.scRightRules(iSC);
        % class constructor function returns staircase object
        staircaseArray(iSC) = staircase(logStimLevels, initStepSize, ...
                                        stepSize, nRightToDescend,   ...
                                        nWrongToAscend, maxNumTrials, ...
                                        maxNumReversals, startLevel, verbose);
    end
    
    % keep looping until all staircases conclude
    while ~all([staircaseArray(:).isFinished])
        
        % select a random unfinished staircase
        iSC = randi(SIM.numStaircases);
        while staircaseArray(iSC).isFinished
            iSC = randi(SIM.numStaircases);
        end
        
        % testing level from current staircase
        logStimLevel = staircaseArray(iSC).curLevel;
        linStimLevel = 10^(logStimLevel/20); % convert to linear units

        % in any real study the stimulus presentation and response collection
        % code would go here, for example purposes we instead take responses
        % from a simulated subject with response behaviour as defined by SIM
        isSimCorrect = do_sim(SIM, linStimLevel); % 0 = incorrect, 1 = correct

        % sending the Boolean (0 or 1) isCorrect to the staircase
        staircaseArray(iSC).doResp(isSimCorrect); % calls "doResp" method of staircase

        % current rough estimate of threshold from staircase reversals
        % n.b. currently uses ALL reversals after the first reversal
        fprintf('Staircase %0.0f: current threshold from reversals: %0.2f +/- %0.2f\n\n', ...
                iSC, staircaseArray(iSC).curReversalThresh, staircaseArray(iSC).curReversalError)

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