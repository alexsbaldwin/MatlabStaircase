function demo_1_minimal_example
% Alex S Baldwin, McGill Vision Research, July 2019
% Demo of object-oriented psychophysical staircase implementation. Requires
% staircase.m class definition file to be in same folder or on Matlab path.
% This minimal example shows use of staircase in a simulated 2AFC task. The
% intensity units are given in dB, as might be used for contrast detection.
% From: https://github.com/alexsbaldwin/MatlabStaircase

logStimLevels   = -24:3:24; % stimulus levels (here in dB logarithmic units)
initStepSize    = 6;        % staircase step size before first reversal
stepSize        = 3;        % staircase step size after first reversal
nRightToDescend = 3;        % number of "correct" responses before descending
nWrongToAscend  = 1;        % number of "incorrect" responses before ascending
maxNumTrials    = 100;      % maximum number of trials for staircase
maxNumReversals = 10;       % maximum number of reversals for staircase
startLevel      = 12;       % starting level of staircase (here in dB units)
verbose         = 1;        % set to 1 for verbose staircase output

% class constructor function returns staircase object
sc = staircase(logStimLevels, initStepSize, stepSize,         ...
               nRightToDescend, nWrongToAscend, maxNumTrials, ...
               maxNumReversals, startLevel, verbose);
           
% setting the simulated noise level of the model subject for this demo
SIM.simNoiseStdDev = sqrt(2); % performance-limiting noise in simulation
           
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

% export data table from staircase object into a csv file
csvFileName = 'staircase_sim_output_table.csv';
csvF = fopen(csvFileName, 'w');
fprintf(csvF, 'logLev,nTrials,nCorrect');
for i = 1:length(sc.levels)
    fprintf(csvF, '\n%0.6f,%0.0f,%0.0f', sc.levels(i), ...
            sc.nTrials(i), sc.nCorrect(i));
end
fclose(csvF);

fprintf('Output saved in: %s\n', csvFileName)

return

function isSimCorrect = do_sim(SIM, linStimLevel)
    % Alex S Baldwin, McGill Vision Research, July 2019
    % Simulate subject behaviour for 2AFC task

    % simulated "noisy internal responses" from the two intervals
    respT = randn * SIM.simNoiseStdDev + linStimLevel; % target interval
    respN = randn * SIM.simNoiseStdDev;                % null interval
    
    isSimCorrect = respT > respN; % correct if noisy target > null

return