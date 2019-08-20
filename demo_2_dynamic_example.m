function demo_2_dynamic_example
% Alex S Baldwin, McGill Vision Research, July 2019
% Demo of object-oriented psychophysical staircase implementation. Requires
% staircase.m class definition file to be in same folder or on Matlab path.
% This example shows the "live" progress of the staircase in a simulated 
% 2AFC task. The intensity units are given in dB, as might be used for 
% contrast detection. From: https://github.com/alexsbaldwin/MatlabStaircase

close all

logStimLevels   = -18:3:24; % stimulus levels (here in dB logarithmic units)
initStepSize    = 6;        % staircase step size before first reversal
stepSize        = 3;        % staircase step size after first reversal
nRightToDescend = 3;        % number of "correct" responses before descending
nWrongToAscend  = 1;        % number of "incorrect" responses before ascending
maxNumTrials    = 80;       % maximum number of trials for staircase
maxNumReversals = 16;       % maximum number of reversals for staircase
startLevel      = 12;       % starting level of staircase (here in dB units)
verbose         = 1;        % set to 1 for verbose staircase output
           
% setting the simulated noise level of the model subject for this demo
SIM.simNoiseStdDev = sqrt(2); % performance-limiting noise in simulation
SIM.nSimulations   = 3;
SIM.trialPauseSec  = 0.01;

figure(1)
figpos = [200 200 600 600];
set(gcf, 'Units', 'pixels','PaperUnits', 'points', 'Position', ...
    figpos, 'PaperPosition', figpos, 'Color', [1 1 1]);
hold on
for iSim = 1:SIM.nSimulations

    % class constructor function returns staircase object
    sc = staircase(logStimLevels, initStepSize, stepSize,     ...
               nRightToDescend, nWrongToAscend, maxNumTrials, ...
               maxNumReversals, startLevel, verbose, 'resetting', 2);
    
    subplot(SIM.nSimulations,1,iSim)
    hold on
    axis([0,maxNumTrials+15,min(logStimLevels),max(logStimLevels)])
    xlabel('Trial number')
    ylabel('Test stimulus intensity')
    iTrial = 0;

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

        iTrial = iTrial + 1;
        if sc.didJustReverse
            if sc.curDirection == 0
                m = 'v'; mSize = 6;
            else
                m = '^'; mSize = 6;
            end
        else
            m = 'o'; mSize = 6;
        end
        
        if isSimCorrect
            plot(iTrial, logStimLevel, 'marker', m, 'markeredgecolor', [0,0,0], ...
                 'markerfacecolor', [1,1,1], 'markersize', mSize)
        else
            plot(iTrial, logStimLevel, 'marker', m, 'markeredgecolor', [0,0,0], ...
                'markerfacecolor', [0,0,0], 'markersize', mSize)
        end
        
        pause(SIM.trialPauseSec);

    end
    
    errorbar(iTrial+2,sc.curReversalThresh,sc.curReversalError, ...
             'marker','s', 'markersize',10,'markeredgecolor',[0,0,0], ...
             'markerfacecolor',[0,0,0],'color',[0,0,0])
    plot([0,iTrial+2],[sc.curReversalThresh,sc.curReversalThresh],...
         'color', [0,0,0], 'linestyle', '--')
    text(iTrial+4,sc.curReversalThresh,sprintf('%0.1f +/- %0.1f', ...
         sc.curReversalThresh,sc.curReversalError))
    
    % export data table from staircase object into a csv file
    csvFileName = sprintf('demo_2_staircase_sim_output_table_sim%0.0f.csv',iSim);
    csvF = fopen(csvFileName, 'w');
    fprintf(csvF, 'logLev,nTrials,nCorrect');
    for i = 1:length(sc.levels)
        fprintf(csvF, '\n%0.6f,%0.0f,%0.0f', sc.levels(i), ...
                sc.nTrials(i), sc.nCorrect(i));
    end
    fclose(csvF);
end

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