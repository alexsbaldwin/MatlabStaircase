classdef staircase < handle
% Alex S Baldwin - McGill Vision Research - April 2014 - Updated August 2017
% Implementation of a staircase class (object-oriented)
% Create a staircase by calling SC = staircase(...with arguments...)
% Arguments are described under the class constructor function below
% Current stimulus level to test is SC.curLevel
% After a response call SC.doResp(isCorrect) with a Boolean (0 or 1)
% Use SC.isFinished to find out whether the staircase is finished
% You can get your data from SC.levels, SC.nTrials, and SC.nCorrect
% Code downloaded from: https://github.com/alexsbaldwin/MatlabStaircase

    properties % properties of a class are like Matlab struct fields
        levels
        initStepSize
        stepSize
        rightRule
        wrongRule
        maxTrials
        maxRevs
        nTrials
        trialCount   % this is only used in staircase tracking
        nCorrect
        revCount
        curLevel
        curDirection % 0 is descending, 1 is ascending
        curRight
        curWrong
        verbose
        reversals
    end

    properties (Dependent) % these are defined dynamically (see below)
        nLevels
        minLevel
        maxLevel
        curIndex
        curStepSize
        isFinished
        curReversalThresh
        curReversalError
    end

    methods % public (callable) methods (class-specific functions)
        % class constructor function is staircase(...)
        function sc = staircase(levels, initStepSize, stepSize, ...
                                rightRule, wrongRule, maxTrials, ...
                                maxRevs, startLevel, verbose)
            if (nargin > 0)
                sc.levels       = levels;       % staircase levels
                sc.initStepSize = initStepSize; % step size before 1st rev
                sc.stepSize     = stepSize;     % subsequent step size
                sc.rightRule    = rightRule;    % n right to make harder
                sc.wrongRule    = wrongRule;    % n wrong to make easier
                sc.maxTrials    = maxTrials;    % max trials to terminate
                sc.maxRevs      = maxRevs;      % max revs to terminate
                sc.verbose      = verbose;      % 1 for console output
                sc.nTrials      = zeros(sc.nLevels,1);
                sc.nCorrect     = zeros(sc.nLevels,1);
                sc.trialCount   = 0;
                sc.revCount     = 0;
                sc.curDirection = 0;
                sc.curRight     = 0; % number of right responses
                sc.curWrong     = 0; % number of wrong responses
                sc.reversals    = [];
                sc.setToNearestLevel(startLevel);
            end
            if sc.verbose
                fprintf('\nStaircase class object initialised.\n')
            end
        end

        function minLevel = get.minLevel(sc) % minimum staircase level
            minLevel = min(sc.levels);
        end

        function maxLevel = get.maxLevel(sc) % maximum staircase level
            maxLevel = max(sc.levels);
        end

        function nLevels = get.nLevels(sc)  % number of staircase levels
            nLevels = length(sc.levels);
        end

        function curIndex = get.curIndex(sc) % current staircase index
            curIndex = find(sc.levels == sc.curLevel);
        end

        function curStepSize = get.curStepSize(sc) % current step size
            if sc.revCount == 0
                curStepSize = sc.initStepSize;
            else
                curStepSize = sc.stepSize;
            end
        end

        function isFinished = get.isFinished(sc) % check if finished
            isDoneRevs   = sc.revCount>=sc.maxRevs;
            isDoneTrials = sc.trialCount>=sc.maxTrials;
        	isFinished   = isDoneRevs||isDoneTrials;
        end

        function curReversalThresh = get.curReversalThresh(sc)
            if length(sc.reversals) < 3
                curReversalThresh = NaN;
            else
                curReversalThresh = mean(sc.reversals(2:end));
            end
        end

        function curReversalError = get.curReversalError(sc)
            if length(sc.reversals) < 3
                curReversalError = NaN;
            else
                curReversalError = std(sc.reversals(2:end)) / sqrt(length(sc.reversals(2:end)));
            end
        end

        function sc = doResp(sc, isCorrect)
           sc.trialInc;
           switch isCorrect
                case 0
                    sc.curWrong = sc.curWrong + 1;
                    if sc.verbose
                        fprintf('Incorrect response, nWrong = %0.0f\n', ...
                                 sc.curWrong)
                    end
                case 1
                    sc.curRight = sc.curRight + 1;
                    sc.nCorrect(sc.curIndex) = sc.nCorrect(sc.curIndex)+1;
                    if sc.verbose
                        fprintf('Correct response, nRight = %0.0f\n', ...
                                 sc.curRight)
                    end
           end
           sc.checkRules;
        end
    end

    methods (Access = private) % private (not callable) methods

        function sc = trialInc(sc)
            sc.nTrials(sc.curIndex) = sc.nTrials(sc.curIndex) + 1;
            if sc.revCount == 0
                if sc.verbose
                    fprintf('No reversals yet so do not count trials\n')
                end
            else
                sc.trialCount = sc.trialCount+1; % increment the counter
                if sc.verbose                    % only after 1st reversal
                    fprintf('Reversal: %0.0f, Trial: %0.0f\n', ...
                             sc.revCount,  sc.trialCount)
                end
            end
        end

        function sc = setToNearestLevel(sc,reqLevel)
            absDiff = abs(sc.levels - reqLevel);
            sc.curLevel = sc.levels(absDiff==min(absDiff));
        end

        function sc = checkRules(sc) % inc/dec staircase if rules exceeded
            if sc.curRight >= sc.rightRule
                if sc.verbose
                    fprintf('Change level: decrease.\n')
                end
                sc.changeLevel(0);

            elseif sc.curWrong >= sc.wrongRule
                if sc.verbose
                    fprintf('Change level: increase.\n')
                end
                sc.changeLevel(1);
            end
        end

        function sc = changeLevel(sc, newDir)

            if sc.curDirection ~= newDir  % if changed direction then
                sc.doReversal(newDir);    % count reversal
            end

            switch sc.curDirection % change staircase level
                case 0
                    sc.curLevel = sc.curLevel - sc.curStepSize;
                case 1
                    sc.curLevel = sc.curLevel + sc.curStepSize;
            end

            if sc.curLevel > sc.maxLevel     % keep staircase level
                sc.curLevel = sc.maxLevel;   % within limits
                if sc.verbose
                    fprintf('Max exceeded, bound at %0.0f.\n', sc.curLevel)
                end
            elseif sc.curLevel < sc.minLevel
                sc.curLevel = sc.minLevel;
                if sc.verbose
                    fprintf('Min exceeded, bound at %0.0f.\n', sc.curLevel)
                end
            else
                if sc.verbose
                    fprintf('Level changed to %0.0f.\n', sc.curLevel)
                end
            end

            sc.curRight = 0;
            sc.curWrong = 0;

        end

        function sc = doReversal(sc, newDir) % count reversal and
        	sc.revCount = sc.revCount + 1;   % change direction
            sc.reversals(end+1) = sc.curLevel;
            if sc.verbose
                fprintf('Reversal, %0.0f to %0.0f, revCount = %0.0f, revThresh = %0.1f +/- %0.1f\n', ...
                        sc.curDirection, newDir, sc.revCount, sc.curReversalThresh, sc.curReversalError)
            end
            sc.curDirection = newDir;
        end

    end

end