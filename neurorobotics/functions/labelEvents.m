function lbl = labelEvents(h, sizes)
% LABELEVENTS Convert the event table into sample-wise label vectors.
%
%   lbl = labelEvents(h, sizes)
%
%   h     : header with h.EVENT.POS/TYP/DUR (e.g. from loadRuns)
%   sizes : [runs x 1] samples of each run (from loadRuns)
%
%   lbl   : struct with [samples x 1] vectors
%       Tk  trial number INSIDE the run (restarts from 1 at each run, 0 outside trials)
%       Fk  fixation (786)
%       Ak  cue (771 both feet, 773 both hands, 783 rest)
%       CFk continuous feedback (781)
%       Xk  outcome (897 hit, 898 miss) - only in online runs
%       Rk  run index (also outside trials)
%
%   A trial goes from its start to the end of the continuous feedback.
%   Trial start = event 1 (offline runs) or fixation 786 when there is no
%   event 1 (online runs).

    CUES = [771 773 783];
    OUTCOMES = [897 898];

    nSamples = sum(sizes);
    Tk  = zeros(nSamples, 1);
    Fk  = zeros(nSamples, 1);
    Ak  = zeros(nSamples, 1);
    CFk = zeros(nSamples, 1);
    Xk  = zeros(nSamples, 1);
    Rk  = zeros(nSamples, 1);

    last_pos = 1;
    offsets = cumsum(sizes);
    for i = 1:numel(sizes)
        Rk(last_pos:offsets(i)) = i;
        last_pos = offsets(i) + 1;
    end

    POS = h.EVENT.POS;
    TYP = h.EVENT.TYP;
    DUR = h.EVENT.DUR;

    trial = 0;
    trialOpen = false; % true between trial start and end of feedback
    startTrialIdx = NaN;

    for i = 1:numel(POS)
        startIdx = POS(i);
        endIdx = startIdx + DUR(i) - 1; % interval includes first and last sample

        assert(endIdx <= nSamples && Rk(startIdx) == Rk(max(endIdx, startIdx)), ...
            'Event %d (type %d) crosses a run boundary', i, TYP(i));

        if i > 1 && Rk(POS(i-1)) ~= Rk(startIdx) % new run -> restart numbering
            trial = 0;
            trialOpen = false;
        end

        switch TYP(i)
            case 1 % start
                trial = trial + 1;
                startTrialIdx = startIdx;
                trialOpen = true;
            case 786 % fixation
                if ~trialOpen % online runs have no event 1
                    trial = trial + 1;
                    startTrialIdx = startIdx;
                    trialOpen = true;
                end
                Fk(startIdx:endIdx) = TYP(i);
            case num2cell(CUES) % cue
                Ak(startIdx:endIdx) = TYP(i);
            case 781 % continuous feedback
                CFk(startIdx:endIdx) = TYP(i);
                assert(trialOpen, 'Feedback at event %d without a trial start', i);
                Tk(startTrialIdx:endIdx) = trial;
                trialOpen = false;
            case num2cell(OUTCOMES) % hit / miss
                Xk(startIdx:endIdx) = TYP(i);
            otherwise
                error('Weird stuff in GDFs: unknown event type %d (event %d)', TYP(i), i);
        end
    end

    lbl = struct('Tk', Tk, 'Fk', Fk, 'Ak', Ak, 'CFk', CFk, 'Xk', Xk, 'Rk', Rk);
end
