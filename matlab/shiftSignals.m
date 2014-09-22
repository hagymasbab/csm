function shifted_signals = shiftSignals(signals)
    % they should be rows
    step = (max(signals(:)) - min(signals(:))) / 100;
    shifted_signals = signals + repmat((0:size(signals,1)-1)*step,size(signals,2),1)';
end