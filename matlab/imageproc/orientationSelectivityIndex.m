function OSI = orientationSelectivityIndex(tuningCurve,FRnonlin)

    if FRnonlin
        % using parameter values from Carandini, 2004, PLoS Biol
        tuningCurve = 10 * max(tuningCurve - 1.9,0).^(1.1);
    end
        
    thetaVals = linspace(0,180,length(tuningCurve));
    [maxResp,maxIndex] = max(tuningCurve);
    prefOrient = tuningCurve(maxIndex);
    orthOrient = prefOrient + 90;
    if orthOrient > 180
        orthOrient = orthOrient - 180;
    end
    [~,orthIndex] = min(abs(thetaVals - orthOrient));
    orthResp = tuningCurve(orthIndex);
    
    OSI = min((maxResp-orthResp)/maxResp,1);
end