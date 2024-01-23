classdef programState
    properties
        currentTime;
        tof_current;
        initial_tof;
        previousTime;
        N;
        initial_DeltaV;
        failedOrbits;
        testedOrbits;
        isRunning;
        axRelative;
        twMap;
        twMapInds;
    end
end

