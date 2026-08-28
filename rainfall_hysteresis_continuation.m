clear;
clc;
clear all;

%% Rainfall sequences

% Manually chosen low-rainfall points followed by 25 mm increments.
forwardRainfall = [2 5 10 14 15 20 25 50:25:300];

% Exactly the same rainfall values in reverse.
backwardRainfall = fliplr(forwardRainfall);

numberForward = numel(forwardRainfall);
numberBackward = numel(backwardRainfall);


%% Storage

forwardBiomass = nan(1,numberForward);
forwardGrass = nan(1,numberForward);
forwardShrub = nan(1,numberForward);

backwardBiomass = nan(1,numberBackward);
backwardGrass = nan(1,numberBackward);
backwardShrub = nan(1,numberBackward);


%% Live hysteresis figure

hysteresisFigure = figure('Color','w');

hold on;
grid on;

forwardLine = plot( ...
    NaN,NaN,'-ob', ...
    'LineWidth',2, ...
    'MarkerFaceColor','b', ...
    'DisplayName','Increasing rainfall');

backwardLine = plot( ...
    NaN,NaN,'-or', ...
    'LineWidth',2, ...
    'MarkerFaceColor','r', ...
    'DisplayName','Decreasing rainfall');

xlabel('Rainfall (mm/year)');
ylabel('Mean final total biomass (g)');
title('Rainfall hysteresis test');

legend('Location','best');

xlim([min(forwardRainfall) max(forwardRainfall)]);


%% Biomass time-series figure

timeSeriesFigure = figure('Color','w');

timeLayout = tiledlayout(2,1);

title( ...
    timeLayout, ...
    'Biomass Through Time for Each Rainfall Simulation');


% Forward sweep panel

forwardTimeAxis = nexttile;

hold(forwardTimeAxis,'on');
grid(forwardTimeAxis,'on');

xlabel(forwardTimeAxis,'Year');
ylabel(forwardTimeAxis,'Total Biomass (g)');

title( ...
    forwardTimeAxis, ...
    'Forward Sweep');


% Backward sweep panel

backwardTimeAxis = nexttile;

hold(backwardTimeAxis,'on');
grid(backwardTimeAxis,'on');

xlabel(backwardTimeAxis,'Year');
ylabel(backwardTimeAxis,'Total Biomass (g)');

title( ...
    backwardTimeAxis, ...
    'Backward Sweep');


%% Forward continuation sweep

state = [];

for k = 1:numberForward
    rainfall = forwardRainfall(k);
        results = run_stewart_rainfall_continuation( ...
                rainfall,state);


    %--------------------------------------------------------------
    % HYSTERESIS VALUES
    %
    % Use final-window means rather than one single final timestep.
    %--------------------------------------------------------------

    forwardBiomass(k) = ...
        results.meanFinalBiomass;

    forwardGrass(k) = ...
        results.meanFinalGrass;

    forwardShrub(k) = ...
        results.meanFinalShrub;

    %--------------------------------------------------------------
    % BIOMASS TIME SERIES
    %--------------------------------------------------------------

    plot( ...
        forwardTimeAxis, ...
        results.time, ...
        results.totalBiomassHistory, ...
        'LineWidth',1.2, ...
        'DisplayName', ...
        sprintf('%.0f mm/year',rainfall));


    legend( ...
        forwardTimeAxis, ...
        'show', ...
        'Location','eastoutside');


    %--------------------------------------------------------------
    % PASS FINAL STATE TO NEXT RAINFALL VALUE
    %--------------------------------------------------------------

    state = results.finalState;


    %--------------------------------------------------------------
    % UPDATE HYSTERESIS CURVE
    %--------------------------------------------------------------

    set( ...
        forwardLine, ...
        'XData',forwardRainfall(1:k), ...
        'YData',forwardBiomass(1:k));


    drawnow;

end

    % Start backward branch from final high-rainfall state.

    backwardState = state;


    %% Backward continuation sweep

    for k = 1:numberBackward
        rainfall = backwardRainfall(k);
            results = run_stewart_rainfall_continuation( ...
                    rainfall,backwardState);

        %----------------------------------------------------------
        % HYSTERESIS VALUES
        %----------------------------------------------------------

        backwardBiomass(k) = ...
            results.meanFinalBiomass;

        backwardGrass(k) = ...
            results.meanFinalGrass;

        backwardShrub(k) = ...
            results.meanFinalShrub;


        %----------------------------------------------------------
        % BIOMASS TIME SERIES
        %----------------------------------------------------------

        plot( ...
            backwardTimeAxis, ...
            results.time, ...
            results.totalBiomassHistory, ...
            'LineWidth',1.2, ...
            'DisplayName', ...
            sprintf('%.0f mm/year',rainfall));


        legend( ...
            backwardTimeAxis, ...
            'show', ...
            'Location','eastoutside');


        %----------------------------------------------------------
        % PASS FINAL STATE TO NEXT LOWER RAINFALL
        %----------------------------------------------------------

        backwardState = results.finalState;


        %----------------------------------------------------------
        % UPDATE HYSTERESIS CURVE
        %----------------------------------------------------------

        set( ...
            backwardLine, ...
            'XData',backwardRainfall(1:k), ...
            'YData',backwardBiomass(1:k));


        drawnow;

    end



%% Put backward values in increasing-rainfall order

backwardBiomassAscending = ...
    fliplr(backwardBiomass);

backwardGrassAscending = ...
    fliplr(backwardGrass);

backwardShrubAscending = ...
    fliplr(backwardShrub);


%% Save results
hysteresisResults = table( ...
    forwardRainfall(:), ...
    forwardBiomass(:), ...
    backwardBiomassAscending(:), ...
    forwardGrass(:), ...
    backwardGrassAscending(:), ...
    forwardShrub(:), ...
    backwardShrubAscending(:), ...
    'VariableNames',{ ...
    'Rainfall_mm_per_year', ...
    'ForwardTotalBiomass_g', ...
    'BackwardTotalBiomass_g', ...
    'ForwardGrass_g', ...
    'BackwardGrass_g', ...
    'ForwardShrub_g', ...
    'BackwardShrub_g'});

writetable( ...
    hysteresisResults, ...
    'rainfall_hysteresis_results.csv');


save( ...
    'rainfall_hysteresis_results.mat', ...
    'hysteresisResults', ...
    'forwardRainfall', ...
    'backwardRainfall', ...
    'forwardBiomass', ...
    'backwardBiomass', ...
    'forwardGrass', ...
    'backwardGrass', ...
    'forwardShrub', ...
    'backwardShrub');


