function results = run_stewart_rainfall_continuation(constantRainfall,initialState)
%RUN_STEWART_RAINFALL_CONTINUATION Run Stewart model with optional continuation.
%
%   results = run_stewart_rainfall_continuation(constantRainfall)
%   results = run_stewart_rainfall_continuation(constantRainfall,initialState)
%
% constantRainfall is the annual water input in mm/year. Nitrogen and all
% other parameters are still read from input.dat.
%
% initialState is optional. When supplied, the simulation starts from the
% previous run's final field_species, mid_resource and deep_resource arrays.
% This allows continuation sweeps for testing hysteresis.
%
% The animation and plot_all_graphs output are disabled so this function can
% be called repeatedly by a rainfall-sweep script.
%
% ADDITIONAL DIAGNOSTICS:
%   - Checks for NaN/Inf after PART ONE and PART TWO.
%   - Records total water before PART ONE, after PART ONE and after PART TWO.
%   - Records equivalent nitrogen totals.
%   - Calculates change over the final 50 years.
%   - Returns final-window means for hysteresis calculations.
%
% IMPORTANT:
%   This version does NOT change the Stewart model equations or add/remove
%   any water-loss processes. It only diagnoses the existing model.

if nargin < 1 || nargin > 2 || ~isnumeric(constantRainfall) || ...
        ~isscalar(constantRainfall) || ~isfinite(constantRainfall) || ...
        constantRainfall < 0

    error('constantRainfall must be one finite, non-negative numeric value.');
end

if nargin < 2
    initialState = [];
end


%FULL VERSION - reads data from input file
%Multidirectional wind option
%
%Note - 'field' refers to the x & y rectangular grid, cells are 1m^2 elements
%
%In this version,
%   species 1 is grass (black gramma)
%   species 2 is shrub (creosotebush)
%   resource 1 is water
%   resource 2 is nitrogen


%---------------DEFINE VARIABLES-------------------------------------------

%use global to pass values to subroutines
%defined variables and arrays

global time fieldsize vectors randp species windir ...
       Growth_lim Bmax QV Maintenance Efficiency fail mortality raindata

global drought Reprod

% calculated variables and arrays

global field_species mid_resource Swater Swind Scow ...
       B_threshold deep_resource

global growth_rate loop counter

global time_series_plant time_series_resource ...
       max_connected_cells ave_connected_cells

% subroutines

global fn_r_and_p fn_change_in_biomass fn_plot_graphs


%---------------READ IN DATA-----------------------------------------------

%Read in parameters and data.
%NOTE - fscanf reads data column by row

fid = fopen('input.dat','r');

if fid == -1
    error('Could not open input.dat. MATLAB current folder is: %s',pwd);
end


fgetl(fid);
fgetl(fid);

time = fscanf(fid,'%d');


fgetl(fid);
fieldsize = fscanf(fid,'%d');


fgetl(fid);
vectors = fscanf(fid,'%d');


fgetl(fid);
randp = fscanf(fid,'%d');


fgetl(fid);
species = fscanf(fid,'%d');


fgetl(fid);
windir = fscanf(fid,'%d');


fgetl(fid);
gradient = fscanf(fid,'%d');


fgetl(fid);
Growth_lim = fscanf(fid,'%g',[species 1]);


fgetl(fid);
fgetl(fid);
Bmax = fscanf(fid,'%g',[species 1]);


fgetl(fid);
fgetl(fid);
QV = fscanf(fid,'%g',[randp 1]);


fgetl(fid);
fgetl(fid);
[M,count] = fscanf(fid,'%g',[species randp]);

Maintenance = M';


fgetl(fid);
fgetl(fid);
E = fscanf(fid,'%g',[species randp]);

Efficiency = E';


fgetl(fid);
fgetl(fid);
fail = fscanf(fid,'%g',[species 1]);


fgetl(fid);
fgetl(fid);
mortality = fscanf(fid,'%g',[species 1]);


fgetl(fid);
fgetl(fid);
drought = fscanf(fid,'%g',[species 1]);


fgetl(fid);
fgetl(fid);
[Reprod,count] = fscanf(fid,'%g',[vectors species]);

Reprod = Reprod';


fgetl(fid);
fgetl(fid);
fgetl(fid);

initial_species = fscanf(fid,'%g',[species 1]);


fgetl(fid);
fgetl(fid);

initial_randp = fscanf(fid,'%g',[randp 1]);


fgetl(fid);
fgetl(fid);

initial_deep_randp = fscanf(fid,'%g',[randp 1]);


fgetl(fid);
fgetl(fid);

Random_field_flag = fscanf(fid,'%d');


fgetl(fid);


%---------------CONSTANT RAINFALL RUN--------------------------------------

% The external rainfall record is not used in this function.

raindata = [];


%---------------RESOURCE SMOOSH MATRICES-----------------------------------

for i = 1:species

    fgetl(fid);

    temp = fscanf(fid,'%g',[vectors vectors]);

    Swater1(:,:,i) = temp';

end


fgetl(fid);
fgetl(fid);

temp = fscanf(fid,'%g',[vectors vectors]);

Swater2(:,:) = temp';


fgetl(fid);


for i = 1:species

    fgetl(fid);

    temp = fscanf(fid,'%g',[vectors vectors]);

    Swind(:,:,i) = temp';

end


fgetl(fid);


for i = 1:species

    fgetl(fid);

    temp = fscanf(fid,'%g',[vectors vectors]);

    Scow(:,:,i) = temp';

end


fgetl(fid);
fgetl(fid);

fclose(fid);


%form water smoosh arrays from input data and gradient
%
%NB Smoosh for wind and cows are independent of gradient

if gradient > 10

    Swater = zeros(3,3,2);

    Swater(3,2,1) = 1;
    Swater(3,2,2) = 1;


elseif gradient == 10

    Swater = Swater1;


elseif gradient == 0

    for i = 1:species
        Swater(:,:,i) = Swater2;
    end


else

    %Linear interpolation between the two smoosh matrices

    for i = 1:species

        Swater(:,:,i) = ...
            (Swater1(:,:,i) - Swater2(:,:)) * (gradient/10) ...
            + Swater2(:,:);

    end

end


%---------------INITIALISE MODEL STATE-------------------------------------

if isempty(initialState)

    % Original field_species initialisation.

    field_species = zeros(fieldsize,fieldsize,species);

    for i = 1:species

        field_species(:,:,i) = initial_species(i);

    end


    if Random_field_flag == 1

        RandomField = csvread('random100.dat',1,0);

        if size(RandomField,1) < fieldsize || ...
           size(RandomField,2) < fieldsize

            error('random100.dat is smaller than the requested field size.');

        end


        for column = 1:fieldsize

            for row = 1:fieldsize

                field_species(row,column,1) = ...
                    field_species(row,column,1) ...
                    * RandomField(row,column);

            end

        end

    end


    % Original resource-store initialisation.

    mid_resource = zeros(fieldsize,fieldsize,randp);

    for i = 1:randp

        mid_resource(:,:,i) = initial_randp(i);

    end


    deep_resource = zeros(fieldsize,fieldsize,randp);

    for i = 1:randp

        deep_resource(:,:,i) = initial_deep_randp(i);

    end


else

    requiredFields = { ...
        'field_species', ...
        'mid_resource', ...
        'deep_resource'};

    for i = 1:numel(requiredFields)

        if ~isfield(initialState,requiredFields{i})

            error( ...
                'initialState is missing the field "%s".', ...
                requiredFields{i});

        end

    end


    expectedSpeciesSize = [fieldsize fieldsize species];
    expectedResourceSize = [fieldsize fieldsize randp];


    if ~isequal(size(initialState.field_species),expectedSpeciesSize)

        error('initialState.field_species has the wrong size.');

    end


    if ~isequal(size(initialState.mid_resource),expectedResourceSize)

        error('initialState.mid_resource has the wrong size.');

    end


    if ~isequal(size(initialState.deep_resource),expectedResourceSize)

        error('initialState.deep_resource has the wrong size.');

    end


    field_species = initialState.field_species;
    mid_resource = initialState.mid_resource;
    deep_resource = initialState.deep_resource;

end


%---------------CHECK INITIAL STATE----------------------------------------

if any(~isfinite(field_species(:)))

    error( ...
        'Initial biomass state contains NaN or Inf at rainfall %.0f.', ...
        constantRainfall);

end


if any(~isfinite(mid_resource(:)))

    error( ...
        'Initial mid-resource state contains NaN or Inf at rainfall %.0f.', ...
        constantRainfall);

end


if any(~isfinite(deep_resource(:)))

    error( ...
        'Initial deep-resource state contains NaN or Inf at rainfall %.0f.', ...
        constantRainfall);

end


%---------------PRELIMINARY CALCULATIONS-----------------------------------

%These calculations are executed ONCE at start of program.
%Section creates arrays that are used in code but not read from input file.

time_series_plant = zeros(time+1,species*2);

time_series_plant(1,1) = 0;


for i = 1:species

    time_series_plant(1,i+1) = ...
        sum(sum(field_species(:,:,i))) ...
        /(fieldsize*fieldsize);

end


time_series_resource = zeros(time+1,randp*3);

time_series_resource(1,1) = 0;


for i = 1:randp

    time_series_resource(1,i+1) = ...
        sum(sum(mid_resource(:,:,i))) ...
        /(fieldsize*fieldsize);

    time_series_resource(1,i+3) = ...
        sum(sum(deep_resource(:,:,i))) ...
        /(fieldsize*fieldsize);

end


for i = 1:species

    %minimum biomass content for cell to be considered vegetated
    %10% of maximum biomass - original Stewart threshold

    B_threshold(i) = Bmax(i)*0.1;

end


%initialise arrays containing transect data for graphs
%transect is taken down center line of grid

growth_rate = zeros(fieldsize,fieldsize,species);


max_connected_cells = zeros(time+1,2);

ave_connected_cells = zeros(time+1,2);


con_count = zeros(fieldsize,2);


for i = 1:fieldsize

    if ((field_species(i,25,1) < B_threshold(1)) && ...
        (field_species(i,25,2) < B_threshold(2)))

        con_count(i,1) = 1;

    end


    if ((i == 1) && (con_count(1,1) == 1))

        con_count(i,2) = 1;

    end


    if ((i >= 2) && (con_count(i,1) == 1))

        con_count(i,2) = ...
            con_count(i-1,2) + con_count(i,1);

    end

end


[C,conn] = max(con_count(:,2));


max_connected_cells(1,1) = 0;

max_connected_cells(1,2) = con_count(conn,2);


ave_connected_cells(1,1) = 0;

ave_connected_cells(1,2) = ...
    (sum(con_count(:,1))/fieldsize)*100;


counter = 0;

clc;


%---------------ALTERNATIVE SMOOSH DESCRIPTIONS----------------------------

%This section was used to generate data for a paper
%can be safely deleted

%Swater=zeros(3,3,2);
%Swind=zeros(3,3,2);
%Scow=zeros(3,3,2);

%Pure advection
%Swater(2,2,1)=1;
%Swind(2,2,1)=1;
%Scow(2,2,1)=1;

%Swater(2,2,2)=1;
%Swind(2,2,2)=1;
%Scow(2,2,2)=1;

%local advection
%Swater(3,2,1)=1;
%Swind(3,2,1)=1;
%Scow(3,2,1)=1;

%Swater(3,2,2)=1;
%Swind(3,2,2)=1;
%Scow(3,2,2)=1;


%---------------OUTPUT HISTORIES FOR RAINFALL SWEEP------------------------

totalGrassHistory = zeros(time,1);

totalShrubHistory = zeros(time,1);

totalBiomassHistory = zeros(time,1);

totalWaterHistory = zeros(time,1);

totalNitrogenHistory = zeros(time,1);


% Water accounting

waterBeforePartOne = zeros(time,1);

waterAfterPartOne = zeros(time,1);

waterAfterPartTwo = zeros(time,1);


% Nitrogen accounting

nitrogenBeforePartOne = zeros(time,1);

nitrogenAfterPartOne = zeros(time,1);

nitrogenAfterPartTwo = zeros(time,1);


% Net changes caused during each section

waterChangePartOne = zeros(time,1);

waterChangePartTwo = zeros(time,1);

waterNetChange = zeros(time,1);


nitrogenChangePartOne = zeros(time,1);

nitrogenChangePartTwo = zeros(time,1);

nitrogenNetChange = zeros(time,1);


%---------------TIME LOOP - CALCULATE CHANGE OF RESOURCE AND BIOMASS-------

for loop = 1:time

    %Time loop - parameters recalculated each time step


    if loop == 1 || mod(loop,50) == 0 || loop == time

        fprintf( ...
            'Rainfall %.0f mm/year: model year %d of %d\n', ...
            constantRainfall,loop,time);

    end


    counter = counter + 1;


    %-----------MODIFY INPUT VARIABLES IF APPROPRIATE----------------------

    % Use the same constant rainfall in every simulated year.

    QV(1) = constantRainfall;


    %Example drought code from original model:

%     if ((loop>50)&&(loop<57))
%
%         QV(1)=25;
%
%     else
%
%         QV(1)=raindata(loop,2);
%
%     end


    %Example disturbance code:

%     if ((loop>220)&&(loop<222))
%
%         field_species(:,:,2)=0.1;
%
%     end


    %-----------RESOURCE TOTAL BEFORE PART ONE-----------------------------

    waterBeforePartOne(loop) = ...
        sum(mid_resource(:,:,1),'all') ...
        + sum(deep_resource(:,:,1),'all');


    nitrogenBeforePartOne(loop) = ...
        sum(mid_resource(:,:,2),'all') ...
        + sum(deep_resource(:,:,2),'all');


    %-----------CALL PART ONE----------------------------------------------

    %PART ONE
    %Calc resource movement under action of vectors.

    [fn_r_and_p] = part_one();


    %-----------CHECK STATE AFTER PART ONE---------------------------------

    if any(~isfinite(mid_resource(:))) || ...
       any(~isfinite(deep_resource(:))) || ...
       any(~isfinite(field_species(:)))

        error( ...
            ['NON-FINITE STATE AFTER PART ONE.\n' ...
             'Rainfall = %.0f mm/year\n' ...
             'Year = %d'], ...
             constantRainfall,loop);

    end


    waterAfterPartOne(loop) = ...
        sum(mid_resource(:,:,1),'all') ...
        + sum(deep_resource(:,:,1),'all');


    nitrogenAfterPartOne(loop) = ...
        sum(mid_resource(:,:,2),'all') ...
        + sum(deep_resource(:,:,2),'all');


    %-----------CALL PART TWO----------------------------------------------

    %PART TWO
    %Calc use of resource by biomass, move remainder to stores.

    [fn_change_in_biomass] = part_two();


    %-----------CHECK STATE AFTER PART TWO---------------------------------

    if any(~isfinite(mid_resource(:))) || ...
       any(~isfinite(deep_resource(:))) || ...
       any(~isfinite(field_species(:)))

        error( ...
            ['NON-FINITE STATE AFTER PART TWO.\n' ...
             'Rainfall = %.0f mm/year\n' ...
             'Year = %d'], ...
             constantRainfall,loop);

    end


    waterAfterPartTwo(loop) = ...
        sum(mid_resource(:,:,1),'all') ...
        + sum(deep_resource(:,:,1),'all');


    nitrogenAfterPartTwo(loop) = ...
        sum(mid_resource(:,:,2),'all') ...
        + sum(deep_resource(:,:,2),'all');


    %-----------CALCULATE RESOURCE CHANGES---------------------------------

    waterChangePartOne(loop) = ...
        waterAfterPartOne(loop) ...
        - waterBeforePartOne(loop);


    waterChangePartTwo(loop) = ...
        waterAfterPartTwo(loop) ...
        - waterAfterPartOne(loop);


    waterNetChange(loop) = ...
        waterAfterPartTwo(loop) ...
        - waterBeforePartOne(loop);


    nitrogenChangePartOne(loop) = ...
        nitrogenAfterPartOne(loop) ...
        - nitrogenBeforePartOne(loop);


    nitrogenChangePartTwo(loop) = ...
        nitrogenAfterPartTwo(loop) ...
        - nitrogenAfterPartOne(loop);


    nitrogenNetChange(loop) = ...
        nitrogenAfterPartTwo(loop) ...
        - nitrogenBeforePartOne(loop);


    %-----------TOTALS OVER FULL GRID--------------------------------------

    grassField = field_species(:,:,1);

    shrubField = field_species(:,:,2);

    combinedBiomass = grassField + shrubField;


    waterField = ...
        mid_resource(:,:,1) ...
        + deep_resource(:,:,1);


    nitrogenField = ...
        mid_resource(:,:,2) ...
        + deep_resource(:,:,2);


    totalGrassHistory(loop) = ...
        sum(grassField,'all');


    totalShrubHistory(loop) = ...
        sum(shrubField,'all');


    totalBiomassHistory(loop) = ...
        sum(combinedBiomass,'all');


    totalWaterHistory(loop) = ...
        sum(waterField,'all');


    totalNitrogenHistory(loop) = ...
        sum(nitrogenField,'all');


    %-----------STORE RESULTS----------------------------------------------

    time_series_plant(loop+1,1) = loop;


    for i = 1:species

        time_series_plant(loop+1,i+1) = ...
            sum(sum(field_species(:,:,i))) ...
            /(fieldsize*fieldsize);

    end


    time_series_resource(loop+1,1) = loop;


    time_series_resource(loop+1,2) = ...
        sum(sum(mid_resource(:,:,1))) ...
        /(fieldsize*fieldsize);


    time_series_resource(loop+1,3) = ...
        sum(sum(mid_resource(:,:,2))) ...
        /(fieldsize*fieldsize);


    time_series_resource(loop+1,4) = ...
        sum(sum(deep_resource(:,:,1))) ...
        /(fieldsize*fieldsize);


    time_series_resource(loop+1,5) = ...
        sum(sum(deep_resource(:,:,2))) ...
        /(fieldsize*fieldsize);


    %-----------CONNECTIVITY CALCULATION-----------------------------------

    con_count = zeros(fieldsize,2);


    for i = 1:fieldsize

        if ((field_species(i,25,1) < B_threshold(1)) && ...
            (field_species(i,25,2) < B_threshold(2)))

            con_count(i,1) = 1;

        end


        if ((i == 1) && (con_count(1,1) == 1))

            con_count(1,2) = 1;

        end


        if ((i >= 2) && (con_count(i,1) == 1))

            con_count(i,2) = ...
                con_count(i-1,2) ...
                + con_count(i,1);

        end

    end


    [C,conn] = max(con_count(:,2));


    max_connected_cells(loop+1,1) = loop;

    max_connected_cells(loop+1,2) = ...
        con_count(conn,2);


    ave_connected_cells(loop+1,1) = loop;

    ave_connected_cells(loop+1,2) = ...
        (sum(con_count(:,1))/fieldsize)*100;


    % plot_all_graphs deliberately disabled during parameter sweeps.

end


%---------------RETURN LONG-TERM RESPONSE----------------------------------

finalWindow = min(50,time);

finalIndices = ...
    (time-finalWindow+1):time;


results.rainfall = constantRainfall;

results.time = (1:time)';


% Complete histories

results.totalGrassHistory = ...
    totalGrassHistory;

results.totalShrubHistory = ...
    totalShrubHistory;

results.totalBiomassHistory = ...
    totalBiomassHistory;

results.totalWaterHistory = ...
    totalWaterHistory;

results.totalNitrogenHistory = ...
    totalNitrogenHistory;


%---------------RESOURCE ACCOUNTING OUTPUT---------------------------------

results.waterBeforePartOne = ...
    waterBeforePartOne;

results.waterAfterPartOne = ...
    waterAfterPartOne;

results.waterAfterPartTwo = ...
    waterAfterPartTwo;

results.waterChangePartOne = ...
    waterChangePartOne;

results.waterChangePartTwo = ...
    waterChangePartTwo;

results.waterNetChange = ...
    waterNetChange;


results.nitrogenBeforePartOne = ...
    nitrogenBeforePartOne;

results.nitrogenAfterPartOne = ...
    nitrogenAfterPartOne;

results.nitrogenAfterPartTwo = ...
    nitrogenAfterPartTwo;

results.nitrogenChangePartOne = ...
    nitrogenChangePartOne;

results.nitrogenChangePartTwo = ...
    nitrogenChangePartTwo;

results.nitrogenNetChange = ...
    nitrogenNetChange;


%---------------FINAL 50 YEAR SERIES---------------------------------------

results.grassSeries = ...
    totalGrassHistory(finalIndices);

results.shrubSeries = ...
    totalShrubHistory(finalIndices);

results.biomassSeries = ...
    totalBiomassHistory(finalIndices);

results.waterSeries = ...
    totalWaterHistory(finalIndices);

results.nitrogenSeries = ...
    totalNitrogenHistory(finalIndices);


%---------------FINAL SINGLE-YEAR VALUES-----------------------------------

results.finalGrass = ...
    totalGrassHistory(end);

results.finalShrub = ...
    totalShrubHistory(end);

results.finalBiomass = ...
    totalBiomassHistory(end);

results.finalWater = ...
    totalWaterHistory(end);

results.finalNitrogen = ...
    totalNitrogenHistory(end);


%---------------FINAL-WINDOW MEANS-----------------------------------------

results.meanFinalGrass = ...
    mean(results.grassSeries);

results.meanFinalShrub = ...
    mean(results.shrubSeries);

results.meanFinalBiomass = ...
    mean(results.biomassSeries);

results.meanFinalWater = ...
    mean(results.waterSeries);

results.meanFinalNitrogen = ...
    mean(results.nitrogenSeries);


%---------------FINAL-WINDOW CHANGE----------------------------------------

% Difference between first and last values of final 50-year window.
%
% Near zero means approximately stationary.
% Large positive/negative value means the variable is still changing.

results.finalGrassChange = ...
    results.grassSeries(end) ...
    - results.grassSeries(1);


results.finalShrubChange = ...
    results.shrubSeries(end) ...
    - results.shrubSeries(1);


results.finalBiomassChange = ...
    results.biomassSeries(end) ...
    - results.biomassSeries(1);


results.finalWaterChange = ...
    results.waterSeries(end) ...
    - results.waterSeries(1);


results.finalNitrogenChange = ...
    results.nitrogenSeries(end) ...
    - results.nitrogenSeries(1);


% Relative final-window changes

results.relativeBiomassChange = ...
    abs(results.finalBiomassChange) ...
    / max(abs(results.meanFinalBiomass),eps);


results.relativeWaterChange = ...
    abs(results.finalWaterChange) ...
    / max(abs(results.meanFinalWater),eps);


results.relativeNitrogenChange = ...
    abs(results.finalNitrogenChange) ...
    / max(abs(results.meanFinalNitrogen),eps);


%---------------FINAL RESOURCE FLUX DIAGNOSTICS----------------------------

results.meanFinalWaterChangePartOne = ...
    mean(waterChangePartOne(finalIndices));


results.meanFinalWaterChangePartTwo = ...
    mean(waterChangePartTwo(finalIndices));


results.meanFinalWaterNetChange = ...
    mean(waterNetChange(finalIndices));


results.meanFinalNitrogenChangePartOne = ...
    mean(nitrogenChangePartOne(finalIndices));


results.meanFinalNitrogenChangePartTwo = ...
    mean(nitrogenChangePartTwo(finalIndices));


results.meanFinalNitrogenNetChange = ...
    mean(nitrogenNetChange(finalIndices));


%---------------PRINT FINAL DIAGNOSTIC-------------------------------------

fprintf('\n');
fprintf('Rainfall %.0f mm/year final diagnostic\n',constantRainfall);

fprintf('   Mean final biomass       = %.6e\n', ...
    results.meanFinalBiomass);

fprintf('   Mean final water         = %.6e\n', ...
    results.meanFinalWater);

fprintf('   Mean final nitrogen      = %.6e\n', ...
    results.meanFinalNitrogen);

fprintf('   Biomass change last %d y = %+.6e\n', ...
    finalWindow,results.finalBiomassChange);

fprintf('   Water change last %d y   = %+.6e\n', ...
    finalWindow,results.finalWaterChange);

fprintf('   Nitrogen change last %d y= %+.6e\n', ...
    finalWindow,results.finalNitrogenChange);

fprintf('   Mean yearly PART 1 dW    = %+.6e\n', ...
    results.meanFinalWaterChangePartOne);

fprintf('   Mean yearly PART 2 dW    = %+.6e\n', ...
    results.meanFinalWaterChangePartTwo);

fprintf('   Mean yearly NET dW       = %+.6e\n', ...
    results.meanFinalWaterNetChange);

fprintf('\n');


%---------------FINAL FIELDS-----------------------------------------------

results.finalGrassField = ...
    field_species(:,:,1);

results.finalShrubField = ...
    field_species(:,:,2);

results.finalBiomassField = ...
    field_species(:,:,1) ...
    + field_species(:,:,2);

results.finalWaterField = ...
    mid_resource(:,:,1) ...
    + deep_resource(:,:,1);

results.finalNitrogenField = ...
    mid_resource(:,:,2) ...
    + deep_resource(:,:,2);


%---------------STATE FOR CONTINUATION-------------------------------------

results.finalState.field_species = ...
    field_species;

results.finalState.mid_resource = ...
    mid_resource;

results.finalState.deep_resource = ...
    deep_resource;


end