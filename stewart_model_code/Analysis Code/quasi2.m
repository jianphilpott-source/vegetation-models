function quasi2(constant_rainfall)
close all;
clc;

%FULL VERSION - reads data from input file 
%Multidirectional wind option
%Note - 'field' refers to the x & y rectangular grid, cells are 1m^2 elements
%In this version, 
%   species 1 is grass (black gramma) and species 2 is shrub (creosotebush)
%   resource 1 is water and resource 2 is nitrogen
%---------------DEFINE VARIABLES-------------------------------------------
%use global to pass values to subroutines
%defined variables and arrays
global time fieldsize vectors randp species windir Growth_lim Bmax QV Maintenance Efficiency fail mortality raindata
global drought Reprod 
% calculated variables and arrays
global field_species mid_resource Swater Swind Scow B_threshold deep_resource
global growth_rate loop counter     
global time_series_plant time_series_resource max_connected_cells ave_connected_cells
% subroutines
global fn_r_and_p fn_change_in_biomass fn_plot_graphs
%---------------READ IN DATA-----------------------------------------------
%Read in parameters and data. NOTE - fscanf reads data column by row 
fid = fopen ('input.dat'); %contains all defined variables and parameters
fgetl(fid);fgetl(fid);%reads line with parameter description and ignores blank lines
time = fscanf (fid, '%d');                          %integer value for years of simulation
fgetl(fid); fieldsize = fscanf (fid, '%d');         %field is fieldsize m^2
fgetl(fid); vectors = fscanf (fid, '%d');           %number of vectors (1=water, 2=wind, 3=cows)
fgetl(fid); randp = fscanf (fid, '%d');             %number of resources
fgetl(fid); species = fscanf (fid, '%d');           %number of plant species
fgetl(fid); windir = fscanf (fid, '%d');            %wind direction 1=downhill, 2=uphill, 3=L-to-R, 4=R-to-L 
fgetl(fid); gradient = fscanf (fid, '%d');          %reads in gradient of hill for smoosh calcs (arbitrary value, where 10 represents a maximum gradient for which the parameters apply)    
fgetl(fid); Growth_lim = fscanf (fid, '%g', [species 1]);    %Limiting growth rate for each species as percent increase in biomass
fgetl(fid);fgetl(fid); Bmax = fscanf (fid, '%g', [species 1]); %Maximum biomass for each species (g/m^2)- derived from 'ecotone'
fgetl(fid);fgetl(fid); QV = fscanf (fid, '%g', [randp 1]);    %Vertical Flux - vertical addition of resource (g/m2.yr)
fgetl(fid);fgetl(fid); [M, count] = fscanf (fid, '%g', [species randp]);     %maintenence term (g resource required to maintain each g of biomass)
    Maintenance=M';      % Matlab reads in arrays column by row - corrected here
fgetl(fid);fgetl(fid); E = fscanf (fid, '%g', [species randp]);     %Efficiency term (g resource per eahc gram of new growth of biomass)
    Efficiency=E';       % Matlab reads in arrays column by row - corrected here
fgetl(fid);fgetl(fid); fail = fscanf (fid, '%g', [species 1]);             %failure rate - percentage reduction of propagules
fgetl(fid);fgetl(fid); mortality = fscanf (fid, '%g', [species 1]);        %plant mortality rate - percentage of biomass that dies each year
fgetl(fid);fgetl(fid); drought = fscanf (fid, '%g', [species 1]);            %drought resistance of species -percentage term to relax quantity of biomass killed by water drought
fgetl(fid);fgetl(fid); [Reprod, count] = fscanf (fid, '%g', [vectors species]);  %percentage of propagules moved by each vector (derived loosely from USFAS website)
    Reprod=Reprod';      % Matlab reads in arrays column by row - corrected here
fgetl(fid);fgetl(fid);fgetl(fid);initial_species = fscanf (fid, '%g', [species 1]); %Initial vegetation cover (g/m2) of each i species (from John and Tony)
fgetl(fid);fgetl(fid);   % Initial values of randp, (for each resource) corresponding to mid soil layer store (model is insensitive to any realistic value set here)
    initial_randp = fscanf (fid, '%g', [randp 1 ]);         % Matlab reads in arrays column by row - corrected here
fgetl(fid);fgetl(fid);   % Initial values of randp, (for each resource) corresponding to deep soil layer store
    initial_deep_randp = fscanf (fid, '%g', [randp 1 ]);    % Matlab reads in arrays column by row - corrected here  (model is insensitive to realistic numbers) 
fgetl(fid);fgetl(fid); Random_field_flag=fscanf (fid, '%d'); %logical flag 1 = generate random species distribution, 0 = uniform distribution
fgetl(fid);
%---------------READ IN WATER INPUT (RAIN) DATA----------------------------
raindata=csvread ('raindat.dat',1,0); %read in yearly raindata (derived from tree-ring, real record, or set own pattern)
%---------------RESOURCE SMOOSH MATRICES-----------------------------------
for i=1:species     % Smoosh matrix for water i=1 for species 1, i=2 for species 2 etc
    fgetl(fid);
    temp = fscanf (fid, '%g', [vectors vectors]);
    Swater1(:,:,i)=temp';         % correct row-column   
end
fgetl(fid);fgetl(fid); %water smoosh for flat field only (gradient = 0)
    temp = fscanf (fid, '%g', [vectors vectors]);
    Swater2(:,:)=temp';           % correct row-column   
fgetl(fid);
for i=1:species     % Smoosh matrix for wind
    fgetl(fid);
    temp = fscanf (fid, '%g', [vectors vectors]);
    Swind(:,:,i)=temp';           % correct row-column  
end
fgetl(fid);
for i=1:species     % Smosh matrix for cows
    fgetl(fid);
    temp = fscanf (fid, '%g', [vectors vectors]);
    Scow(:,:,i)=temp';             % correct row-column  
end
fgetl(fid);fgetl(fid);
%form water smoosh arrays from input data and gradient - NB Smoosh for wind and cows are independent of gradient
if (gradient>10)
    Swater=zeros(3,3,2); %resource smoosh
    Swater(3,2,1)=1;Swater(3,2,2)=1; %resource smoosh
elseif (gradient==10), Swater=Swater1; %Smoosh array for water is dependent on gradient (0=flat, 10=represents gradient that corresponds to defininition of smoosh matrix) 
elseif (gradient==0)
    for i=1:species,Swater(:,:,i)=Swater2;end,
else        %Linear interpolation between the two smoosh matrices
    for i=1:species,Swater(:,:,i)=(Swater1(:,:,i)-Swater2(:,:))*(gradient/10)+Swater2(:,:);end
end
%---------------POPULATE FIElD_SPECIES ARRAY-------------------------------
field_species=zeros(fieldsize,fieldsize,species); %3 dimensional array - x and y are field size, each z is biomass of EACH species
for i=1:species %populate field_species array
    field_species(:,:,i)=initial_species(i); %field_species contains biomass (g/m2) for each species in cell
end                                          
if (Random_field_flag==1)  %Calculation of 'RandomField' data was done ouside of this code - its a matrix of random numbers approximating white noise
    RandomField=csvread ('random100.dat',1,0); %stored data will only work for fields of up to 100 x 100
    %for i=1:species  % Use this loop if initial distribution of species etc is also to be random, a pointless excersise if intital value of
    %shrubs represents a seed source and not established plants(i.e.biomass<B_threshold)
        for column=1:fieldsize
            for row=1:fieldsize
               field_species(row,column,1)=field_species(row,column,1)*RandomField(row,column); %multiply initial species by random value (random value is a %)
            end %end of row loop
        end %end of column loop
    %end %end of species for loop
end %end of if test  
%---------------POPULATE mid_resource ARRAY---------------------------------
%3 dimensional array - x and y are field size, each z is mass of EACH resource mid_resource in mid layer resource store, deep_resource is deep soil layer
%store for resources and layers are defined to match the two plant species being considered.  Resources are added to and remove from these stores by
%the action of the plants (see part_two()) and NOT by vectors.  
mid_resource=zeros(fieldsize,fieldsize,randp); %Mid layer resource store (name is a legacy that I am too lazy to change)
for i=1:randp
    mid_resource(:,:,i)=initial_randp(i);  % mid_resource contains mass (g/m2) of resource in each cell
end                                       
deep_resource=zeros(fieldsize,fieldsize,randp); %Deep layer resource store (deep is a name that means what is says it means)
for i=1:randp
    deep_resource(:,:,i)=initial_deep_randp(i);  % deep store contains mass (g/m2) of resource
end
%---------------PRELIMINARY CALCULATIONS- ---------------------------------
%These calculations are executed ONCE at start of program
%Section creates arrays that are used in code but not read from input file
time_series_plant=zeros(time+1,species*2); time_series_plant(1,1)=0; %used to store average biomass cover through time
for i=1:species %Sum both species in field for output graph - gives total biomass of each species in 'field' per itteration
    time_series_plant(1,i+1)=sum(sum(field_species(:,:,i)))/(fieldsize*fieldsize);
end
time_series_resource=zeros(time+1,randp*3); time_series_resource(1,1)=0; %used to store average remaining resource through time
for i=1:randp %sum resource for graph
    time_series_resource(1,i+1)=sum(sum(mid_resource(:,:,i)))/(fieldsize*fieldsize); %mid layer store
    time_series_resource(1,i+3)=sum(sum(deep_resource(:,:,i)))/(fieldsize*fieldsize);  %deep layer store
end
for i=1:species
    B_threshold(i)=Bmax(i)*0.1; %minimum biomass content for cell to be considered vegetated (i.e. 10% of max biomass - arbitrary number)
end
%initialise arrays containing transect data for graphs - transect is taken down center line of the grid
growth_rate=zeros(fieldsize,fieldsize,species);%initialise array to hold data for redistribution in biomass for each time step
max_connected_cells=zeros(time+1,2); ave_connected_cells=zeros(time+1,2);%rough connectivity data for mid line transect
con_count=zeros(fieldsize,2); %used later for counting number of connectied cells
%Note, column 1 corresponds to itteration (year number) and is set to zero here as the calculation is made before the start of the simulation
for i=1:fieldsize
    if ((field_species(i,25,1)<B_threshold(1))&&(field_species(i,25,2)<B_threshold(2)))
            con_count(i,1)=1; %labels cell as 1 if unvegetated
    end  
    if ((i==1)&&(con_count(1,1)==1))
        con_count(i,2)=1; %Connectivity count for first cell only
    end
    if ((i>=2)&&(con_count(i,1)==1)) %connectivity count for the rest of the cells
            con_count(i,2)=con_count(i-1,2)+con_count(i,1); %sums number of connected cells
    end
end
%matlab preserves values from previous executions in its arrays - hence the above 'zero' statements)
[C conn]=max(con_count(:,2));  %find the maximum number of connected cells
max_connected_cells(1,1)=0; max_connected_cells(1,2)=con_count(conn,2); %data arrays for connectivity counts
ave_connected_cells(1,1)=0; ave_connected_cells(1,2)=(sum(con_count(:,1))/fieldsize)*100;
counter=0; clc, %initialise couter and clear screen (watch out though - this doesn't clear previous graphs if they are still open)
%---------------ALTERNATIVE SMOOSH DESCRIPTIONS----------------------------
%This section was used to generate data for a paper - can be safely deleted
%Swater=zeros(3,3,2);Swind=zeros(3,3,2); Scow=zeros(3,3,2);
%Pure advection
%Swater(2,2,1)=1;Swind(2,2,1)=1; Scow(2,2,1)=1;
%Swater(2,2,2)=1;Swind(2,2,2)=1; Scow(2,2,2)=1;
%local advection
%Swater(3,2,1)=1;Swind(3,2,1)=1; Scow(3,2,1)=1;
%Swater(3,2,2)=1;Swind(3,2,2)=1; Scow(3,2,2)=1;

%% ---------------- INITIALISE ANIMATION ----------------

animationEvery = 1;

animationFigure = figure( ...
    'Name','Biomass and Resources', ...
    'NumberTitle','off', ...
    'Color','w');

animationLayout = tiledlayout(animationFigure,1,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% Initial data
combinedBiomass = sum(field_species,3);
totalWater = mid_resource(:,:,1) + deep_resource(:,:,1);
totalNitrogen = mid_resource(:,:,2) + deep_resource(:,:,2);

% Biomass axes
axBiomass = nexttile(animationLayout,1);
imageBiomass = imagesc(axBiomass,combinedBiomass);
axis(axBiomass,'image');
set(axBiomass,'YDir','normal');
colorbar(axBiomass);
title(axBiomass,'Combined biomass: year 0');
xlabel(axBiomass,'x cell');
ylabel(axBiomass,'y cell');
clim(axBiomass,[0 sum(Bmax)]);

% Water axes
axWater = nexttile(animationLayout,2);
imageWater = imagesc(axWater,totalWater);
axis(axWater,'image');
set(axWater,'YDir','normal');
colorbar(axWater);
title(axWater,'Total water: year 0');
xlabel(axWater,'x cell');
ylabel(axWater,'y cell');

% Nitrogen axes
axNitrogen = nexttile(animationLayout,3);
imageNitrogen = imagesc(axNitrogen,totalNitrogen);
axis(axNitrogen,'image');
set(axNitrogen,'YDir','normal');
colorbar(axNitrogen);
title(axNitrogen,'Total nitrogen: year 0');
xlabel(axNitrogen,'x cell');
ylabel(axNitrogen,'y cell');

title(animationLayout,'Vegetation and Resource Dynamics');

drawnow;

%% ---------------- QUASI-STEADY-STATE TRACKING ----------------

% Record combined mean biomass per grid cell from year 0 onwards
mean_biomass_history = zeros(time+1,1);

mean_biomass_history(1) = ...
    sum(field_species(:))/(fieldsize^2);

% Length of each comparison window
window_length = 30;

% Allowed relative change in the window mean and standard deviation
quasi_steady_tolerance = 0.01;   % 5%

% Number of consecutive successful comparisons required
required_stable_windows = 10;

stable_window_count = 0;
quasi_steady_year = NaN;
%---------------TIME LOOP - CALCULATE CHANGE OF RESOURCE AND BIOMASS-------
for loop=1:time  
    
%Time loop - parameters recalcualted each time step
    loop 
    counter=counter+1; %counter (seperate from 'loop' value) to control number of output graphs
    %-----------MODIFY INPUT VARIABLES IF APPROPRIATE----------------------
    %Rainfall variability
    QV(1) = constant_rainfall;
    %Example of code that can be used to specify a drought (delete % signs to execute this section)
%    if ((loop>50)&&(loop<57)) %drought between years 51 and 56
%        QV(1)=25;  %i.e. very low rainfall
%    else
%        QV(1)=raindata(loop,2); %value from input data
%    end  
    %Example of code that can be used to remove some proportion of either species - Disturbances
    %if ((loop>220)&&(loop<222))  %during itteration (year)221
    %    field_species(:,:,2)=0.1; %all shrubs die leaving a seed bank
    %end
    %-----------CALL SUBROUTINES-------------------------------------------
    %PART ONE - calc randp movement (top soil layer) under action of vectors 
    [fn_r_and_p]=part_one(); %move resource according to biomass distribution
    %PART TWO - Calc use of resource by biomass, move remainder to mid and deep stores
    [fn_change_in_biomass]=part_two(); %Distribute biomass according to resource levels

    %% Record mean biomass

    mean_biomass_history(loop+1) = ...
        sum(field_species(:))/(fieldsize^2);


 
    %% ---------------- UPDATE ANIMATION ----------------

    if mod(loop,animationEvery) == 0 || loop == time

        combinedBiomass = sum(field_species,3);
        totalWater = mid_resource(:,:,1) + deep_resource(:,:,1);
        totalNitrogen = mid_resource(:,:,2) + deep_resource(:,:,2);

        % Only update existing image data
        imageBiomass.CData = combinedBiomass;
        imageWater.CData = totalWater;
        imageNitrogen.CData = totalNitrogen;

        title(axBiomass,sprintf('Combined biomass: year %d',loop));
        title(axWater,sprintf('Total water: year %d',loop));
        title(axNitrogen,sprintf('Total nitrogen: year %d',loop));

        title(animationLayout, ...
            sprintf('Vegetation and Resource Dynamics — Year %d of %d', ...
            loop,time));

        % Keep resource scales valid
        waterMin = min(totalWater(:));
        waterMax = max(totalWater(:));

        if isfinite(waterMin) && isfinite(waterMax) && waterMax > waterMin
            clim(axWater,[waterMin waterMax]);
        end

        nitrogenMin = min(totalNitrogen(:));
        nitrogenMax = max(totalNitrogen(:));

        if isfinite(nitrogenMin) && isfinite(nitrogenMax) && ...
                nitrogenMax > nitrogenMin
            clim(axNitrogen,[nitrogenMin nitrogenMax]);
        end

        drawnow;
    end
    %-----------STORE RESULTS----------------------------------------------
    time_series_plant(loop+1,1)=loop;   %Data for graphs - calculate average biomass of each species in field
    for i=1:species 
        time_series_plant(loop+1,i+1)=sum(sum(field_species(:,:,i)))/(fieldsize*fieldsize);
    end
    time_series_resource(loop+1,1)=loop; %Data for graphs - calculate average of each resource in field   
    time_series_resource(loop+1,2)=(sum(sum(mid_resource(:,:,1))))/(fieldsize*fieldsize);    
    time_series_resource(loop+1,3)=sum(sum(mid_resource(:,:,2)))/(fieldsize*fieldsize);
    time_series_resource(loop+1,4)=(sum(sum(deep_resource(:,:,1))))/(fieldsize*fieldsize); 
    time_series_resource(loop+1,5)=(sum(sum(deep_resource(:,:,2))))/(fieldsize*fieldsize); 
    con_count=zeros(fieldsize,2); %Rough connectivity calculation - applies only to species one
    for i=1:fieldsize     %connectivity count loop
        if ((field_species(i,25,1)<B_threshold(1))&&(field_species(i,25,2)<B_threshold(2)))
            con_count(i,1)=1;
        end  
        if ((i==1)&&(con_count(1,1)==1))
            con_count(i,2)=1;
        end
        if ((i>=2)&&(con_count(i,1)==1))
            con_count(i,2)=con_count(i-1,2)+con_count(i,1);
        end
    end %end of connectivity count for-loop
    [C conn]=max(con_count(:,2)); %find maximum number of continuously connected cells on transect
    max_connected_cells(loop+1,1)=loop; max_connected_cells(loop+1,2)=con_count(conn,2); %store values
    ave_connected_cells(loop+1,1)=loop; ave_connected_cells(loop+1,2)=(sum(con_count(:,1))/fieldsize)*100;
    
end %end of time loop

%% ---------------- FIND QUASI-STEADY BIOMASS ----------------

% Smooth out year-to-year fluctuations
smooth_window = 30;

smoothed_biomass = movmean(mean_biomass_history, smooth_window);

% Estimate the final quasi-steady biomass from the final 50 years
final_window = 50;

final_biomass = mean(smoothed_biomass(end-final_window+1:end));

biomass_tolerance = 0.002;   %Percentage Tolerance
required_years = 50;   %Number of time steps required to stay within the tolerance

for k = 2:(length(smoothed_biomass)-required_years+1)

    biomass_window = ...
        smoothed_biomass(k:k+required_years-1);

    relative_error = ...
        abs(biomass_window-final_biomass) ./ ...
        abs(final_biomass);

    if all(relative_error < biomass_tolerance)

       quasi_steady_year = k + required_years - 2; %Approximate year that steady state is reached
        break
    end
end

figure

plot(0:time,mean_biomass_history, 'DisplayName','Mean biomass')

hold on

plot(0:time,smoothed_biomass, 'LineWidth',2, 'DisplayName','Smoothed mean biomass')

yline(final_biomass,'--', 'DisplayName','Final quasi-steady mean biomass')

if ~isnan(quasi_steady_year)

    xline(quasi_steady_year,'--', sprintf('Quasi-steady year = %d',quasi_steady_year), 'LineWidth',1.5, 'DisplayName', 'Quasi-Steady Year');

end

xlabel('Year')
ylabel('Mean biomass per cell (g m^{-2})')
title('Mean biomass quasi-steady-state detection')

legend('Location','southeast')
grid on
end
