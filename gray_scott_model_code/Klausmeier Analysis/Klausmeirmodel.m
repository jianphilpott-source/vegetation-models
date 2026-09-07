clear
clc
close all


%% =========================================================
%                    GRID
% ==========================================================

% x = across slope
% y = slope direction
%
% Increasing y = UPHILL
% Decreasing y = DOWNHILL

Nx = 30;
Ny = 30;
L = 20;

dx = L/Nx;
dy = L/Ny;

x = (0:Nx-1)*dx;
y = (0:Ny-1)*dy;

N = Nx*Ny;


%% =========================================================
%                  INITIAL CONDITIONS
% ==========================================================

rng(1);

w0 = 0.1 + 1.0*rand(Nx,Ny);
b0 = 0.1 + 1.0*rand(Nx,Ny);

u0 = [ ...
    w0(:);
    b0(:)];


%% =========================================================
%                     PARAMETERS
% ==========================================================

a = 2;
nu = 10;
beta = 0.540;
Dw = 1;
Db = 0.01;


%% =========================================================
%                    COLOUR MAPS
% ==========================================================

ncol = 256;


%% ---------------------------------------------------------
% Biomass colour map
% ----------------------------------------------------------

brown = [0.36 0.20 0.08];
tan   = [0.72 0.52 0.25];
light = [0.75 0.82 0.40];
green = [0.25 0.55 0.20];
dark  = [0.05 0.30 0.10];

p = [0 0.25 0.55 0.80 1];

c = [ ...
    brown;
    tan;
    light;
    green;
    dark];

biomassMap = interp1( ...
    p, ...
    c, ...
    linspace(0,1,ncol));


%% ---------------------------------------------------------
% Water colour map
% ----------------------------------------------------------

waterMap = parula(ncol);


%% =========================================================
%                  OUTPUT TIMES
% ==========================================================

% ode15s chooses its own internal time steps.
%
% These are the times stored for plotting/animation.

tEnd = 200;

dtOutput = 0.5;

tspan = 0:dtOutput:tEnd;


%% =========================================================
%                   SOLVER OPTIONS
% ==========================================================

options = odeset( ...
    'RelTol',1e-4, ...
    'AbsTol',1e-6);


%% =========================================================
%                    SOLVE MODEL
% ==========================================================

fprintf('Solving Klausmeier model...\n');

tic

[t,u] = ode15s( ...
    @klausmeir, ...
    tspan, ...
    u0, ...
    options);

elapsed = toc;

fprintf('Simulation complete.\n');
fprintf('Runtime: %.2f seconds\n\n',elapsed);


%% =========================================================
%                    FINAL STATE
% ==========================================================

uf = u(end,:)';

wf = reshape( ...
    uf(1:N), ...
    Nx,Ny);

bf = reshape( ...
    uf(N+1:2*N), ...
    Nx,Ny);


%% =========================================================
%                BIOMASS SPACE-TIME DATA
% ==========================================================

% Average biomass ACROSS the slope.
%
% Therefore the remaining spatial coordinate is y,
% which is the uphill/downhill direction.

Bprofile = zeros(length(t),Ny);

for n = 1:length(t)

    un = u(n,:)';

    bn = reshape( ...
        un(N+1:2*N), ...
        Nx,Ny);

    % Average across x

    Bprofile(n,:) = mean(bn,1);

end


%% =========================================================
%                WAVESPEED CALCULATION
% ==========================================================

% Ignore the initial pattern-forming transient.

startTime = 40;

startIndex = find( ...
    t >= startTime, ...
    1);

numberSteps = ...
    length(t)-startIndex;

shiftDistance = ...
    zeros(numberSteps,1);

waveSpeedInstant = ...
    zeros(numberSteps,1);


for n = startIndex:length(t)-1

    %% -----------------------------------------------------
    % Successive biomass profiles
    % ------------------------------------------------------

    profile1 = Bprofile(n,:);
    profile2 = Bprofile(n+1,:);


    %% -----------------------------------------------------
    % Remove mean biomass
    % ------------------------------------------------------

    profile1 = ...
        profile1 - mean(profile1);

    profile2 = ...
        profile2 - mean(profile2);


    %% -----------------------------------------------------
    % Periodic cross-correlation
    % ------------------------------------------------------

    correlation = real( ...
        ifft( ...
        fft(profile2).*conj(fft(profile1)) ...
        ));


    [~,shiftIndex] = ...
        max(correlation);


    %% -----------------------------------------------------
    % Sub-grid interpolation of correlation maximum
    %
    % This gives a much smoother estimate of wavespeed than
    % allowing only integer grid-cell shifts.
    % ------------------------------------------------------

    im1 = shiftIndex - 1;
    ip1 = shiftIndex + 1;

    if im1 < 1
        im1 = Ny;
    end

    if ip1 > Ny
        ip1 = 1;
    end


    cLeft   = correlation(im1);
    cCentre = correlation(shiftIndex);
    cRight  = correlation(ip1);


    denominator = ...
        cLeft ...
        - 2*cCentre ...
        + cRight;


    if abs(denominator) > 1e-12

        deltaPeak = ...
            0.5*(cLeft-cRight) ...
            / denominator;

    else

        deltaPeak = 0;

    end


    %% -----------------------------------------------------
    % Convert peak to periodic shift
    % ------------------------------------------------------

    shiftCells = ...
        (shiftIndex-1) ...
        + deltaPeak;


    if shiftCells > Ny/2

        shiftCells = ...
            shiftCells - Ny;

    end


    %% -----------------------------------------------------
    % Physical displacement
    % ------------------------------------------------------

    displacement = ...
        shiftCells*dy;


    dt = ...
        t(n+1)-t(n);


    index = ...
        n-startIndex+1;


    shiftDistance(index) = ...
        displacement;


    waveSpeedInstant(index) = ...
        displacement/dt;

end


%% =========================================================
%                MEAN NUMERICAL WAVE SPEED
% ==========================================================

numericalWaveSpeed = ...
    mean(waveSpeedInstant);


fprintf('==========================================\n');
fprintf('NUMERICAL BIOMASS WAVE SPEED\n');
fprintf('==========================================\n');

fprintf( ...
    'Signed wave speed = %.4f\n', ...
    numericalWaveSpeed);

fprintf( ...
    'Absolute wave speed = %.4f\n', ...
    abs(numericalWaveSpeed));

fprintf('\n');

fprintf( ...
    'Positive speed corresponds to increasing y (uphill).\n');

fprintf('==========================================\n\n');


%% =========================================================
%                     CREATE ANIMATION
% ==========================================================

fprintf('Creating animation...\n');


animationFigure = figure( ...
    'Color','w', ...
    'Position',[100 100 1000 520]);


%% ---------------------------------------------------------
% Water
% ----------------------------------------------------------

subplot(1,2,1)

imWater = imagesc( ...
    x, ...
    y, ...
    w0');

set(gca,'YDir','normal');

axis image;

colormap(gca,waterMap);

cb1 = colorbar;

cb1.Label.String = ...
    'Water';

xlabel('Across slope');

ylabel('Slope position');

title('Water');

clim([0 3]);


%% ---------------------------------------------------------
% Add slope labels
% ----------------------------------------------------------

text( ...
    -1, ...
    L*0.92, ...
    'UPHILL', ...
    'FontWeight','bold', ...
    'HorizontalAlignment','center');


text( ...
    -1, ...
    L*0.08, ...
    'DOWNHILL', ...
    'FontWeight','bold', ...
    'HorizontalAlignment','center');


%% ---------------------------------------------------------
% Biomass
% ----------------------------------------------------------

subplot(1,2,2)

imBiomass = imagesc( ...
    x, ...
    y, ...
    b0');

set(gca,'YDir','normal');

axis image;

colormap(gca,biomassMap);

cb2 = colorbar;

cb2.Label.String = ...
    'Biomass';

xlabel('Across slope');

ylabel('Slope position');

title('Biomass');

clim([0 4]);


sgtitle( ...
    'Klausmeier model: t = 0');


drawnow;


%% =========================================================
%                     VIDEO WRITER
% ==========================================================

video = VideoWriter( ...
    'klausmeier_uphill_bands.mp4', ...
    'MPEG-4');


video.FrameRate = 20;

video.Quality = 100;

open(video);


%% =========================================================
%                  PLAY AND SAVE FRAMES
% ==========================================================

% 1 = every stored state
% 2 = every second stored state
% etc.

frameSkip = 1;


for n = 1:frameSkip:length(t)

    %% -----------------------------------------------------
    % Retrieve current solution
    % ------------------------------------------------------

    un = u(n,:)';


    %% -----------------------------------------------------
    % Water
    % ------------------------------------------------------

    wn = reshape( ...
        un(1:N), ...
        Nx,Ny);


    %% -----------------------------------------------------
    % Biomass
    % ------------------------------------------------------

    bn = reshape( ...
        un(N+1:2*N), ...
        Nx,Ny);


    %% -----------------------------------------------------
    % Update images
    % ------------------------------------------------------

    set( ...
        imWater, ...
        'CData', ...
        wn');


    set( ...
        imBiomass, ...
        'CData', ...
        bn');


    %% -----------------------------------------------------
    % Update title
    % ------------------------------------------------------

    sgtitle( ...
        sprintf( ...
        'Klausmeier model: t = %.1f', ...
        t(n)));


    %% -----------------------------------------------------
    % Draw frame
    % ------------------------------------------------------

    drawnow;


    %% -----------------------------------------------------
    % Capture frame
    % ------------------------------------------------------

    frame = ...
        getframe(animationFigure);


    %% -----------------------------------------------------
    % Write frame to MP4
    % ------------------------------------------------------

    writeVideo( ...
        video, ...
        frame);

end


%% =========================================================
%                    CLOSE VIDEO
% ==========================================================

close(video);


fprintf('Animation complete.\n');
fprintf('Video saved as:\n');
fprintf('klausmeier_uphill_bands.mp4\n\n');


%% =========================================================
%                   WAVE SPEED PLOT
% ==========================================================

speedTimes = ...
    t(startIndex:end-1);


figure('Color','w');


plot( ...
    speedTimes, ...
    waveSpeedInstant, ...
    'LineWidth',1.5);


hold on


yline( ...
    numericalWaveSpeed, ...
    '--', ...
    sprintf( ...
    'Mean c = %.3f', ...
    numericalWaveSpeed), ...
    'LineWidth',1.5);


xlabel('Time');

ylabel('Wave speed');


title( ...
    'Numerically measured biomass wave speed');


grid on;


%% =========================================================
%                   SPACE-TIME PLOT
% ==========================================================

figure( ...
    'Color','w', ...
    'Position',[200 100 700 600]);


imagesc( ...
    y, ...
    t, ...
    Bprofile);


set(gca,'YDir','normal');


xlabel( ...
    'Slope position');


ylabel( ...
    'Time');


title( ...
    'Biomass space-time diagram');


cb = colorbar;

cb.Label.String = ...
    'Mean biomass';


colormap(biomassMap);


%% =========================================================
%                     FINAL PLOT
% ==========================================================

figure( ...
    'Color','w', ...
    'Position',[100 100 950 500]);


%% ---------------------------------------------------------
% Final water
% ----------------------------------------------------------

subplot(1,2,1)


imagesc( ...
    x, ...
    y, ...
    wf');


set(gca,'YDir','normal');

axis image;


colormap(gca,waterMap);


cb1 = colorbar;

cb1.Label.String = ...
    'Water';


xlabel('Across slope');

ylabel('Slope position');


title( ...
    sprintf( ...
    'Final water, t = %.1f', ...
    t(end)));


clim([0 3]);


%% ---------------------------------------------------------
% Final biomass
% ----------------------------------------------------------

subplot(1,2,2)


imagesc( ...
    x, ...
    y, ...
    bf');


set(gca,'YDir','normal');

axis image;


colormap(gca,biomassMap);


cb2 = colorbar;

cb2.Label.String = ...
    'Biomass';


xlabel('Across slope');

ylabel('Slope position');


title( ...
    sprintf( ...
    'Final biomass, t = %.1f', ...
    t(end)));


clim([0 4]);


%% =========================================================
%                  KLAUSMEIER MODEL
% ==========================================================

function dudt = klausmeir(~,u)


%% =========================================================
%                     PARAMETERS
% ==========================================================

a = 2;
nu = 10;
beta = 0.540;
Dw = 1;
Db = 0.01;


%% =========================================================
%                         GRID
% ==========================================================

Nx = 30;
Ny = 30;

L = 20;

dx = L/Nx;
dy = L/Ny;

N = Nx*Ny;


%% =========================================================
%                   RESHAPE VARIABLES
% ==========================================================

w = reshape( ...
    u(1:N), ...
    Nx,Ny);


b = reshape( ...
    u(N+1:2*N), ...
    Nx,Ny);


%% =========================================================
%              PERIODIC X-BOUNDARIES
%
%                  ACROSS SLOPE
% ==========================================================

wRight = ...
    circshift(w,[-1 0]);


wLeft = ...
    circshift(w,[1 0]);


bRight = ...
    circshift(b,[-1 0]);


bLeft = ...
    circshift(b,[1 0]);


%% =========================================================
%              PERIODIC Y-BOUNDARIES
%
%                  SLOPE DIRECTION
% ==========================================================

wUp = ...
    circshift(w,[0 -1]);


wDown = ...
    circshift(w,[0 1]);


bUp = ...
    circshift(b,[0 -1]);


bDown = ...
    circshift(b,[0 1]);


%% =========================================================
%                    X LAPLACIAN
% ==========================================================

lapw_x = ...
    ( ...
    wRight ...
    - 2*w ...
    + wLeft ...
    )/dx^2;


lapb_x = ...
    ( ...
    bRight ...
    - 2*b ...
    + bLeft ...
    )/dx^2;


%% =========================================================
%                    Y LAPLACIAN
% ==========================================================

lapw_y = ...
    ( ...
    wUp ...
    - 2*w ...
    + wDown ...
    )/dy^2;


lapb_y = ...
    ( ...
    bUp ...
    - 2*b ...
    + bDown ...
    )/dy^2;


%% =========================================================
%                  TOTAL LAPLACIAN
% ==========================================================

lapw = ...
    lapw_x + lapw_y;


lapb = ...
    lapb_x + lapb_y;


%% =========================================================
%                 WATER GRADIENT ALONG SLOPE
% ==========================================================

% y increases UPHILL.
%
% Therefore downhill water movement is toward decreasing y.
%
% For
%
%     w_t = + nu w_y
%
% the corresponding transport velocity is -nu,
% i.e. toward decreasing y.

gradwy = ...
    (wUp-wDown)/(2*dy);


%% =========================================================
%                     WATER PDE
% ==========================================================

dw = ...
    a ...
    - w ...
    - w.*b.^2 ...
    + nu*gradwy ...
    + Dw*lapw;


%% =========================================================
%                    BIOMASS PDE
% ==========================================================

db = ...
    w.*b.^2 ...
    - beta*b ...
    + Db*lapb;


%% =========================================================
%                     RETURN
% ==========================================================

dudt = [ ...
    dw(:);
    db(:)];

end