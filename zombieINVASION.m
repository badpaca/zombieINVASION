function [nalive, V, nimmu, iter, mini, minpercent] = zombieINVASION

percentvect = linspace(0,0.98,100); %the percents iterated over
beforecostvect = zeros(1,100); %will fill this in thru the loop (it's the 
                          %number of people cured in the end for each
beforecostvect2 = zeros(1,100); %will fill this in thru the loop (it's the
                            %cost of vaccinating the population chosen
for pvaci = 1:100 %THE LOOP!!!!
    tic 
    %% -----------------
    %% Definitions
    %% -----------------

    N = 128;        % grid size: N*N
    T = 15;         % hours remaining before immunized/recovery
    Nitermax = 5000;    % max iterations

    density = 0.5;  % density of population
    prob = 0.16;     % probability of infection
    pvac = percentvect(pvaci);     % fraction of walkers given vaccine 
                                    % (loop over this)

    M = floor(density*N^2);     % number of walkers
    beforecostvect2(pvaci) = M*pvac*30; %cost of vaccinating
    V = floor(pvac*M);

    plot_flag = 1;

    %% -----------------
    %% Initialization
    %% -----------------

    infect = zeros(M,1);    
    % infect takes the values: 0 (healthy); 1 (zombie); 2 (immunized); 3 (treated)
    hours = T*ones(M,1);     % hours gets initialized to T upon infection and 
    % decreases by 1 every iteration until hours==0, which implies immunization

    % Initialize the walkers' starting positions randomly
    x = randi(N,M,1);
    y = randi(N,M,1);

    % Vaccination campaign for V walkers
    infect(randperm(M,V)) = 2;

    % Infect a random walker
    z = randi(M); infect(z) = 1;

    nimmu = 50;    % number of people recovered initially
    nzombies = 1;   % first infected walker
    iter = 1;       % about to start first iteration


    %% -----------------
    %% Random walk
    %% -----------------

    while nzombies > 0 && iter < Nitermax

        % Count infection countdowns and deaths before taking a random step
        hours(infect==1) = hours(infect==1) - 1;

        % Does the zombie get immunized?
        infect(hours==0) = 3;     
        nimmu = nimmu + sum(infect==3);

        % Random walk

        rand_walk = rand(M,1);
        walk_x = rand_walk < 0.5;
        walk_y = rand_walk >= 0.5;

            % Walk in the x direction
            deltax = 2*floor(2*rand(size(x))) - 1;  % random directions
            deltax(infect==3) = 0;  % the dead don't move
            deltax(walk_y) = 0;     % don't walk in y!
            xold = x;
            xnew = x + deltax; 
            % Periodic boundary conditions
            xnew(xnew==0) = N;
            xnew(xnew==N+1) = 1;
            x = xnew;

            % Walk in the y direction
            deltay = 2*floor(2*rand(size(x))) - 1;  % random directions
            deltay(infect==3) = 0;  % the dead don't move
            deltay(walk_x) = 0;     % don't walk in x!
            yold = y;
            ynew = y + deltay;
            % Periodic boundary conditions
            ynew(ynew==0) = N;
            ynew(ynew==N+1) = 1;
            y = ynew;

        % Contamination step 
        H_coord = zeros(M,2);     % because we need to keep track of positions
        pos_H = infect==0;
        H_coord(pos_H,:) = [x(pos_H), y(pos_H)];  % positions of healthy people

        pos_I = infect==1;
        I_coord = [x(pos_I), y(pos_I); xold(pos_I), yold(pos_I)];  
        % positions of zombies, before and after random step

        indH = find( ismember(H_coord,I_coord,'rows') == 1);
        % indH are the "j-positions" of the healty walkers who are
        % getting infected (infect(j) = 0 -> 1 for j in indH)

        infection_status = (rand(size(indH)) < prob); % true = turns into zombie
        infect(indH) = infection_status;
        nzombies = sum(infect==1);

        % Increment
        iter = iter + 1;
        if iter == Nitermax
            display('Maximum number of iterations reached; increase Nitermax')
        end

        % Visualization in 2D
        if plot_flag == 1
            colorgrid = 4*ones(N);
            map = [ 0, 1, 0;    % healthy = green
                    1, 0, 0;    % zombies = red
                    0, 0, 1;    % immunized = blue
                    0, 0, 0;    % dead = black
                    1, 1, 1];   % empty = white
            for i=1:M
                colorgrid(x(i),y(i)) = infect(i);
            end
            g = pcolor(colorgrid); colormap(map); axis square; 
            set(g,'LineStyle','none'); 
            set(gca,'XDir','normal')
            drawnow
        end

    end

    nalive = M - nzombies - nimmu - V; % alive and not immunized
    running_time = toc;
    beforecostvect(pvaci) = nimmu; 
end

%disp(horzcat('number of people cured: ',num2str(nimmu)));

            %PLOTTING STUFF (optional)

% costvect = beforecostvect.*100 + beforecostvect2; %The total amount of $$
% plot(percentvect.*100, costvect, 'LineWidth',2);
% ylim([0,200000]);
% xlabel('Percent Vaccinated'); ylabel('Total Cost in $'); 
% title('Plot for 0.5 Density');
% [mini,index] = min(costvect); %To find the index of the min total cost
% minpercent = percentvect(index)*100; %The percent at which $ is a min
% disp('Cheapest total cost');
% disp(num2str(mini));
% disp('Ideal % of people vaccinated');
% disp(num2str(minpercent));
end