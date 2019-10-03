clear all
clc
%%%%%%%%%%%%%%%%%%%%%%% Stochastic Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Authors: Ding Yuan (dy248) and Ian Paul (ijp9)
%
%   Based on our research, there are approximately 60k people living in
%   Ithaca. Hence create a map (a matrix A) with dim(A) = 244 by 244, each 
%   grid represent one person so that total number of people would be
%   roughly 60k. Since based on our research, there is no H3N2 case in
%   Ithaca right now, so we will simulate outbreak as follow: one person in
%   Cornell region is infected and start to spread out the disease. (We 
%   choose starting at Cornell region because it make more sense since 
%   students or professors often travel outside Ithaca. When they get back
%   to campus, they probably bring disease into Ithaca.)
%
%   We divide total population into four groups:
%
%               S: Susceptible population (denote as state 1)
%               I: Infected population    (denote as state 2)
%               R: Recovered population   (denote as state 3)
%               D: Dead population        (denote as state 4)
%
%   Assumption:
%               1. 30 days/month, 360 days/year
%               2. People get vaccinated at the
%                  begninning of a month
%               3. The number of days in flu season in 150
%               4. 


str_x = 0;
str_y = 1;
num_runs = 500;
dead_std_devs = double.empty(length(str_x), 0);
dead_means = double.empty(length(str_x), 0);
infected_std_devs = double.empty(length(str_x), 0);
infected_means = double.empty(length(str_x), 0);

count = 1;
while str_x <= 1
    
    num_infected_runs = double.empty(num_runs, 0);
    num_dead_runs = double.empty(num_runs, 0);
    
    str_x = str_x + 0.1;
    str_y = 1 - str_x;
   
    for i = 1:1:num_runs
        [num_infected_runs(i), num_dead_runs(i)] = simulation(str_x, str_y, 120, 100, 150);
    end
    
    infected_means(count) = mean(num_infected_runs) % mean number of infected people over num_runs with a str_x % vax alloaction to students and 1-str_x % vax allocation to elderly
    infected_std_devs(count) = std(num_infected_runs) % std dev. of num infected people over num_runs under str_x%, 1-str_x% vax allocation to students,elders respectively
    
    dead_means(count) = mean(num_dead_runs);
    dead_std_devs(count) = std(num_dead_runs);
    count = count + 1;
end

figure(1)
scatter(0:10:10*length(infected_means)-10, infected_means)
title('Number of Infected People vs. % Vax Given to Small Children & Students / Month')
ylabel('Number of Infected People Over Flu Season')
xlabel('% Vax Given to Small Children & Students / Month')
grid

figure(2)
scatter(0:10:10*length(dead_means)-10, dead_means)
title('Number of Dead People vs. % Vax Given to Small Children & Students / Month')
ylabel('Number of Dead People Over Flu Season')
xlabel('% Vax Given to Small Children & Students / Month')
grid

figure(3)
scatter(0:10:10*length(infected_std_devs)-10, infected_std_devs)
title('Std. Dev. of Number of Infected People vs. % Vax Given to Small Children & Students / Month')
ylabel('Std. Dev. of Number of Infected People Over Flu Season')
xlabel('% Vax Given to Small Children & Students / Month')
grid

figure(4)
scatter(0:10:10*length(dead_std_devs)-10, dead_std_devs)
title('Standard Deviation of Number of Dead People vs. % Vax Given to Small Children & Students / Month')
ylabel('Standard Deviation of Number of Dead People Over Flu Season')
xlabel('% Vax Given to Small Children & Students / Month')
grid

function [infected, dead] = simulation(str_x, str_y, patient_zero_x, patient_zero_y, flu_szn_days )
    str = [str_x, str_y];  % Strategy vector: give students str_x %, elder str_y %. 

    %%% Initial state matrix for everyone in Ithaca
    A = ones(244,244);

    A(patient_zero_x, patient_zero_y) = 2;    % the person who bring the disease into Ithaca
    %A(200,150) = 2;
    %A(1,100) = 2;

    % Status matrix for vaccine
    %             B(i,j) = 1 if vaccinated
    %             B(i,j) = 0 if not vaccinated
    B = zeros(244,244);  

    % Based on research, 21.5% of total ppl 
    num_elders = 0.215 * 244^2;

    % Infect rate for students & general population
    ifrstu = 0.3;
    ifreld = 0.1;

    % death rate
    mryoung = 0.00001;
    mrstu = 0.00002;
    mrgeneral = 0.0002;
    mrold = 0.00135;

    % Counter for # of dead
    num_dead = 0;

    %%% Set up infection rate matrix 
    IFR = zeros(244,244);   
    for i = 1:1:244
        for j = 1:1:244
            if (i >= 110) && (j >= 90)
                IFR(i,j) = ifrstu;
            elseif (i>=60) && (i<110) && (j>=90) && (j<= 140)
                IFR(i,j) = ifrstu;
            elseif (i>=60) && (i<=110) && (j>=180) && (j<=220)
                IFR(i,j) = ifrstu;
            else
                IFR(i,j) = ifreld;
            end
        end
    end


    p1 = makepoint(patient_zero_x, patient_zero_y, randi([7 10000]));
    InfecPpl = [p1];      % initialize a structure array for infected ppl (list of anyone infected at any time during flu season)

    %%% Start Simulation

    for t = 1 : 1 : flu_szn_days
        if rem(t,30) == 1    % distribute 3800 vaccines, assuming 30 days / month
            % distribute to students
            for v =1:1:floor(3800*str(1))
                x = 109+ceil(135*rand);
                y = 89+ceil(155*rand);
                while B(x,y)==1 
                    x = 109+ceil(135*rand);
                    y = ceil(244*rand);
                end
                B(x,y) = 1;
            end

            % distribute to elders and general
            for vv = 1:1:floor(3800*str(2))
                xe = ceil(109*rand);
                ye = ceil(244*rand);
                while B(xe,ye) ==1
                    xe = ceil(109*rand);
                    ye = ceil(109*rand);
                end
                B(xe,ye) = 1;
            end
        end

        n = length(InfecPpl);
        counter_n = n;
        for e = 1:1:n
            
            cur_x = InfecPpl(e).x;
            cur_y = InfecPpl(e).y;
            
            InfecPpl(e).t = InfecPpl(e).t - 1;
            if InfecPpl(e).t <= 0       % amount of time needed to recovery passed
                A(cur_x, cur_y) = 3;    % person is no longer in infected state and in recovered state     
                continue
            end
            
            direc = rand;
            delx = 0;
            dely = 0;
            if direc <= 0.25            % up
                dely = 1;
            elseif direc <= 0.5          % down
                dely = -1;
            elseif direc <= 0.75        % left
                delx = -1;
            else                        % right
                delx = 1;
            end
            new_x = cur_x + delx;
            new_y = cur_y + dely;
            new_t = InfecPpl(e).t - 1;
            
            if A(new_x, new_y) == 1     % if susceptible
                check = rand;
                if new_x >= 110 && new_y >= 90  % in school region
                    if (check <= ifrstu) && (B(new_x,new_y)~=1)
                        A(new_x, new_y) = 2; % change person's state to infected
              
                        InfecPpl(counter_n+1) = makepoint(new_x,new_y,new_t);
                        if rand <= mrstu
                            num_dead = num_dead + 1;
                            A(new_x, new_y) = 4; % change person's state to dead
                        end
                    end
                elseif new_x>= 60 && new_y >=90 && new_y <= 140 && (B(new_x,new_y)~=1) % in downtown
                    if check <= ifrstu
                        A(new_x, new_y) = 2; % change person's state to infected
                        InfecPpl(counter_n+1) = makepoint(new_x,new_y,new_t);
                        if rand <= mrgeneral
                            num_dead = num_dead + 1;
                            A(new_x, new_y) = 4; % change person's state to dead
                        end
                    end
                elseif new_x >= 60 && new_y >=180 && new_y <= 220 && (B(new_x,new_y)~=1)
                    if check <= ifrstu
                        A(new_x, new_y) = 2; % change person's state to infected
                        InfecPpl(counter_n+1) = makepoint(new_x,new_y,new_t);
                        if rand <= mrstu
                            num_dead = num_dead + 1;
                            A(new_x, new_y) = 4; % change person's state to dead
                        end                 
                    end
                else
                    if check <= ifreld
                        A(new_x, new_y) = 2; % change person's state to infected
                        InfecPpl(counter_n+1) = makepoint(new_x,new_y,new_t);
                        if rand <= (0.55*mrgeneral + 0.45*mrold)
                            num_dead = num_dead + 1;
                            A(new_x, new_y) = 4; % change person's state to dead
                        end    
                    end
                end        
            end
        end
    end
    
infected = length(InfecPpl);
dead = num_dead;

end







