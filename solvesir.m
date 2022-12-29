% Ryker Dolese, CAAM 210, FALL 2022, Project 5 on Infectious Disease Modeling
% solvesir.m
% This script is for the 6th Project on Canvas (Pledged)
% Last Modified: October 12th, 2022


function solvesir
M=7900000; % set our initial population
I0=10; % initial number of infections
initialval = [M-I0, 0]; % suceptibles is population - infected and recovered is 0
simpleSIR (M, 0.7, 0.1, initialval, 150); % call simple SIR function
mu = 1./(76*365); % define mu for variableSIR
variableSIR (0.5, 0.33, 0.01, mu, 7900000, 0, 10, 150); % call variable SIR
omega = 1./365; % define omega for variableimmSIR
variableimmSIR (0.5, 0.33, 0.01, mu, omega, 7900000, 0, 10, 1460);
% call variableimmSIR function, which accounts for loss of immunity

end

function [S, R, I] = simpleSIR (M, alpha, beta, initialval, Tfinal)
% Inputs:
% - M (total population)
% - alpha (fixed numbercontacts per day, sufficient to spread the disease)
% - beta (fraction of infected group that will recover per day)
% - initialval (vector of initial # of susceptibles and recovered)
% - Tfinal (total time that we track populations)
% Outputs:
% - S (# of susceptibles)
% - R (# of recovered individuals)
% - I (# of infected individuals)
% This function uses the SIR model to predict how an infectious disease
% spreads throughout a population given a set number of inputs.
t = 1:Tfinal; % time vector
S(1) = initialval(1); % intitial # of susceptibles
R(1) = initialval(2); % intitial # of Recovered
I(1) = M-S(1)-R(1); % intitial # of Infected

for i = 1:length(t)-1 % run simulation for given amount of time
    S(i+1)=S(i)-(alpha/M)*S(i)*I(i);
    R(i+1)=R(i)+beta*I(i);
    I(i+1)=M-S(i+1)-R(i+1);
    % calculating our different populations for each time iteration
end

figure % produce figure
hold on
plot(t, S, 'b', 'linewidth', 2) % plot susceptibles
plot(t, R, 'r', 'linewidth', 2) % plot recovered
plot(t, I, 'y', 'linewidth', 2) % plot infected
legend ('Susceptibles', 'Recovered', 'Infected') % produce legend
xlabel('# of Days')
ylabel('Population Sizes')
% label our axes accordingly
title('Evolution of the population groups over 150 days') % give our plot a great title

% Question 1 - part 2
% Since we increased alpha (the number of contacts an infected individual
% has per day), we see that the infections increase exponentially much
% sooner. Additionally, since beta (the recovery rate) decreased, we can
% see that the amount of recovered people tends to be delayed until
% after the peak of infections. Due to these factors, nearly everyone is
% infected and the susceptible count approaches zero.

end

function variableSIR (alpha,beta,gamma,mu, S0, R0, I0, Tfinal)
% Inputs:
% - alpha (fixed numbercontacts per day, sufficient to spread the disease)
% - beta (fraction of infected group that will recover per day)
% - gamma (loss of individuals from the infected group)
% - mu (deaths from other causes)
% - S0 (initial # of susceptibles)
% - R (initial # of recovered individuals)
% - I (initial # of infected individuals)
% - Tfinal (total time that we track populations)
% Outputs:
% This function uses the SIR model to predict how an infectious disease
% spreads throughout a population given a set number of inputs. This model
% takes into consideration the fluctuation of total population and deaths
% associated with the infection.
t = 1:Tfinal; % time vector
S(1) = S0; % assign initial # of susceptibles
R(1) = R0; % assign initial # of recovered individuals)
I(1) = I0; % assign initial # of infected individuals
M(1) = R0+S0+I0; % total population = susceptibles + recovered + infected
for i = 1:length(t)-1
    S(i+1)=S(i)-(alpha/M(i))*S(i)*I(i)+(mu*M(i))-(mu*S(i));
    R(i+1)=R(i)+beta*I(i)-mu*R(i);
    I(i+1)=I(i)+(alpha/M(i))*S(i)*I(i)-beta*I(i)-(mu+gamma)*I(i);
    M(i+1)= M(i)-gamma*I(i);

    % calculating our different populations for each time iteration; this
    % time we include total population, since it fluctuates due to deaths
end

figure % produce figure
hold on
plot(t, S, 'b', 'linewidth', 2) % plot susceptibles
plot(t, R, 'r', 'linewidth', 2) % plot recovered
plot(t, I, 'y', 'linewidth', 2) % plot infected
plot(t, M, 'c', 'linewidth', 2) % plot of total population
legend ('Susceptibles', 'Recovered', 'Infected', 'Total Population') % produce legend
xlabel('# of Days')
ylabel('Population Sizes')
% label our axes accordingly
title('Evolution of the population groups over 150 days') % give our plot a great title
hold off

figure % produce 2nd figure
plot(S, I, 'linewidth', 3) % plot susceptibles (x axis) and infections (y axis)
xlabel('# of Susceptibles')
ylabel('# of Infected')
% label the axes
title('Infections as a function of Susceptibles') % give it a title
hold off

end

function variableimmSIR (alpha,beta,gamma,mu, omega, S0, R0, I0,Tfinal)
% Inputs:
% - alpha (fixed numbercontacts per day, sufficient to spread the disease)
% - beta (fraction of infected group that will recover per day)
% - gamma (loss of individuals from the infected group)
% - mu (deaths from other causes)
% - omega (fraction of recovered population that become susceptible again)
% - S0 (initial # of susceptibles)
% - R (initial # of recovered individuals)
% - I (initial # of infected individuals)
% - Tfinal (total time that we track populations)
% Outputs:
% This function uses the SIR model to predict how an infectious disease
% spreads throughout a population given a set number of inputs. This model
% takes into consideration the fluctuation of total population and deaths
% associated with the infection. Additionally, this model accounts for loss
% of immunity from the recovered population.
t = 1:Tfinal; % time vector
S(1) = S0; % assign initial # of susceptibles
R(1) = R0; % assign initial # of recovered individuals)
I(1) = I0; % assign initial # of infected individuals
M(1) = S0+R0+I0; % total population = susceptibles + recovered + infected
for i = 1:length(t)-1
    S(i+1)=S(i)-(alpha/M(i))*S(i)*I(i)+(mu*M(i))-(mu*S(i))+omega*R(i);
    R(i+1)=R(i)+beta*I(i)-mu*R(i)-omega*R(i);
    I(i+1)=I(i)+(alpha/M(i))*S(i)*I(i)-beta*I(i)-(mu+gamma)*I(i);
    M(i+1)= M(i)-gamma*I(i);
    % calculating our different populations for each time iteration; this
    % time we subtract omega*R(i) from the # of recovered and add this to
    % the # of susceptibles. This takes into consideration the loss of
    % immunity.
end
figure % produce figure
hold on
plot(t, S, 'b', 'linewidth', 2) % plot susceptibles
plot(t, R, 'r', 'linewidth', 2) % plot recovered
plot(t, I, 'y', 'linewidth', 2) % plot infected
plot(t, M, 'c', 'linewidth', 2) % plot of total population
legend ('Susceptibles', 'Recovered', 'Infected', 'Total Population')
% make a legend
xlabel('# of Days')
ylabel('Population Sizes')
% label our axes
title('Evolution of the population groups over 4 years')
% make a new title (this plot is over 4 years)
hold off

figure % produce new figure
plot(S, I, 'linewidth', 3) % plot infections v. susceptibles
xlabel('# of Susceptibles')
ylabel('# of Infected')
% label the axes
title('Phase Plane of Epidemic Wave') % title our plot accordingly
hold off
end



