% Wrapper code to run simulation study in reported in the paper. Outputs
% Table 1 (up to random seed).

PA = [0.5; 0.5]; % Mutation model at locus A
PB = [0.5; 0.5]; % Mutation model at locus B

R = 100; % Number of replicates
N = 1e6; % Number of steps
dt = 1e-6; % Stepsize

x0a = [0.75 0.10;
      0.05 0.10]; % Initial haplotype frequencies

x0b = [0.4 0.2;
      0.2 0.2]; % Initial haplotype frequencies

% Main experiments on varying rho
% First, tA = tB = 5
tA = 5; tB = 5;

rho = 0; % Note the relevant parameter is rho, not rho/2!
simresults = diffusionpathEuler(tA,tB,PA,PB,rho,R,N,dt,x0b);
writetable(simresults,'results_tA5e0_rho0_N1e6_x0b.txt');

rho = 0.1;
simresults = diffusionpathEuler(tA,tB,PA,PB,rho,R,N,dt,x0b);
writetable(simresults,'results_tA5e0_rho0.1_N1e6_x0b.txt');

rho = 1;
simresults = diffusionpathEuler(tA,tB,PA,PB,rho,R,N,dt,x0b);
writetable(simresults,'results_tA5e0_rho1_N1e6_x0b.txt');

rho = 2.5;
simresults = diffusionpathEuler(tA,tB,PA,PB,rho,R,N,dt,x0b);
writetable(simresults,'results_tA5e0_rho2.5_N1e6_x0b.txt');

rho = 5;
simresults = diffusionpathEuler(tA,tB,PA,PB,rho,R,N,dt,x0b);
writetable(simresults,'results_tA5e0_rho5_N1e6_x0b.txt');

rho = 10;
simresults = diffusionpathEuler(tA,tB,PA,PB,rho,R,N,dt,x0b);
writetable(simresults,'results_tA5e0_rho10_N1e6_x0b.txt');

rho = 25;
simresults = diffusionpathEuler(tA,tB,PA,PB,rho,R,N,dt,x0b);
writetable(simresults,'results_tA5e0_rho25_N1e6_x0b.txt');

chisq_critical_90 = chi2inv(0.9,1);
rhov = [0 0.1 1 2.5 5 10 25];
table_for_paper = zeros(length(rhov),9); % rho mean(rhohatMLE) bias(rhohatMLE) var(rhohatMLE) 5th-percentile median 95th-percentile Iinfs-freq Power
table_for_paper(:,1) = rhov';
for rv=1:length(rhov)
    filename = ['results_tA5e0_rho',num2str(rhov(rv)),'_N1e6_x0b.txt'];
    simresults = readtable(filename);
    rhohat = simresults.rhohat;
    rhoerror = simresults.rhoerror;
    Iinf = simresults.Iinf;
    rhohat(Iinf == 1) = rhov(rv);
    rhoerror(Iinf == 1) = 0;
    rhohatMLE = max(rhohat,0);
    rhohatMLEerror = rhohatMLE - rhov(rv);

    table_for_paper(rv,2) = mean(rhohatMLE);
    table_for_paper(rv,3) = mean(rhohatMLEerror);
    table_for_paper(rv,4) = var(rhohatMLE);
    table_for_paper(rv,5) = prctile(rhohatMLE,5);
    table_for_paper(rv,6) = median(rhohatMLE);
    table_for_paper(rv,7) = prctile(rhohatMLE,95);
    table_for_paper(rv,8) = sum(Iinf)/R;

    % Reconstruct LRT statistic
    LRT = rhohatMLE.*simresults.rhohatN;
    if (rhov(rv) == 0)
        LRT(isnan(LRT)) = 0;
    else
        LRT(isnan(LRT)) = Inf;
    end
    table_for_paper(rv,9) = sum(LRT > chisq_critical_90)/R;
end

table_for_paper
writematrix(table_for_paper,'table_tA5e0_N1e6_x0b.txt', 'Delimiter','|');

% Next tA = tB = 1
tA = 1; tB = 1;

rho = 0;
simresults = diffusionpathEuler(tA,tB,PA,PB,rho,R,N,dt,x0b);
writetable(simresults,'results_tA1e0_rho0_N1e6_x0b.txt');

rho = 0.1;
simresults = diffusionpathEuler(tA,tB,PA,PB,rho,R,N,dt,x0b);
writetable(simresults,'results_tA1e0_rho0.1_N1e6_x0b.txt');

rho = 1;
simresults = diffusionpathEuler(tA,tB,PA,PB,rho,R,N,dt,x0b);
writetable(simresults,'results_tA1e0_rho1_N1e6_x0b.txt');

rho = 2.5;
simresults = diffusionpathEuler(tA,tB,PA,PB,rho,R,N,dt,x0b);
writetable(simresults,'results_tA1e0_rho2.5_N1e6_x0b.txt');

rho = 5;
simresults = diffusionpathEuler(tA,tB,PA,PB,rho,R,N,dt,x0b);
writetable(simresults,'results_tA1e0_rho5_N1e6_x0b.txt');

rho = 10;
simresults = diffusionpathEuler(tA,tB,PA,PB,rho,R,N,dt,x0b);
writetable(simresults,'results_tA1e0_rho10_N1e6_x0b.txt');

rho = 25;
simresults = diffusionpathEuler(tA,tB,PA,PB,rho,R,N,dt,x0b);
writetable(simresults,'results_tA1e0_rho25_N1e6_x0b.txt');

rhov = [0 0.1 1 2.5 5 10 25];
table_for_paper2 = zeros(length(rhov),9); % rho mean(rhohatMLE) bias(rhohatMLE) var(rhohatMLE) 5th-percentile median 95th-percentile Iinfs-freq Power
table_for_paper2(:,1) = rhov';
for rv=1:length(rhov)
    filename = ['results_tA1e0_rho',num2str(rhov(rv)),'_N1e6_x0b.txt'];
    simresults = readtable(filename);
    rhohat = simresults.rhohat;
    rhoerror = simresults.rhoerror;
    Iinf = simresults.Iinf;
    rhohat(Iinf == 1) = rhov(rv);
    rhoerror(Iinf == 1) = 0;
    rhohatMLE = max(rhohat,0);
    rhohatMLEerror = rhohatMLE - rhov(rv);

    table_for_paper2(rv,2) = mean(rhohatMLE);
    table_for_paper2(rv,3) = mean(rhohatMLEerror);
    table_for_paper2(rv,4) = var(rhohatMLE);
    table_for_paper2(rv,5) = prctile(rhohatMLE,5);
    table_for_paper2(rv,6) = median(rhohatMLE);
    table_for_paper2(rv,7) = prctile(rhohatMLE,95);
    table_for_paper2(rv,8) = sum(Iinf)/R;

     % Reconstruct LRT statistic
    LRT = rhohatMLE.*simresults.rhohatN;
    if (rhov(rv) == 0)
        LRT(isnan(LRT)) = 0;
    else
        LRT(isnan(LRT)) = Inf;
    end
    table_for_paper2(rv,9) = sum(LRT > chisq_critical_90)/R;
end

table_for_paper2
writematrix(table_for_paper2,'table_tA1e0_N1e6_x0b.txt', 'Delimiter','|');