%Kevin M. Rivera
%Parity plots for comparing the peng robinson theoretical
%data to the experimental data found on NIST.
%
%
%
%%%%%%%%%%%%%%%%%Comparing saturation pressures%%%%%%%%%%%%%%%%%%%
figure;
plot(CO2_PR_P, CO2_NIST_P, '.');
line([0 80],[0 80]);
xlabel('PR pressure (bar)');
ylabel('NIST pressure (bar)');
title('Saturation pressure parity plot for CO2');
grid;

figure;
plot(H2O_PR_P, H2O_NIST_P, '.');
line([0 250],[0 250]);
xlabel('PR pressure (bar)');
ylabel('NIST pressure (bar)');
title('Saturation pressure parity plot for H2O'); 
grid;

figure;
plot(nC6H14_PR_P, nC6H14_NIST_P, '.');
line([0 30],[0 35]);
xlabel('PR pressure (bar)');
ylabel('NIST pressure (bar)');
title('Saturation pressure parity plot for nC6H14');
grid;

figure;
plot(EtOH_PR_P, EtOH_NIST_P, '.');
line([0 70],[0 90]);
xlabel('PR pressure (bar)');
ylabel('NIST pressure (bar)');
title('Saturation pressure parity plot for EtOH');
grid;

%%%%%%%%%%%%%%%%%Comparing saturated liquid densities%%%%%%%%%%%%%
figure;
plot(CO2_PR_rho, CO2_NIST_rho,'.');
line([8 22],[10 28]);
xlabel('PR \rho (mol/L)');
ylabel('NIST \rho (mol/L)');
title('Liquid density parity plot for CO2');
grid;

figure;
plot(H2O_PR_rho, H2O_NIST_rho,'.');
line([10 24],[20 55]);
xlabel('PR \rho (mol/L)');
ylabel('NIST \rho (mol/L)');
title('Liquid density parity plot for H2O');
grid;

figure;
plot(nC6H14_PR_rho, nC6H14_NIST_rho,'.');
line([2 5.5],[3 7.5]);
xlabel('PR \rho (mol/L)');
ylabel('NIST \rho (mol/L)');
title('Liquid density parity plot for nC6H14');
grid;

figure;
plot(EtOH_PR_rho, EtOH_NIST_rho,'.');
line([4.5 9.5],[10 24]);
xlabel('PR \rho (mol/L)');
ylabel('NIST \rho (mol/L)');
title('Liquid density parity plot for EtOH');
grid;