%Kevin M. Rivera
%
%
%
clear,clc

%fluid = Project_2(Temperature, Pressure, critP, critT, w, critDensity, [Cp polynomial coefficients], name)
tic 
CO2 = Project_2(216, 5, 73.86593, 304.25, 0.228, 10.8, [19.8 0.07344 -5.602E-5 -1.715E-8], 'CO2');
H2O = Project_2(273.15, 0.0061, 220.60, 647.0, 0.344, 17.9, [7.243E1 1.039E-2 -1.497E-06 0], 'H2O');
nC6H14 = Project_2(178.0, 1.23E-5, 30.2, 507.6, 0.305, 2.71, [-4.413 0.528 -3.119E-4 6.494E-8], 'nC6H14');
EtOH = Project_2(150.0, 4.3E-9,63, 514.0, 0.637, 6.0, [9.014 0.2141 -8.39E-5 1.373E-9], 'EtOH'); 
toc