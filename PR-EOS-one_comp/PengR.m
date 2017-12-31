classdef Project_2 < handle
    
    %This is a class definition file that acts as a template
    %for a pure fluid being evaluated under a Peng Robinson
    %equation of state. The class handles all of the data that
    %is collected and performs various functions in order to
    %manipulate data.
    
    properties
        
        T;              %Triple point temperature of the pure fluid
        P;              %Triple point pressure of the pure fluid
        T_trip;         %Permanent holder of the triple point temperature
        P_trip;         %Permanent holder of the triple point pressure
        P_c;            %Critical pressure of the pure fluid
        T_c;            %Critical temperature of the pure fluid
        rho_c           %Critical density of the pure fluid
        omega;          %Acentric factor of the pure fluid
        a;              %a component of the PR EOS
        b;              %b component of the PR EOS
        k;              %k component of the PR EOS
        A;              %A component of the PR EOS
        B;              %B component of the PR EOS
        T_r;            %Reference temperature of the pure fluid
        P_r;            %Reference pressure of the pure fluid
        alpha;          %alpha component of the PR EOS
        Z;              %Compressibility factor matrix of the pure fluid
        Z_l;            %Liquid portion of the compressibility factor
        Z_v;            %Vapor portion of the compressibility factor
        Z_c             %Critical compressibility factor
        f_l;            %Fugacity of the liquid portion
        f_v;            %Fugacity of the vapor portion
        V;              %The vapor and liquid volumes of the fluid
        V_lNT;          %Untranslated liquid volume of the pure fluid
        V_l;            %Translated liquid volume of the pure fluid
        l_rho;          %Liquid density of the fluid
        v_rho;          %Vapor density of the fluid
        V_v;            %Vapor volume of the pure fluid
        P_sat;          %The saturated pressure of the fluid
        T_sat;          %The corresponding temperature for a saturation pressure
        Hvap;           %The enthalpy corresponding to vaporization
        H_l;            %The approximation of the real liquid enthalpy for the fluid
        H_v;            %The approximation of the real vapor enthalpy for the fluid\
        H_liquid;       %Matrix to hold all of the liquid enthalpies
        H_vapor;        %Matrix to hold all of the vapor enthalpies
        Cp;             %The heat capacity for the pure fluid
        name;           %The name of the fluid
        
    end
    
    properties (Constant)
        
        R = 0.083144;   %(L bar)/(mol K)

    end
    
    methods
        
        %The main constructor of this class used to initialize
        %some values when an object is created
        function fluid = Project_2(Temperature, Pressure, critP, critT, w, critDensity, heatCap, component)
            
            fluid.T = Temperature;        %K
            fluid.T_trip = Temperature;   %K  
            fluid.T_c = critT;            %K
            fluid.P = Pressure;           %bar
            fluid.P_trip = Pressure;      %bar
            fluid.P_c = critP;            %bar
            fluid.omega = w;              %unitless
            fluid.rho_c = critDensity;    %mol/L
            fluid.Cp = heatCap;           %J/mol
            fluid.name = component;       %string
            
            %Call to the iteration function which calls all
            %of the other essential functions
            generate_plots(fluid);

            
            
        end
        
        function fluid = get_EOS_Roots(fluid)
            
            %Multiple calls to helper functions that calculate
            %basic values that presumably do not change much
            get_a(fluid);
            get_b(fluid);
            get_k(fluid);
            get_Tr(fluid);
            get_alpha(fluid);
            get_A(fluid);
            get_B(fluid);
            imaginary_root_check(fluid);
            get_critZ(fluid);
            get_all_Vroots(fluid);

            %Finds the value of the untranslated liquid volume
            %as well as the vapor volume. The translated liquid 
            %volume is then calculated
            fluid.V_lNT = min(fluid.V);
            fluid.V_v = max(fluid.V);
            fluid.V_l = fluid.V_lNT + ((fluid.R*fluid.T_c)/fluid.P_c)*(0.3074 - fluid.Z_c);
 
        end
        
        function fluid = get_EOS_Fug(fluid)
            
            %A call to the roots function in order to update the
            %objects public properties
            fluid.get_EOS_Roots();
            
            %Retrieves the liquid and vapor compressibility factors
            get_real_compFactors(fluid);
            
            sq = sqrt(2);
            fluid.f_l = exp(fluid.Z_l-1-log(fluid.Z_l-fluid.B)-(fluid.A/(2*sq*fluid.B))*log((fluid.Z_l+(1+sq)*fluid.B)./(fluid.Z_l+(1-sq)*fluid.B)));
            fluid.f_v = exp(fluid.Z_v-1-log(fluid.Z_v-fluid.B)-(fluid.A/(2*sq*fluid.B))*log((fluid.Z_v+(1+sq)*fluid.B)./(fluid.Z_v+(1-sq)*fluid.B)));
     
        end
        
        function fluid = get_fluid_Enthalpy(fluid)
           
            fluid.P_r = fluid.P/fluid.P_c;
            get_k(fluid);
            get_alpha(fluid);
            get_A(fluid);
            get_Hvap(fluid);
            r = 8.3144;

            Cp_ = fluid.Cp(1) + fluid.Cp(2)*fluid.T + fluid.Cp(3)*fluid.T^2 + fluid.Cp(4)*fluid.T^3;
            H_ideal = Cp_*(fluid.T - fluid.T_trip);
            
            sq = sqrt(2);
            sq8 = sqrt(8);
            fluid.H_l = H_ideal + r*fluid.T*(fluid.Z_l-1 - log((fluid.Z_l+(1+sq)*fluid.B)/(fluid.Z_l+(1-sq)*fluid.B))*(fluid.A/(fluid.B*sq8))*(1+(fluid.k*sqrt(fluid.T_r))/(sqrt(fluid.alpha))));
            fluid.H_v = H_ideal + r*fluid.T*(fluid.Z_v-1 - log((fluid.Z_v+(1+sq)*fluid.B)/(fluid.Z_v+(1-sq)*fluid.B))*(fluid.A/(fluid.B*sq8))*(1+(fluid.k*sqrt(fluid.T_r))/(sqrt(fluid.alpha))))+fluid.Hvap;

        end
        
        function fluid = generate_plots(fluid)
            
            Iterate_P(fluid);
            test_plot(fluid);
            
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%Helper functions%%%%%%%%%%%%%%%%%%%%%%%%

function fluid = get_a(fluid)

    %Calculates the a component of the PR EOS
    fluid.a = (0.45724 * fluid.R^2 * fluid.T_c^2)/fluid.P_c;

end

function fluid = get_b(fluid)

    %Calculates the b component of the PR EOS
    fluid.b = (0.07780 * fluid.R * fluid.T_c) / fluid.P_c;
    
end

function fluid = get_k(fluid)

    %Calculates the k component of the PR EOS
    fluid.k = 0.37464 + 1.54226*fluid.omega - 0.26992*fluid.omega^2;

end

function fluid = get_Tr(fluid)

    %Calculates the reference temperature
    fluid.T_r = fluid.T/fluid.T_c;

end

function fluid = get_alpha(fluid)

    %Calculates the alpha component of the PR EOS
    fluid.alpha = (1 + fluid.k.*(1- fluid.T_r.^.5)).^2;

end

function fluid = get_A(fluid)

    %Calculates the A component of the PR EOS
    fluid.A = (fluid.a*fluid.alpha*fluid.P)/(fluid.R^2 * fluid.T.^2);

end

function fluid = get_B(fluid)

    %Calculates the B component of the PR EOS
    fluid.B = (fluid.b*fluid.P)./(fluid.R.*fluid.T);

end

function fluid = get_critZ(fluid)

    %Calculates the critical compressibility 
    %factor of the PR EOS for volume translation
    fluid.Z_c = fluid.P_c/(fluid.R*fluid.T_c*fluid.rho_c);
    
end

function fluid = get_all_Vroots(fluid)

    %Calculates the matrix of possible roots which 
    %includes the vapor and liquid volumes.
    fluid.V = (fluid.Z.*fluid.T.*fluid.R)./ fluid.P;

end

function fluid = imaginary_root_check(fluid)

    %Finds imaginary roots of the PR EOS and sets them to 0
    %and also reports the number of real roots the cubic
    %function contains

    %Sets up the cubic form of the PR EOS in a matrix
    cubic = [ 1 (fluid.B-1) (fluid.A - 2.*fluid.B - 3*fluid.B.^2) (-fluid.A.*fluid.B + fluid.B.^2 + fluid.B.^3)];
    Zroots = roots(cubic);
                                                %References Carl Lira
    index = find(imag(Zroots)== 0);              %code in order to 
    if length(index)>1                          %remove any imaginary
    Zreal=real(Zroots(index));                   %roots that may appear
    fluid.Z = [max(Zreal), min(Zreal)];         %in the matrix that is
    else                                        %from the root of the
    Zreal=real(Zroots(index));                   %cubic form of the PR
    fluid.Z = [max(Zreal), min(Zreal)];         %EOS. These roots cause
    end                                         %lots of issues so this 
                                                %step is neccessary
end

function fluid = get_real_compFactors(fluid)

    %Using the base formula of the compressibility factor
    %(PV/RT) in order to calculate both the liquid and
    %compressibility factors
    fluid.Z_l = (fluid.P*fluid.V_l)./(fluid.R.*fluid.T);
    fluid.Z_v = (fluid.P*fluid.V_v)./(fluid.R.*fluid.T);
    
end

function fluid = Iterate_P(fluid)

    %Error tolerance for the objective function
    error = 1E-4;
    
    %Loop counter;
    i = 1;
    
    %Infinitessimal pressure change
    dP = 1E-10;
    
    while fluid.T < fluid.T_c                       %Open outer loop
        
        %Call to the fugacity function in order to
        %get the initial liquid and vapor fugacities
        %before entering the inner loop
        get_EOS_Fug(fluid);
  
        while abs(fluid.f_l - fluid.f_v) > error    %Open inner loop
            
            %Addition of an infinitessimal pressure change
            %in order to make a derivative approximation
            fluid.P = fluid.P + dP;
            
            %A call to the fugacity method in order to update
            %the objects public properties
            get_EOS_Fug(fluid);
            fxna = (fluid.f_l - fluid.f_v).^2;
            
            %Subtraction of an the same infinitessimal pressure
            %change in order to revert the pressure back to its
            %original value
            fluid.P = fluid.P - dP;
            
            %A call to the fugacity method in order to update
            %the objects public properties
            get_EOS_Fug(fluid);
            fxnb = (fluid.f_l - fluid.f_v).^2;
            
            %Derivative approximation
            dfxn = (fxna-fxnb)./dP;             
            
            %Newton's method 
            fluid.P = fluid.P - (fxnb./ dfxn); 

        end                                         %Close inner loop
            
            get_fluid_Enthalpy(fluid);
            
            %Matrices that are generated as this 
            %function keeps looping. These are needed
            %for the graphs that are generated
            fluid.P_sat(i) = fluid.P;
            fluid.T_sat(i) = fluid.T;
            fluid.l_rho(i) = 1/fluid.V_l;
            fluid.v_rho(i) = 1/fluid.V_v;
            fluid.H_liquid(i) = fluid.H_l;
            fluid.H_vapor(i) = fluid.H_v;
            
            %An increment to the temperature in order 
            %prevent an infinite loop and also an 
            %increment to i in order to continue 
            %the dynamic expansion of the matrices
            %generated above
            fluid.T = fluid.T + 1;
            i = i + 1;
            
    end                                             %Close outer loop
              
end 

function fluid = get_Hvap(fluid)
    
    %Calculation of the heat of vaporization for a fluid
    r = 8.3144;       %J/K mol
    fluid.Hvap = -(log(fluid.P/fluid.P_r)*r)/((1/fluid.T)+(1/fluid.T_r));

end

function fluid = test_plot(fluid)
    
    %Plot for the saturation curve of a fluid
    figure;
    subplot(3,1,1);
    plot(fluid.T_sat, fluid.P_sat);
    xlabel('K');
    ylabel('Bar');
    title(strcat(fluid.name,' saturated pressure'));
    grid;
    
    %Plot for the liquid density of a fluid
    subplot(3,1,2);
    plot(fluid.T_sat,fluid.l_rho,fluid.T_sat,fluid.v_rho);
    xlabel('K');
    ylabel('mol/L');
    title(strcat(fluid.name,' saturated liquid density'));
    legend('Liquid','Vapor');
    grid;
     
    %Plot for the liquid and vapor enthalpies of a fluid
    subplot(3,1,3);
    semilogy(fluid.H_liquid, fluid.P_sat, fluid.H_vapor, fluid.P_sat);
    xlabel('J/mol');
    ylabel('Bar');
    title(strcat(fluid.name,' enthalpy of saturated liq and vap'));
    legend('Liquid','Vapor');
    grid;
    
end