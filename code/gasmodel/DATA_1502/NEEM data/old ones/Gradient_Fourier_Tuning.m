
% -- Tuning parameters
Tuning_step = 0.01;
Iterations = 220;
Nrfunctions = 31;

% -- Initial guess for tortuosity profile

Calc_profiles;

for teller_all=1:Iterations

    Misfit = Calc_misfit(z,c_profiles,gas_samples,gases);


    Cost = 0;
    for teller = 1: length(gases)
        Cost = Cost + sqrt(mean((Misfit{teller}(:,2)).^2));
    end

    Coefficients = zeros(1,Nrfunctions);
    Cost_coefficients = zeros(1,Nrfunctions);
    for teller_GFT = 1:Nrfunctions;
        Coefficients =  zeros(1,Nrfunctions);
        Coefficients(teller_GFT) = 1;

        InvTort_temp = min(1,InvTort.* (1+Tuning_step*Tort_modifiers(z,Lockin,Coefficients)));
        InvTort_temp(z>Lockin)=0;
        InvTort_temp = efilter(InvTort_temp).*(z<74);

        for teller = 1:length(gases)
            gas = gases{teller};
            [alpha,beta,gamma] = Crank_abc(N,z,dz,dt,Acc,rho,rho_ice,s_open,s_closed,M_air,M_gas(strcmp(gases,gas)),D_gas(strcmp(gases,gas)),decay_gas(strcmp(gases,gas)),g,R,T,D_0eddy, H_eddy,teller_co,InvTort_temp, Diff_m);
            c_gases = CrankNic(alpha,beta,gamma,M,N,P,P2,teller_co,gas_history(strcmp(gases,gas),:),N_over);
            c_profiles(teller,:) = c_gases(:,end);
        end
        % display('  Profiles calculated')


        if (sum(strcmp(gases,'C14_CO2'))|sum(strcmp(gases,'C13_CO2')))
            if sum(strcmp(gases,'C14_CO2')) == 0
                C13_profile = 0;
            end
            c_profiles = Conc_to_delta(c_profiles,gases,C13_profile);
        end

        Misfit = Calc_misfit(z,c_profiles,gas_samples,gases);

        for teller = 1: length(gases)
            Cost_coefficients(teller_GFT) = Cost_coefficients(teller_GFT) + sqrt(mean((Misfit{teller}(:,2)).^2));
        end

    end
    Cost_coefficients = (Cost-Cost_coefficients);
    Cost_coefficients = (Cost_coefficients)./(max(max(abs(Cost_coefficients))));
    
    InvTort = min(1,InvTort.*(1+Tuning_step*Tort_modifiers(z,Lockin,Cost_coefficients)));
    InvTort(z>Lockin)=0;
    InvTort = efilter(InvTort).*(z<74);

    for teller = 1:length(gases)
        gas = gases{teller};
        [alpha,beta,gamma] = Crank_abc(N,z,dz,dt,Acc,rho,rho_ice,s_open,s_closed,M_air,M_gas(strcmp(gases,gas)),D_gas(strcmp(gases,gas)),decay_gas(strcmp(gases,gas)),g,R,T,D_0eddy, H_eddy,teller_co,InvTort_temp, Diff_m);
        c_gases = CrankNic(alpha,beta,gamma,M,N,P,P2,teller_co,gas_history(strcmp(gases,gas),:),N_over);
        c_profiles(teller,:) = c_gases(:,end);
    end
    % display('  Profiles calculated')


    if (sum(strcmp(gases,'C14_CO2'))|sum(strcmp(gases,'C13_CO2')))
        if sum(strcmp(gases,'C14_CO2')) == 0
            C13_profile = 0;
        end
        c_profiles = Conc_to_delta(c_profiles,gases,C13_profile);
    end

end



