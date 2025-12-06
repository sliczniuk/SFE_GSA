function D = Compute_Diffusion_for_GSA(X)

    T = X(:,1); P = X(:,2); F = X(:,3);

    Parameters_table = readtable('Parameters.csv') ;                     % Table with prameters
    parameters       = num2cell(Parameters_table{:,3});                  % Parameters within the model + (m_max), m_ratio, sigma

    Z                = Compressibility(T, P, parameters);
    RHO              = rhoPB_Comp(T, P, Z, parameters);   
    MU               = Viscosity(T,RHO);
    VELOCITY         = Velocity(F, RHO, parameters);
    dp               = parameters{5};                                    % Paritcle diameter
    Re               = dp .* RHO .* VELOCITY ./ MU .* 1.3;

    a                = parameters{44};
    b                = parameters{45};
    c                = parameters{46};

    D                =  a -  b * Re + c  * F * 10^5;
    %D                = max(D,0) ;
end