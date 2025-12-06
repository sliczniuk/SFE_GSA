function Upsilon = Compute_Gamma_for_GSA(X)

    T = X(:,1); P = X(:,2); F = X(:,3);

    Parameters_table = readtable('Parameters.csv') ;                     % Table with prameters
    parameters       = num2cell(Parameters_table{:,3});                  % Parameters within the model + (m_max), m_ratio, sigma

    Z                = Compressibility(T, P, parameters);
    RHO              = rhoPB_Comp(T, P, Z, parameters);   
    MU               = Viscosity(T,RHO);
    VELOCITY         = Velocity(F, RHO, parameters);
    dp               = parameters{5};                                    % Paritcle diameter
    Re               = dp .* RHO .* VELOCITY ./ MU .* 1.3;

    a = parameters{47};%3.158;
    b = parameters{48};%11.922;
    c = parameters{49};%0.6868;

    Upsilon  = a + b * Re - c * F * 10^5;
end