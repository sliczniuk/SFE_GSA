function Sat_coe = Compute_kinetic_for_GSA(X, Csolid_percentage_left)

    Di    = Compute_Diffusion_for_GSA(X);
    gamma = Compute_Gamma_for_GSA(X);

    Sat_coe                = Saturation_Concentration(Csolid_percentage_left, gamma, Di);        % Inverse logistic is used to control saturation. Close to saturation point, the Sat_coe goes to zero.
end