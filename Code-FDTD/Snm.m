function [Q, Snm_film_excited_silicon] = Snm(m, f_n, gama_n, omega_n, wp, delta_t, N_points_x_all, N_points_y_all, N_points_z_all, N_points_x_film, N_points_y_main, N_points_z_main, x_start_PML_low, x_end_PML_low, N_points_PML, ...
    x_start_Substract, x_end_Substract, x_start_SiO2, x_end_SiO2, x_start_film, x_end_film, x_start_air, x_end_air, x_start_PML_high, x_end_PML_high, order_Snm)
global exp_neg_alfa_delta_t exp_neg_alfa_delta_t_ cos_beta_delta_t sin_beta_delta_t_gama_delta_t
Snm_air = 0;

alfa = gama_n / 2;
beta_squared = omega_n^2 - gama_n^2 / 4;
gama_beta = f_n * wp.^2;
if order_Snm == 1 %%% first order talor expansion
    exp_neg_alfa_delta_t = 1 / (1 + alfa * delta_t);
    exp_neg_alfa_delta_t_ = 1 - alfa * delta_t;
    cos_beta_delta_t = 1 - (beta_squared * delta_t^2) / 2;
    sin_beta_delta_t_gama_delta_t = gama_beta * delta_t^2;
elseif order_Snm == 2 %%% second order talor expansion
    exp_neg_alfa_delta_t = 1 / (1 + alfa * delta_t + (alfa^2 * delta_t^2) / 2);
    exp_neg_alfa_delta_t_ = 1 - alfa * delta_t + (alfa^2 * delta_t^2) / 2;
    cos_beta_delta_t = 1 - (beta_squared * delta_t^2) / 2 + (beta_squared * delta_t^2)^2 / 24;
    sin_beta_delta_t_gama_delta_t = gama_beta * delta_t^2 - (beta_squared * gama_beta * delta_t^4) / 6;
elseif order_Snm == 3 %%% third order talor expansion
    exp_neg_alfa_delta_t = 1 / (1 + alfa * delta_t + (alfa^2 * delta_t^2) / 2 + (alfa^3 * delta_t^3) / 6);
    exp_neg_alfa_delta_t_ = 1 - alfa * delta_t + (alfa^2 * delta_t^2) / 2 - (alfa^3 * delta_t^3) / 6;
    cos_beta_delta_t = 1 - (beta_squared * delta_t^2) / 2 + (beta_squared * delta_t^2)^2 / 24 - (beta_squared * delta_t^2)^3 / 720;
    sin_beta_delta_t_gama_delta_t = gama_beta * delta_t^2 - (beta_squared * gama_beta * delta_t^4) / 6 + (beta_squared^2 * gama_beta * delta_t^6) / 120;
end
if m == 1
    Snm_film_excited_silicon = (2 * exp_neg_alfa_delta_t * cos_beta_delta_t) * ones(N_points_x_film+1, N_points_y_main+1, N_points_z_main+1);
elseif m == 2
    Snm_film_excited_silicon = (exp_neg_alfa_delta_t * exp_neg_alfa_delta_t_) * ones(N_points_x_film+1, N_points_y_main+1, N_points_z_main+1);
elseif m == 3
    Snm_film_excited_silicon = exp_neg_alfa_delta_t * sin_beta_delta_t_gama_delta_t;
end
Snm_glass = 0;
Snm_substract = 0;

Q = zeros(N_points_x_all, N_points_y_all, N_points_z_all);
%%% set up the lorentz-Drude term S11 for the PML layers
Q(x_start_PML_low:x_end_PML_low, 1:N_points_PML+1, :) = Snm_substract; %%% the left part of PML layer
Q(x_start_Substract:x_end_Substract, 1:N_points_PML+1, :) = Snm_substract; %%% the left part of PML layer
Q(x_start_SiO2:x_end_SiO2, 1:N_points_PML+1, :) = Snm_glass; %%% the left part of PML layer
Q(x_start_film:x_end_film, 1:N_points_PML+1, N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = repmat_Y ...
    (reshape(Snm_film_excited_silicon(:, 1, :), N_points_x_film+1, N_points_z_main+1), N_points_PML + 1);%%% the left part of PML layer
Q(x_start_air:x_end_air, 1:N_points_PML+1, :) = Snm_air; %%% the left part of PML layer
Q(x_start_PML_high:x_end_PML_high, 1:N_points_PML+1, :) = Snm_air; %%% the left part of PML layer

Q(x_start_PML_low:x_end_PML_low, N_points_y_all-(N_points_PML):N_points_y_all, :) = Snm_substract; %%% the right part of PML layer
Q(x_start_Substract:x_end_Substract, N_points_y_all-(N_points_PML):N_points_y_all, :) = Snm_substract; %%% the right part of PML layer
Q(x_start_SiO2:x_end_SiO2, N_points_y_all-(N_points_PML):N_points_y_all, :) = Snm_glass; %%% the right part of PML layer
Q(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = ...
    repmat_Y(reshape(Snm_film_excited_silicon(:, N_points_y_main+1, :), N_points_x_film+1, N_points_z_main+1), N_points_PML+1);%%% the right part of PML layer
Q(x_start_air:x_end_air, N_points_y_all-(N_points_PML):N_points_y_all, :) = Snm_air; %%% the right part of PML layer
Q(x_start_PML_high:x_end_PML_high, N_points_y_all-(N_points_PML):N_points_y_all, :) = Snm_air; %%% the right part of PML layer

Q(N_points_x_all-(N_points_PML):N_points_x_all, :, :) = Snm_air; %%% the upper part of PML layer
Q(1:N_points_PML+1, :, :) = Snm_substract; %%% the lower part of PML layer

%%% set up S11 for the PML layers, back
Q(x_start_film:x_end_film, N_points_PML+1+1:N_points_y_all-(N_points_PML + 1), 1:N_points_PML+1) = repmat_Z ...
    (Snm_film_excited_silicon(:, :, 1), N_points_PML + 1);%%% the back part of PML layer
Q(x_start_film:x_end_film, 1:N_points_PML+1, 1:N_points_PML+1) = repmat_YZ ...
    (Snm_film_excited_silicon(:, 1, 1), N_points_PML + 1, N_points_PML + 1);
Q(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, 1:N_points_PML+1) = repmat_YZ ...
    (Snm_film_excited_silicon(:, N_points_y_main+1, 1), N_points_PML + 1, N_points_PML + 1);
%%% set up S11 for the PML layers, front
Q(x_start_film:x_end_film, N_points_PML+1+1:N_points_y_all-(N_points_PML + 1), N_points_z_all-(N_points_PML):N_points_z_all) = repmat_Z ...
    (Snm_film_excited_silicon(:, :, N_points_z_main+1), N_points_PML + 1);%%% the front part of PML layer
Q(x_start_film:x_end_film, 1:N_points_PML+1, N_points_z_all-(N_points_PML):N_points_z_all) = repmat_YZ ...
    (Snm_film_excited_silicon(:, 1, N_points_z_main+1), N_points_PML + 1, N_points_PML + 1);
Q(x_start_film:x_end_film, N_points_y_all-(N_points_PML):N_points_y_all, N_points_z_all-(N_points_PML):N_points_z_all) = repmat_YZ ...
    (Snm_film_excited_silicon(:, N_points_y_main+1, N_points_z_main+1), N_points_PML + 1, N_points_PML + 1);
%%%

%%% set up S11 for other parts of the structure
Q(x_start_Substract:x_end_Substract, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), :) = Snm_substract; %%% Substract
Q(x_start_SiO2:x_end_SiO2, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), :) = Snm_glass; %%% SiO2
Q(x_start_film:x_end_film, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), N_points_PML+1+1:N_points_z_all-(N_points_PML + 1)) = ...
    Snm_film_excited_silicon;%%% film
Q(x_start_air:x_end_air, (N_points_PML + 1 + 1):N_points_y_all-(N_points_PML + 1), :) = Snm_air; %%% air