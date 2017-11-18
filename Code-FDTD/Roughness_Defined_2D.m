%% 2D Gaussian random rough surface with Gaussian autocovariance function
function [f, x, y] = Roughness_Defined_2D(N, rL, h, lx, ly)
x = linspace(-rL/2, rL/2, N);
y = linspace(-rL/2, rL/2, N);
[X, Y] = meshgrid(x, y);

% uncorrelated Gaussian random rough surface distribution with mean 0 and standard deviation 1
%Z = h.*randn(N,N);
Z = randn(N, N);

% isotropic surface lx = ly
if lx == ly
    % Gaussian filter
    F = exp(-((X.^2 + Y.^2) * (lx^2 / 4)));
    
    % correlation of surface including convolution, inverse Fourier transform and normalizing prefactors
    f = 1e-5 * ((lx^2 * h^2) / (4 * pi)) * ifft2(fft2(Z).*fft2(F));
    
    % non-isotropic surface lx ~= ly
else
    % Gaussian filter
    F = exp(-((X.^2 * lx^2) / 4)-((Y.^2 * ly^2) / 4));
    
    % correlated surface generation including convolution and inverse Fourier transform
    f = 1e-5 * ((lx * ly * h^2) / (4 * pi)) * ifft2(fft2(Z).*fft2(F));
end
end