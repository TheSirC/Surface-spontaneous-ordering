%% 1D Gaussian random rough surface with Gaussian autocovariance function
function [f, x] = Roughness_Defined_1D(N, rL, h, l)
x = linspace(-rL/2, rL/2, N);

% uncorrelated Gaussian random rough surface distribution with mean = 0 and standard deviation = 1
% Z = h.*randn(1,N);
Z = randn(1, N);

% Gaussian filter
% F = exp(-x.^2/(l^2*2)); % wikipedia
F = exp(-(x.^2 * l^2)/4); %publication Kuga & Phu

% correlation of surface using convolution, inverse Fourier transform and normalizing prefactors
% f = sqrt(1/(2*pi))*(1/l)*ifft(fft(Z).*fft(F));  % wikipedia
f = sqrt(1/(4 * pi)) * (h^2 * l) * ifft(fft(Z).*fft(F)); %publication Kuga & Phu
end