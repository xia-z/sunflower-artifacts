function [Hx, Hy] = H_partial_sum_of_sines(xx, yy, k, D, A, omega, phi, t)
%   xx: M x 1 vector, x-coordinates of the positions
%   yy: M x 1 vector, y-coordinates of the positions
%   k: scalar, distortion constant
%   D: N x 1 vector, directions of the sine waves in radians
%   A: N x 1 vector, amplitudes of the sine waves
%   omega: N x 1 vector, determine the wavelengths of the sine waves
%   phi: N x 1 vector, determine the moving speeds of the sine waves
%   t: T x 1 vector, sampling timings
%
%   Hx: M x T matrix, x-partial derivative of the height of the wave at
%       each position and time
%   Hy: M x T matrix, y-partial derivative of the height of the wave at
%       each position and time
    M = size(xx, 1);
    N = size(D, 1);
    T = size(t, 1);
    SPATIAL_PROJECTION = zeros(M, 1, N);
    SPATIAL_PROJECTION(:, 1, :) = [xx yy] * [cos(D) sin(D)]';
    SPATIAL_PROJECTION = repmat(SPATIAL_PROJECTION, 1, T, 1);
    OMEGA = zeros(1, 1, N);
    OMEGA(1, 1, :) = omega;
    OMEGA = repmat(OMEGA, M, T, 1);
    PHASE = zeros(1, T, N);
    PHASE(1, :, :) = t * phi';
    PHASE = repmat(PHASE, M, 1, 1);
    AMPLITUDE = zeros(1, 1, N);
    AMPLITUDE(1, 1, :) = A;
    AMPLITUDE = repmat(AMPLITUDE, M, T, 1);
%     angle = [xx yy] * [cos(D) sin(D)]' .* repmat(omega', M, 1)  + ...
%         repmat(phase', M, 1);
    angle = SPATIAL_PROJECTION .* OMEGA + PHASE;
    tmp = k * OMEGA .* AMPLITUDE .* ...
        (((sin(angle) + 1) / 2) .^ (k - 1)) .* cos(angle);
    Dx = zeros(1, 1, N);
    Dy = zeros(1, 1, N);
    Dx(1, 1, :) = cos(D);
    Dy(1, 1, :) = sin(D);
    Dx = repmat(Dx, M, T, 1);
    Dy = repmat(Dy, M, T, 1);
    Hx = tmp .* Dx;
    Hy = tmp .* Dy;
    Hx = sum(Hx, 3);
    Hy = sum(Hy, 3);
end