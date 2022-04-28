function H = H_sum_of_sines(xx, yy, k, D, A, omega, phi, t)
%   xx: M x 1 vector, x-coordinates of the positions
%   yy: M x 1 vector, y-coordinates of the positions
%   k: N x 1, distortion constant
%   D: N x 1 vector, directions of the sine waves in radians
%   A: N x 1 vector, amplitudes of the sine waves
%   omega: N x 1 vector,wavenumber, determine wavelength the of the sine waves
%   2Pi/Labmda
%   phi: N x 1 vector, determine the moving speeds of the sine waves
%   t: T x 1 vector, sampling timings
%
%   H: M x T matrix, height of the wave at each position and time
    T = size(t, 1);
    M = size(xx, 1);
    N = size(D, 1);
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
    ORDER = zeros(1, 1, N);
    ORDER(1, 1, :) = k;
    ORDER = repmat(ORDER, M, T, 1);
%     angle = [xx yy] * [cos(D) sin(D)]' .* repmat(omega', M, 1)  + ...
%         repmat(phase', M, 1);
    angle = SPATIAL_PROJECTION .* OMEGA + PHASE;
    H = 2 * AMPLITUDE .* ((((sin(angle) + 1) / 2) .^ ORDER) - 1/2);
    H = sum(H, 3);
end