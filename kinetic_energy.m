function [ ke ] = kinetic_energy( vel, mass )
%   vel = 3 x n matrix of initial velocity of atoms (m/s)
%   mass = n length vector of mass of atoms in molecule (kg)
%   ----------------------------------------------
%   ke = current kinetic energy of molecule (J)

% calculate square velocities
vel_square = sum(vel .^ 2, 1);

% calculate ke using formula: ke = 1/2mv^2
ke = 0.5 * sum(mass .* vel_square);

end

