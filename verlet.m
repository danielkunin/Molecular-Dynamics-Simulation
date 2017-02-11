function [ pos, vel ] = verlet( pos, vel, mass, charge, connect, step, k0 )
%   pos = 3 x n matrix of current positions (m)
%   vel =  3 x n matrix of current velocity of atoms (m/s)
%   charge = n length vector of atom charge (elementary charge)
%   mass = n length vector of atom mass (kg)
%   connect = n x n array representing initial bond distance
%   step = time step of simulation (s)
%   k0 = bond strength constant (N/M)
%   ----------------------------------------------
%   pos = 3 x n matrix of updated positions (m)
%   vel =  3 x n matrix of updated velocity of atoms (m/s)

% calculate current force vector
f1 = force(pos, charge, connect, k0);

% update position
pos = pos + (step * vel) + 0.5 * step ^ 2 * f1 ./ repmat(mass, 3, 1);

% calculate next force vector
f2 = force(pos, charge, connect, k0);

% update velocity
vel = vel + 0.5 * step * (f1 + f2) ./ repmat(mass, 3, 1);

end

