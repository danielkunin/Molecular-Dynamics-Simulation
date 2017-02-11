function [ vel ] = initial_velocity( mass, temp )
%   mass = n length vector of mass of atoms in molecule (kg)
%   temp = the current temperature (kelvin)
%   ----------------------------------------------
%   vel = 3 x n matrix of initial velocity of atoms (m/s)

% number of atoms in molecule
n = size(mass,2);

% Boltzmann's constant
kb = physconst('Boltzmann');

% vector of variances
var = 3 * kb * temp ./ mass;

% vector of speeds
speed = abs(normrnd(0,sqrt(var)));

% matrix of 3D vectors
vect = randn(3,n);

% matrix of 3D unit vectors
unit = vect ./ repmat(sqrt(sum(abs(vect).^2,1)), 3, 1);

% matrix of 3D unit vectors scaled by speeds
vel =  unit .* repmat(speed, 3, 1);
end

