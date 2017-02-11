function [ pe ] = potential_energy( pos, charge, connect, k0 )
%   pos = 3 x n matrix of current positions (m)
%   charge = n length vector of atom charge (elementary charge)
%   connect = n x n array representing initial bond distance
%   k0 = bond strength constant (N/M)
%   ----------------------------------------------
%   pe = current potential energy of molecule (J)

e_bond = energy_bond(pos, connect, k0);
e_electro = energy_electro(pos, charge);
%e_vdw = energy_vdw(pos);
pe = e_bond + e_electro;% + e_vdw;

end



function [ e_bond ] = energy_bond( pos, connect, k0 )
%   pos = 3 x n matrix of current positions (m)
%   connect = n x n array representing initial bond distance
%   ----------------------------------------------
%   e_bond = current bond energy of molecule (J)

n = length(connect);
e_bond = 0;
for i = 1:n
    for j = 1:n
        if connect(i,j)
            rij = norm(pos(:,i) - pos(:,j));
            e_bond = e_bond + k0 * (rij - connect(i,j))^2;
        end
    end
end

end


function [ e_electro ] = energy_electro( pos, charge )
%   pos = 3 x n matrix of current positions (m)
%   charge = n length vector of atom charge (elementary charge)
%   ----------------------------------------------
%   e_electro = current electro potential energy of molecule (J)

n = size(pos,2);
coulomb = 1 / (4 * pi * 8.85e-12);
e_electro = 0;
for i = 1:n
    for j = 1:n
        if i ~= j
            rij = norm(pos(:,i) - pos(:,j));
            e_electro = e_electro + coulomb * charge(i) * charge(j) / rij;
        end
    end
end

end



function [ e_vdw ] = energy_vdw( pos )
%   pos = 3 x n matrix of current positions (m)
%   mass = n length vector of atom mass (kg)
%   charge = n length vector of atom charge (elementary charge)
%   connect = n x n array representing initial bond distance
%   ----------------------------------------------
%   pe = current potential energy of molecule (J)

n = length(pos);
eij = 8.6; % for H2 only (not sure if this is correct)
sigij = 3.01e-10; % for H2 only
e_vdw = 0;
for i = 1:n
    for j = 1:n
        if i ~= j
            rij = norm(pos(:,i) - pos(:,j));
            e_vdw = e_vdw + eij * ((sigij / rij)^12 - 2 * (sigij / rij)^6);
        end
    end
end

end