function [ f ] = force( pos, charge, connect, k0 )
%   pos = 3 x n matrix of current positions (m)
%   charge = n length vector of atom charge (elementary charge)
%   connect = n x n array representing initial bond distance
%   k0 = bond strength constant (N/M)
%   ----------------------------------------------
%   f = 3 x n matrix of current force vectors (N)

f_bond = force_bond(pos, connect, k0);
f_electro = force_electro(pos, charge);
% f_vdw = force_vdw(pos);
f = f_bond + f_electro;% + f_vdw;

end

function [ f_bond ] = force_bond( pos, connect, k0 )
%   pos = 3 x n matrix of current positions (m)
%   connect = n x n array representing initial bond distance
%   k0 = bond strength constant (N/M)
%   ----------------------------------------------
%   f_bond = current force vector on molecule by bond (N)

n = size(connect,2);
f_bond = zeros(3,n);
for k = 1:n
    e_bond = 0;
    for i = 1:k-1
        if connect(i,k)
            rik = norm(pos(:,i) - pos(:,k));
            e_prime = 2 * k0 * (rik - connect(i,k));
            e_bond = e_bond - e_prime / rik;
        end
    end
    for i = k+1:n
        if connect(k,i)
            rki = norm(pos(:,k) - pos(:,i));
            e_prime = 2 * k0 * (rki - connect(k,i));
            e_bond = e_bond + e_prime / rki;
        end
    end
    f_bond(:,k)  = -1 * e_bond * pos(:,k);
end

end


function [ f_electro ] = force_electro( pos, charge )
%   pos = 3 x n matrix of current positions (m)
%   charge = n length vector of atom charge (elementary charge)
%   ----------------------------------------------
%   f_bond = current force vector on molecule by electro (N)

n = size(pos,2);
f_electro = zeros(3,n);
coulomb = 1 / (4 * pi * 8.85e-12);
for k = 1:n
    e_electro = 0;
    for i = 1:k-1
        rik = norm(pos(:,i) - pos(:,k));
        e_prime = -1 * coulomb * charge(i) * charge(k) / rik^2;
        e_electro = e_electro - e_prime / rik;
    end
    for i = k+1:n
        rki = norm(pos(:,k) - pos(:,i));
        e_prime = -1 * coulomb * charge(i) * charge(k) / rki^2;
        e_electro = e_electro + e_prime / rki;
    end
    f_electro(:,k)  = -1 * e_electro * pos(:,k);
end

end



function [ f_vdw ] = force_vdw( pos )
%   pos = 3 x n matrix of current positions (m)
%   mass = n length vector of atom mass (kg)
%   charge = n length vector of atom charge (elementary charge)
%   connect = n x n array representing initial bond distance
%   ----------------------------------------------
%   f_bond = current force vector on molecule by vdw (N)

n = size(pos,2);
eij = .0658; % for H2 only (not sure if this is correct)
sigij = 3.01e-9; % for H2 only
e_vdw = 0;
for i = 1:n
    for j = i+1:n
        rij = norm(pos(:,i) - pos(:,j));
        e_vdw = e_vdw + eij * ((sigij / rij)^12 - 2 * (sigij / rij)^6);
    end
end

end