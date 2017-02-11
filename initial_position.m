function [ pos, mass, charge, connect, k0 ] = initial_position( molecule )
%   molecule =  name of molecule
%   ----------------------------------------------
%   pos = 3 x n matrix of intial positions (m)
%   mass = n length vector of atom mass (kg)
%   charge = n length vector of atom charge (elementary charge)
%   connect = n x n array representing initial bond distance
%   k0 = bond strength constant (N/M)

% path to PDB file
path = strcat('molecules/', molecule,'.pdb.txt');

% read molecule PDB file
pdbstruct = pdbread(path);

% Get array of Atoms 
atoms = pdbstruct.Model.Atom;
n = length(atoms);
pos = zeros(3,n);
mass = zeros(1,n);
charge = zeros(1,n);

% loop through atoms and get pos, mass, charge 
angstrom_to_meter = 1e-9;
elem_charge = 1.60217656535e-19;
for i = 1:n
    pos(1,i) = atoms(i).X * angstrom_to_meter;
    pos(2,i) = atoms(i).Y * angstrom_to_meter;
    pos(3,i) = atoms(i).Z * angstrom_to_meter;
    mass(i) = get_mass(atoms(i).AtomName);
    charge(i) = str2double(atoms(i).charge) * elem_charge;
end

% create connectivity matrix
structure = pdbstruct.Connectivity;
num = length(structure);
connect = zeros(n,n);

% loop through structure adding bonds
for i = 1:num
    atom_1 = structure(i).AtomSerNo;
    bonds = structure(i).BondAtomList;
    for j = 1:length(bonds)
        atom_2 = bonds(j);
        if (atom_1 > 0) && (atom_2 > 0)
            r0 = norm(pos(:,atom_1) - pos(:,atom_2));
            connect(atom_1,atom_2) = r0;
            connect(atom_2,atom_1) = r0;
        end
    end
end

k0 = get_bond_strength(molecule);
    
end


function [ mass ] = get_mass( atom )
%   atom = char of atomic symbol
%   ----------------------------------------------
%   mass = mass of atom (kg)

% Object of molecular symbols
symbol = {'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V',...
    'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Yt','Zr','Nb','Mo','Tc','Ru','Rh',...
    'Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho'...
    'Er','Tm','Yt','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',...
    'Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','#'};

% Object of molecular mass
mole = [1.0079,4.0026,6.941,9.01218,10.81,12.011,14.0067,15.9994,18.9984,20.179,22.98977,24.305,...
    26.98154,28.086,30.97376,32.06,35.453,39.948,39.098,40.08,44.9559,47.9,50.9414,51.996,...
    54.938,55.847,58.9332,58.7,63.546,65.38,69.72,72.59,74.9216,78.96,79.904,83.8,85.4678,...
    87.62,88.9059,91.22,92.9064,95.94,98.9062,101.07,102.9055,106.4,107.868,112.4,114.82,...
    118.69,121.75,127.6,126.9045,131.3,132.9054,137.34,138.9055,140.12,140.9077,144.24,147,...
    150.4,151.96,157.25,158.9254,162.5,164.9304,167.26,168.9342,173.04,174.97,178.49,180.9479,...
    183.85,186.207,190.2,192.22,195.09,196.9665,200.59,204.37,207.2,208.9804,210,210,222,...
    223,226.0254,227,232.0381,231.0359,238.029,237.0482,244,243,247,247,251,254,257,0,...
    0,0,0,0,0,0,0,0,0];

% Create struct of atom information
molar = containers.Map(symbol,mole);

% Obtain mass and convert moles to kg
avos_num = 6.0221e23;
mass = molar(atom) / 1000 / avos_num;

end

function [ k0 ] = get_bond_strength( molecule )
%   molecule = name of molecule
%   ----------------------------------------------
%   k0 = bond strength constant (N/M)

% Object of molecular symbols
name = {'hydrogen','nitrogen','oxygen','flourine','chlorine',...
          'bromine', 'iodine', 'nitride', 'monoxide', 'nitric',...
          'flouride', 'chloride'};

% Object of constants
value = [488.09, 2510.53, 1322.84, 498.33, 337.34, 231.24, 167.76,...
         861.11, 2041.15, 1791.51, 929.06, 578.68];


% Create struct of molecule information
constants = containers.Map(name,value);

% Obtain bond strength constant
k0 = constants(molecule);

end