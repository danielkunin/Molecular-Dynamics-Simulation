function [ke, pe, te, davg] = simulation( molecule, temp, time)
%   molecule = name of input molecule
%   temp = temperature of simulation (kelvin)
%   time = length of time to run simulation (s)
%   ----------------------------------------------
%   ke = kinetic energy vector for each time step (J)
%   pe = potential energy vector for each time step (J)
%   te = total energy vector for each time step (J)


% initial parameters
[pos, mass, charge, connect, k0] = initial_position(molecule);
vel = initial_velocity(mass, temp);

% set step size and energy vectors
step = 1e-16;%-15
num = round(time / step);
ke = zeros(1,num);
pe = zeros(1,num); 
te = zeros(1,num);
davg = zeros(1,num);

% iterative simulation
for i = 1:num
    ke(i) = kinetic_energy(vel, mass);
    pe(i) = potential_energy(pos, charge, connect, k0);
    te(i) = ke(i) + pe(i);
    davg(i) = avg_dist(pos, connect);
    [pos, vel] = verlet( pos, vel, mass, charge, connect, step, k0 );
end

% plot figures of potential & kinetic energy
figure
plot((0:step:(num-1)*step),ke,'LineWidth',3)
hold on
plot((0:step:(num-1)*step),pe,'LineWidth',3)
xlabel('Time (sec)')
ylabel('Energy (Joules)')
legend('Kinetic Energy', 'Potential Energy')
set(gca,'fontsize',14)

% plot figures of total energy & bond length
figure
yyaxis left
plot((0:step:(num-1)*step),te,'LineWidth',3)
ylabel('Energy (Joules)')
hold on
yyaxis right
plot((0:step:(num-1)*step),davg,'LineWidth',3)
ylabel('Distance (m)')
hold off
xlabel('Time (sec)')
legend('Total Energy', 'Average Bond Length')
set(gca,'fontsize',14)


end

