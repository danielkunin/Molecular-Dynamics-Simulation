%% Make Stem Plot of Frequency and Amplitude

molecules = {'hydrogen','nitrogen','oxygen','flourine','chlorine',...
             'bromine', 'iodine'};
         
n = length(molecules);
frequency = zeros(1,n);
amplitude = zeros(1,n);
step = 1e-16;

for i=1:n
    [ke,pe,te,davg] = simulation_threshold(char(molecules(i)),100);
    [frequency(i),amplitude(i)] = find_frequency(davg,step,0);
end

figure;
hold on
for i=1:n
    stem(frequency(i),amplitude(i),'LineWidth',1.5);
end
hold off
legend('hydrogen','nitrogen','oxygen','flourine','chlorine','bromine', 'iodine')
xlabel('Frequency (1/sec)')
ylabel('Amplitude (m)')
set(gca,'fontsize',14)

%% Make Stem Plot of Frequency and Amplitude as function of Temp

n = 5000;         
frequency = zeros(1,n/10 + 1);
amplitude = zeros(1,n/10 + 1);
step = 1e-16;

ind = 0;
for i=0:10:5000
    ind = ind + 1
    [~,~,~,davg1] = simulation_threshold('hydrogen',i);
    [freq1,ampl1] = find_frequency(davg1,step,0);
    close ALL
    
    [~,~,~,davg2] = simulation_threshold('hydrogen',i);
    [freq2,ampl2] = find_frequency(davg2,step,0);
    close ALL
    
    [~,~,~,davg3] = simulation_threshold('hydrogen',i);
    [freq3,ampl3] = find_frequency(davg3,step,0);
    close ALL
    
    frequency(ind) = (freq1 + freq2 + freq3) / 3;
    amplitude(ind) = (ampl1 + ampl2 + ampl3) / 3;
end

figure;
plot((0:10:5000),frequency,'LineWidth',3,'Color',[ 0.8500    0.3250    0.0980 ])
ylabel('Frequency (1/sec)')
xlabel('Temperature (K)')
set(gca,'fontsize',14)


figure;
plot((0:10:5000),amplitude,'LineWidth',3,'Color',[ 0.4660    0.6740    0.1880 ])
ylabel('Amplitude (m)')
xlabel('Temperature (K)')
set(gca,'fontsize',14)


%% Make Stem Plot of Frequency and Amplitude as function of Temp

n = 5000;         
iterations = zeros(1,n/10 + 1);
step = 1e-16;

ind = 0;
for i=0:10:5000
    ind = ind + 1
    [~,~,~,davg1] = simulation_threshold('hydrogen',i);
    close ALL
    
    [~,~,~,davg2] = simulation_threshold('hydrogen',i);
    close ALL
    
    [~,~,~,davg3] = simulation_threshold('hydrogen',i);
    close ALL
    
    iterations(ind) = (length(davg1) + length(davg2) + length(davg3)) / 3;
end

figure;
plot((0:10:5000),iterations * step,'LineWidth',3)
ylabel('Time (sec)')
xlabel('Temperature (K)')
set(gca,'fontsize',14)

%% Make Stem Plot of Frequency and Amplitude

% molecules = {'hydrogen','nitrogen','oxygen','flourine','chlorine',...
%              'bromine', 'iodine'};
         
molecules = {'nitride', 'monoxide', 'nitric','flouride', 'chloride'};
         
n = length(molecules);
m = 10;
frequency = zeros(1,m*n);
amplitude = zeros(1,m*n);
step = 1e-16;

for i=0:n-1
    i
    for j=1:m
        [ke,pe,te,davg] = simulation(char(molecules(i+1)),100,1e-13);
        [frequency(i * m + j),amplitude(i * m + j)] = find_frequency(davg,step,0);
    end
    close ALL
end

figure;
hold on
for i=0:n-1
    plot(frequency(i*m + 1:i*m+m),amplitude(i*m+1:i*m+m),'.','markersize',30)
end
hold off
legend('nitride', 'monoxide', 'nitric','flouride', 'chloride')
xlabel('Frequency (1/sec)')
ylabel('Amplitude (m)')
grid
set(gca,'fontsize',14)