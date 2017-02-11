function [ freq, amp ] = find_frequency( data, step, plot )
%   data = sequence data to find frequency and amplitude (m or J)
%   step = time step of sequence data (s)
%   plot = boolean to plot or not plot data and peaks
%   ----------------------------------------------
%   freq = average frequency of data (1/s)
%   amp = average amplitude of data (m or J)

    
% find peaks and locations
[PKS,LOCS] = findpeaks(data);

% get average period in second
period = mean(LOCS(2:end) - LOCS(1:end - 1)) * step;

% find average frequency in 1/second
freq = 1/period;

% find average amplitude
amp = mean(PKS) - mean(data);

% plot data and peaks
if plot
    n = length(data);
    figure
    plot((0:step:(n-1)*step),data,'LineWidth',3)
    hold on
    plot(LOCS*step,PKS,'o')
    ylabel('Energy (Joules)')
    xlabel('Time (sec)')
    set(gca,'fontsize',14)
end

end

