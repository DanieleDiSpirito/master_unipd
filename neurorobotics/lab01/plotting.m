function plotting(filename)
    [data, h] = sload(filename);
    means = mean(data(:, :));
    std = var(data(:, :)).^0.5;
    subplot(2,1,1)
    hold on;
    for i = 1:size(data, 2)-1
        ksdensity(data(:, i))
    end
    subplot(2,1,2)
    hold on;
    xlabel('Mean')
    ylabel('Std Dev')
    plot(means, std, 'X')
end