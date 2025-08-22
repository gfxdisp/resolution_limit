%% Required pixel resolution for picture heights (H)

clear all; close all;

ds = readtable("../data/resolution_limit_median_data.csv", 'Delimiter', ',');

show_exponential_fit = true;
if show_exponential_fit
    load('Watson2018_params.mat');
end         

plot_population_var = true;

legendHandle = [];

% colorbar
cmin = 0;
cmax = 1;
colormap('parula');  % Choose a colormap that suits your preference

clf,
gap = [0.05 ,0.05];
marg_h = [0.2, 0.01]; % [lower upper]
marg_w = [0.18, 0.01]; % [left right]
ha = tight_subplot(1,1, gap, marg_h, marg_w);

set(gcf, 'defaultAxesFontName', 'Arial');
set(gcf, 'defaultTextFontName', 'Arial');
FontSize = 8;

fig = gcf;
fig_size = [7.5 4.5]*0.49*0.93;
fig.Units = 'inches';
fig.Position = [2, 2, fig_size];
fig.PaperUnits = 'inches';
fig.PaperSize = fig_size;

hold on,

axes(ha(1));
hh = [];

hold on,

ppd_low = 0;
ppd_high = 120;

px_low = 500;
px_high = 18000;

pic_heights = linspace(0.1, 6.2, 500);

cc = 1; % achormatic only
foveal_ppd_pred = Secc_watson2018(pars_all(cc, 1), 0, pars_all(cc, 2:4));

if show_exponential_fit

    if plot_population_var
        px_bins = linspace(px_low, px_high, 1000);
    end

    % Generate pdf
    P_all = [];

    ds_cc = ds((ds.color_direction == cc) & (ds.eccentricity == 0), :);
    ds_sigma = (ds_cc.pi_high_fppd-ds_cc.pi_low_fppd)./(2*1.96);

    for i = 1:length(pic_heights)
        ppd_bins = px_bins./(2*rad2deg(atan(16./(2.*9.*pic_heights(i)))));

        x_vals = -ppd_to_fppd(ppd_bins);
        central_value = -ppd_to_fppd(foveal_ppd_pred);
        gaussian = @(standard_deviation) (1 / (sqrt(2 * pi) * standard_deviation))... 
            * exp(-((x_vals - central_value).^2 / (2 * standard_deviation^2)));
    
        pdf_gauss = gaussian(ds_sigma);

        cdf_gauss = cumsum(pdf_gauss);
        cdf_gauss = 1-(cdf_gauss./max(cdf_gauss));
        P_all = [P_all, (cdf_gauss)'];

    end

    imagesc(pic_heights, (px_bins), P_all);
    clim([cmin, cmax]);  % Apply the global color scale
    
    prctile_vals = [0.95, 0.5];
    ls = {':', '--', '--'};
    
    hold on,
    if 1 % plot lines of intervals
        percentile_all = 1-P_all;
            
        [~, index_zero] = min(percentile_all);

        % Create a logical index for each column
        [col_indices, row_indices] = meshgrid(1:size(percentile_all, 2), 1:size(percentile_all, 1));
        
        lw = 1;
        for pp = 1:length(prctile_vals)

            logical_index = row_indices >= index_zero;
            percentile_all_temp = percentile_all;
            percentile_all_temp(~logical_index) = NaN;
            [~, index_high] = min(abs(percentile_all_temp-prctile_vals(pp)));
            index_high(index_high == 1) = NaN;
            
            plot(pic_heights(~isnan(index_high)),... 
                px_bins(index_high(~isnan(index_high))), ls{pp},... 
                'Color', 'white', 'LineWidth', lw, 'HandleVisibility','off');

            if prctile_vals(pp) == 0.5
                plt_name = sprintf('%dth percentile (median)', round(prctile_vals(pp)*100));
            else
                plt_name = sprintf('%dth percentile', round(prctile_vals(pp)*100));
            end
            
            hh(pp) = plot(NaN, NaN, ls{pp}, 'color', 'k', 'LineWidth',lw,...
                'DisplayName', plt_name);
            
        end

    end


end

if 1 % plot standard res lines
    standard_px = [1920, 3840, 7680, 15360];
    px_text = {'FHD', '4k', '8k', '16k'};
    itu_ranges = [3.2, 3.2; 1.6, 3.2; 0.8, 3.2; 0.8, 3.2];
    for i = 1:length(standard_px)
        plot([0, 20], [standard_px(i), standard_px(i)], ':r', 'HandleVisibility', 'off');
        text(max(pic_heights)-0.8, standard_px(i)+1000, px_text{i}, 'Color','red','FontSize',FontSize);
        if i == 1
            scatter(itu_ranges(i, 1), standard_px(i), 30, 'r', 'filled', 's')
        elseif i < 4
            plot([itu_ranges(i, 1), itu_ranges(i, 2)], [standard_px(i), standard_px(i)], '-r',...
                'LineWidth', 2, 'HandleVisibility', 'off');
        end
    end
    
    
end
    
xlabel({'Viewing distance in display heights (H)'});
ylabel({'Horizontal display resolution'});
ax = gca;
yLabel = ax.YLabel;

grid on
ax = gca; % Get the current axes to modify properties
ax.GridColor = [1 1 1].*0.5;  % Choose grid lines color
ax.GridLineStyle = '-'; % Choose grid lines style
ax.GridAlpha = 0.4; % Choose grid lines transparency
ax.Layer = 'top'; % This places the grid on top of the image
ax.MinorGridColor = [1 1 1].*0.5;  % Set minor grid lines color
ax.MinorGridLineStyle = '-'; % Set minor grid lines style
ax.MinorGridAlpha = 0.2; % Set minor grid lines transparency
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';

ylim([px_low, px_high]);
xlim([min(pic_heights), max(pic_heights)]);

ytick = round(linspace(px_low*1.1, px_high*0.9, 5), 2, 'significant');
set(gca, 'YTick', ytick);
set(gca, 'YTickLabel', ytick);

xtick = 1:1:10;
set(gca, 'XTick', xtick);
set(gca, 'XTickLabel', xtick);

% Find all text objects in the figure
textObjects = findall(gcf, 'Type', 'text');
set(textObjects, 'FontSize', FontSize);

axesObjects = findall(gcf, 'Type', 'axes');
set(axesObjects, 'FontSize', FontSize);
