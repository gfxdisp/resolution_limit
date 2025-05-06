%% Plot prediction interval from resolution limit data

clear all; close all;

eccs = [0, 10, 20];
ds = readtable("../data/resolution_limit_median_data.csv", 'Delimiter', ',');

show_exponential_fit = true;
if show_exponential_fit
    load('Watson2018_params.mat');
end

plot_population_var = true;

log_axis = true;
if log_axis
    axis_string = 'log';
else
    axis_string = 'lin';
end

plot_cdf = true;
if plot_cdf
    prb_string = 'cdf';
else
    prb_string = 'pdf';
end

if log_axis
    ppd_low = 0.3;
    ppd_high = 200;
else
    ppd_low = 0;
    ppd_high = 120;
end

legendHandle = [];

% colorbar
cmin = 0;
cmax = 1;
colormap('viridis'); 

% figure,
gap = [0.02 ,0.04];
marg_h = [0.135, 0.15]; % [lower upper]
marg_w = [0.09, 0.015]; % [left right]
ha = tight_subplot(1,3, gap, marg_h, marg_w);

set(gcf, 'defaultAxesFontName', 'Arial');
set(gcf, 'defaultTextFontName', 'Arial');
FontSize = 8;

fig = gcf;
fig_size = [7.5 2.7]*0.93;
fig.Units = 'inches';
fig.Position = [2, 2, fig_size];
fig.PaperUnits = 'inches';
fig.PaperSize = fig_size;

hold on,

labels = {'Achromatic', 'Red-green', 'Yellow-violet'};

for cc = 1:3
    axes(ha(cc));
    hh = [];
    hold on,

    if show_exponential_fit
    
        eccentricity_plot = linspace(0, 50, 2000);

        fittedCurve = Secc_watson2018(pars_all(cc, 1), eccentricity_plot, pars_all(cc, 2:4));
        
        if plot_population_var
            
            if log_axis
                ppd_bins = logspace(log10(ppd_low), log10(ppd_high), 1000);
            else
                ppd_bins = linspace((ppd_low), (ppd_high), 1000);
            end

            % Generate pdf
            P_all = [];
            % Calculate the prediction intervals for each point in 'eccentricity_plot'
            for i = 1:length(eccentricity_plot)

                ds_cc = ds(ds.color_direction == cc, :);
                ds_eccs = ds_cc.eccentricity;
                ds_sigma = (ds_cc.pi_high_fppd-ds_cc.pi_low_fppd)./(2*1.96);

                x_vals = -ppd_to_fppd(ppd_bins);
                central_value = -ppd_to_fppd(fittedCurve(i));

                gaussian = @(standard_deviation) (1 / (sqrt(2 * pi) * standard_deviation)) * exp(-((x_vals - central_value).^2 / (2 * standard_deviation^2)));

                pdf_gauss = gaussian(interp1(ds_eccs, ds_sigma, eccentricity_plot(i)));

                pdf_gauss = pdf_gauss./max(pdf_gauss);
                if plot_cdf
                    cdf_gauss = cumsum(pdf_gauss);
                    cdf_gauss = 1-(cdf_gauss./max(cdf_gauss));
                    P_all = [P_all, (cdf_gauss)'];
                else
                    P_all = [P_all, (pdf_gauss)'];
                end

            end

        end
        
        if log_axis
            intensity = log10(ppd_bins);
        else
            intensity = ppd_bins;
        end
        
        imagesc(eccentricity_plot, (intensity), P_all);
        clim([cmin, cmax]);  % Apply the global color scale
        
        prctile_vals = [0.95, 0.75, 0.5];
        ls = {':', '-.', '--'};
        
        hold on,
        if 1 % plot lines of intervals

            percentile_all = 1-P_all;
            
            [~, index_zero] = min(percentile_all);

            % Create a logical index for each column
            [col_indices, row_indices] = meshgrid(1:size(percentile_all, 2), 1:size(percentile_all, 1));
            
            lw = 1;
            for pp = 1:length(prctile_vals)

                if ~plot_cdf
                    logical_index = row_indices <= index_zero;
                    percentile_all_temp = percentile_all;
                    percentile_all_temp(~logical_index) = NaN;
                    [~, index_low] = min(abs(percentile_all_temp-prctile_vals(pp)));
                    index_low(index_low == 1) = NaN;

                    plot(eccentricity_plot(~isnan(index_low)),... 
                        intensity(index_low(~isnan(index_low))), ls{pp},... 
                        'Color', 'white', 'LineWidth', lw, 'HandleVisibility', 'off');
                end
    
                logical_index = row_indices >= index_zero;
                percentile_all_temp = percentile_all;
                percentile_all_temp(~logical_index) = NaN;
                [~, index_high] = min(abs(percentile_all_temp-prctile_vals(pp)));
                index_high(index_high == 1) = NaN;
                
                
                plot(eccentricity_plot(~isnan(index_high)),... 
                    intensity(index_high(~isnan(index_high))), ls{pp},... 
                    'Color', 'white', 'LineWidth', lw, 'HandleVisibility','off');
                
                if pp == 3
                    plt_name = sprintf('%dth percentile (median)', round(prctile_vals(pp)*100));
                else
                    plt_name = sprintf('%dth percentile', round(prctile_vals(pp)*100));
                end
                
                hh(pp) = plot(NaN, NaN, ls{pp}, 'color', 'k', 'LineWidth',lw,...
                    'DisplayName', plt_name);
            end

        end

    end

    title(labels{cc}, 'Fontweight', 'Normal');
    
    if cc == 2
        xlabel('Retinal eccentricity (visual degrees)')
    end

    if cc == 1
        ylabel({'Resolution limit', '(pixels per degree)'});
    end

    if cc == 3
        legend boxoff
        legend(hh, 'NumColumns', 3,... 
            'Position', [0.1913    0.8824    0.6279    0.1117]);
    end

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

    ax.YAxis.MinorTickValues = log10([1:1:9, 10:10:90, 100:100:500]); % Set minor y-ticks

    xlim([0, 20]);
    xtick = [0, 5, 10, 15, 20];
    set(gca, 'XTick', xtick);
    set(gca, 'XTickLabel', xtick);
    
    if log_axis
        
        if plot_cdf
            switch cc
                case 1
                    y_lim = (log10([14, 120]));
                case 2
                    y_lim = (log10([5, 113]));
                case 3
                    y_lim = (log10([3, 85]));
            end
        else

            switch cc
                case 1
                    y_lim = (log10([10, 130]));
                case 2
                    y_lim = (log10([3, 130]));
                case 3
                    y_lim = (log10([1, 100]));
            end
        end
        ylim(y_lim);
        ytick = unique(round(logspace(y_lim(1)+0.1, y_lim(2)-0.1, 5), 1, 'significant'));
        set(gca, 'YTick', log10(ytick));

    else
        ytick = [0:20:150];
        set(gca, 'YTick', ytick);
        switch cc
            case 1
                y_lim = (([10, 150]));
            case 2
                y_lim = (([3, 150]));
            case 3
                y_lim = (([1, 100]));
        end
    end
        
    if 1% cc == 1
        set(gca, 'YTickLabel', ytick);
    end

end

% Find all text objects in the figure
textObjects = findall(gcf, 'Type', 'text');
set(textObjects, 'FontSize', FontSize);

axesObjects = findall(gcf, 'Type', 'axes');
set(axesObjects, 'FontSize', FontSize);
