%% Required pixel resolution for viewing distance

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
marg_w = [0.15, 0.01]; % [left right]
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

log_axis = true;
log_x_axis = true;

m_to_inch = 39.37;

if log_x_axis
    vd = logspace(log10(0.1), log10(5), 500);
else
    vd = linspace(0.1, 5, 500);
end

d_in = 27;
d_m = d_in./m_to_inch;

if log_axis
    ppi_low = 20;
    ppi_high = 700;
else
    ppi_low = 20;
    ppi_high = 600;
end

cc = 1; % achormatic only
foveal_ppd_pred = Secc_watson2018(pars_all(cc, 1), 0, pars_all(cc, 2:4));

if show_exponential_fit

    pred_ppi = (foveal_ppd_pred*2.*rad2deg(atan(d_m./(2.*vd))))./(m_to_inch*d_m);

    if plot_population_var
        if log_axis
            ppi_bins = logspace(log10(ppi_low), log10(ppi_high), 1000);
        else
            ppi_bins = linspace(ppi_low, ppi_high, 1000);
        end
    end

    % Generate pdf
    P_all = [];
    
    ds_cc = ds((ds.color_direction == cc) & (ds.eccentricity == 0), :);
    ds_sigma = (ds_cc.pi_high_fppd-ds_cc.pi_low_fppd)./(2*1.96);

    for i = 1:length(vd)
        ppd_bins = ppi_bins.*m_to_inch.*d_m./(2*rad2deg(atan(d_m./(2.*vd(i)))));

        x_vals = -ppd_to_fppd(ppd_bins);
        central_value = -ppd_to_fppd(foveal_ppd_pred);
        gaussian = @(standard_deviation) (1 / (sqrt(2 * pi) * standard_deviation))... 
            * exp(-((x_vals - central_value).^2 / (2 * standard_deviation^2)));
    
        pdf_gauss = gaussian(ds_sigma);

        cdf_gauss = cumsum(pdf_gauss);
        cdf_gauss = 1-(cdf_gauss./max(cdf_gauss));
        P_all = [P_all, (cdf_gauss)'];

    end

    if log_axis
        intensity = log10(ppi_bins);
    else
        intensity = ppi_bins;
    end
    
    if log_x_axis
        imagesc(log10(vd), (intensity), P_all);
    else
        imagesc(vd, (intensity), P_all);
    end
    
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
            
            if log_x_axis
                vd_lines = log10(vd(~isnan(index_high)));
            else
                vd_lines = vd(~isnan(index_high))
            end

            if log_axis
                plot(vd_lines,... 
                    log10(ppi_bins(index_high(~isnan(index_high)))), ls{pp},... 
                    'Color', 'white', 'LineWidth', lw, 'HandleVisibility','off');
            else
                plot(vd_lines,... 
                    ppi_bins(index_high(~isnan(index_high))), ls{pp},... 
                    'Color', 'white', 'LineWidth', lw, 'HandleVisibility','off');
            end

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
    standard_ppi = [460, 136];
    px_text = {'iPhone 15', {'Samsung Neo','QLED 65"'}};
    offset = 1;
    x_pos = 1.7;
    if log_axis
        standard_ppi = log10(standard_ppi);
        offset = 0.13;
    end
    if log_x_axis
        x_pos = log10(x_pos);
    end
    for i = 1:length(standard_ppi)
        plot([-20, 20], [standard_ppi(i), standard_ppi(i)], '-r', 'HandleVisibility', 'off');
        text(x_pos, standard_ppi(i)+offset, px_text{i},... 
            'Color','red','FontSize',FontSize, 'VerticalAlignment', 'top');
    end
    
    
end
    
xlabel({'Viewing distance (m)'});
ylabel({'Pixel-per-inch (ppi) resolution'});
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

if log_axis
    ylim(log10([ppi_low, ppi_high]));
    ytick = ([1 3 10 30 100 300 1000]);
    set(gca, 'YTick', log10(ytick));
    set(gca, 'YTickLabel', ytick);
else
    ylim([ppi_low, ppi_high]);
    ytick = 0:100:1000;
    set(gca, 'YTick', ytick);
    set(gca, 'YTickLabel', ytick);
end

if log_x_axis
    xlim(log10([min(vd), max(vd)]));
    xtick = [0.1 0.2 0.4 1 2 4 10];
    set(gca, 'XTick', log10(xtick));
    set(gca, 'XTickLabel', xtick);
else
    xlim([min(vd), max(vd)]);
    xtick = 1:1:10;
    set(gca, 'XTick', xtick);
    set(gca, 'XTickLabel', xtick);
end

% Find all text objects in the figure
textObjects = findall(gcf, 'Type', 'text');
set(textObjects, 'FontSize', FontSize);

axesObjects = findall(gcf, 'Type', 'axes');
set(axesObjects, 'FontSize', FontSize);

if 0
    name = sprintf('plot_min_ppi_heatmap');
     print(fig, '-vector', name,'-dpdf');
    print(fig, '-vector', name,'-dpng', '-r600');
end
