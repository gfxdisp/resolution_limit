%% Plot resolution limit measured data with model fits

clear all; close all;
addpath('utils');

file_name = '../data/resolution_limit_data.csv';
ds = readtable(file_name, 'Delimiter',',');

text_file_name = '../data/text_resolution_limit_data.csv';
ds_text = readtable(text_file_name, 'Delimiter', ',');

contrast_vals = [0.913353, 0.134855, 0.487391];

central_val_model = 2; % 1: mean, 2: median

pars_all = [];

log_axis = true;
show_violin = true;
show_exponential_fit = true;

ci_factor = 95; % confidence interval

% figure,
gap = [0.01 ,0.01];
marg_h = [0.12, 0.2]; % [lower upper]
marg_w = [0.08, 0.01]; % [left right]
ha = tight_subplot(1,4, gap, marg_h, marg_w); % Create four subplots

% Adjust the positions
pos1 = get(ha(1), 'Position'); % Get position of the first subplot
pos2 = get(ha(2), 'Position'); % Get position of the second subplot

% Expand the first subplot to take more horizontal space
pos1(3) = pos1(3) + pos2(3);
set(ha(1), 'Position', pos1);

% Delete the second subplot
delete(ha(2));

% Reassign the remaining subplots
ha(2) = ha(3);
ha(3) = ha(4);

% Remove the now-unused fourth subplot handle
ha = ha(1:3);

set(gcf, 'defaultAxesFontName', 'Arial');
set(gcf, 'defaultTextFontName', 'Arial');
FontSize = 8;

fig = gcf;
fig_size = [7.5 3]*1;
fig.Units = 'inches';
fig.Position = [2, 2, fig_size];
fig.PaperUnits = 'inches';
fig.PaperSize = fig_size;

colors = [50, 50, 50;...
            255, 50, 0;...
            200, 100, 255]./255;

ls = {':', '-', '-.', '--'};
mk = {'<', 'd', 'o', 'v', 'p', 's', '^', 'h'};

hold on,

labels = {'Achromatic', 'Red-green', 'Yellow-violet'};
obs = unique(ds.observer);
obs_cols = lines(length(obs));

cols = [0.5 0.5 0.5;...
    1 0.5 0.7;...
    1 1 0.5];

ds_all_mean = [];
ds_all_median = [];
ds_all_pi = [];

for cc = 1:3
    axes(ha(cc));
    hh = [];
    hold on,

    alphaValue = 0.3;

    next_plt_idx = 2;
    ds_cc = ds((ds.color_direction == cc), :);
    
    ds_cc_ee_mean = [];
    ds_cc_ee_median = [];
    ds_cc_ee_pi = [];
    ds_cc_clean = [];

    if cc == 1
        eccs = [-20; -10; unique(ds_cc.eccentricity_deg)];
    else
        eccs = unique(ds_cc.eccentricity_deg);
    end

    for ee = 1:length(eccs)

        ds_cc_ee_temp = ds_cc(ds_cc.eccentricity_deg == eccs(ee), :);
        ds_cc_ee_temp.S = ones(height(ds_cc_ee_temp), 1)./contrast_vals(cc);
        pars_all(cc, 1) = 1./contrast_vals(cc);

        ds_cc_ee_temp.inv_cpd = 1./(ds_cc_ee_temp.threshold_ppd./2);

        switch eccs(ee)

            case -20 % white on black text
                allData = ds_text.threshold_ppd(ds_text.color_direction == -1);
                allDataScaled = ds_text.threshold(ds_text.color_direction == -1);
            case -10 % black on white text
                allData = ds_text.threshold_ppd(ds_text.color_direction == 1);
                allDataScaled = ds_text.threshold(ds_text.color_direction == 1);
            otherwise % Gabors
                allData = ds_cc.threshold_ppd(ds_cc.eccentricity_deg == eccs(ee));
                allDataScaled = ds_cc.threshold(ds_cc.eccentricity_deg == eccs(ee));
        end

        if ee == 1
            ds_cc_clean = ds_cc_ee_temp;
        else
            ds_cc_clean = [ds_cc_clean; ds_cc_ee_temp];
        end
        
        nBootstrap = 10000;  % Number of bootstrap samples
        bootstrapCenter = zeros(nBootstrap, 1);
        
        if 1 % confidence interval of mean/median
            for i = 1:nBootstrap
                % Sample with replacement from the aggregated data
                bootstrapSample = randsample(allData, length(allData), true);
                
                % Calculate the mean of this bootstrap sample
                bootstrapCenter_mean(i) = mean(bootstrapSample);
                bootstrapCenter_median(i) = median(bootstrapSample);
                
            end
            
            ciLower_mean = prctile(bootstrapCenter_mean, (100-ci_factor)/2);
            ciUpper_mean = prctile(bootstrapCenter_mean, ci_factor+(100-ci_factor)/2);
            
            ds_cc_ee_mean = [ds_cc_ee_mean; cc, eccs(ee), mean(bootstrapCenter_mean),... 
                ciLower_mean, ciUpper_mean];

            ciLower_median = prctile(bootstrapCenter_median, (100-ci_factor)/2);
            ciUpper_median = prctile(bootstrapCenter_median, ci_factor+(100-ci_factor)/2);
            
            ds_cc_ee_median = [ds_cc_ee_median; cc, eccs(ee), median(bootstrapCenter_median),... 
                ciLower_median, ciUpper_median];

        end

        if 1 % prediction interval from sample

            switch central_val_model
                case 1
                    % Calculate the mean of this bootstrap sample
                    bootstrap_samples = bootstrp(nBootstrap, @(x) [mean(x),... 
                        prctile(x,[(100-ci_factor)/2, ci_factor+(100-ci_factor)/2])],... 
                        allData);
                    plot_name = sprintf('Mean data (%d%% PI)', ci_factor);
                case 2
                    bootstrap_samples = bootstrp(nBootstrap, @(x) [median(x),... 
                        prctile(x,[(100-ci_factor)/2, ci_factor+(100-ci_factor)/2])],... 
                        allData);
                    plot_name = sprintf('Median data (%d%% PI)', ci_factor);
            end
            ds_cc_ee_pi = [ds_cc_ee_pi; cc, eccs(ee), mean(bootstrap_samples)];

        end

        if show_violin
%         hh(1) = Violin({allData}, 1);
            if log_axis
                vs = Violin({log10(allData)}, eccs(ee), 'Width', 2,...
                    'ViolinColor', {cols(cc, :)},...
                    'MedianMarkerSize', 50, 'BoxColor', [0 0 0],...
                    'ShowNotches', false,... 
                    'ShowMedian', true, 'ShowMean', true,...
                    'ShowWhiskers',true,...
                    'WhiskerData', log10(ds_cc_ee_median(ee, 4:5)),...
                    'ShowData', false,...
                    'QuartileStyle', 'shadow',...
                    'QuartileData', log10([ds_cc_ee_pi(ee, 4:5)]));
            else
                vs = Violin({(allData)}, eccs(ee), 'Width', 2,...
                    'ViolinColor', {cols(cc, :)},...
                    'MedianMarkerSize', 50, 'BoxColor', [0 0 0],...
                    'ShowNotches', false,... 
                    'ShowMedian', true, 'ShowMean', true,...
                    'ShowWhiskers',true,...
                    'WhiskerData', (ds_cc_ee_median(ee, 4:5)),...
                    'ShowData', false,...
                    'QuartileStyle', 'shadow',...
                    'QuartileData', ([ds_cc_ee_pi(ee, 4:5)]));
            end
        end

    end

    ds_all_mean = [ds_all_mean; ds_cc_ee_mean];
    ds_all_median = [ds_all_median; ds_cc_ee_median];
    ds_all_pi = [ds_all_pi; ds_cc_ee_pi];
   
  
    if show_exponential_fit
    
        eccentricity_plot = linspace(0, 50, 100);

        % Optimise for Watson's formula
        initialPars = double([1.85 -0.0525 0.393]);

        % Fit the non-linear model
        warning('off', 'stats:nlinfit:IllConditionedJacobian');

        % Specify options for fitnlm, changing the optimization algorithm
        options = statset('fitnlm');
        options.RobustWgtFun = 'cauchy';  % Example of robust fitting

        modelFunc = @(pars, X) 1.*power(Secc_watson2018(X(:,1), X(:,2),... 
            pars)./2, 1/3);
        nlModel = fitnlm([ds_cc_clean.S, ds_cc_clean.eccentricity_deg], -1.*ds_cc_clean.threshold, ...
            modelFunc, initialPars , 'Options', options);

        parsOptimized = nlModel.Coefficients.Estimate;
        pars_all(cc, 2:4) = parsOptimized';
        
        se = sqrt(diag(nlModel.CoefficientCovariance));
        delta = se * tinv(1-0.05/2,nlModel.DFE);
        
        parsOptimized_ci = [nlModel.Coefficients.Estimate-delta,... 
            nlModel.Coefficients.Estimate+delta];

        fittedCurve = Secc_watson2018(ds_cc_clean.S(1), eccentricity_plot,... 
            parsOptimized);
        
        % Get the covariance matrix of the estimated parameters
        covMatrix = nlModel.CoefficientCovariance;
        
        % Calculate the condition number
        condNumber = cond(covMatrix);
        % Large condition number suggests ill-conditioning

        % Check for normality of residuals
        [h, pValueNormality] = kstest(nlModel.Residuals.Standardized);
        
        if 1 % stats for latex table
            labels_short = {'Ach', 'RG', 'YV'};

            fprintf('Results for %s stimulus:\n', labels{cc});
            fprintf('Hypothesis test result (KS test): %d \n\n', h);
            fprintf('Residual normality test (KS test) p-value: %f\n\n', pValueNormality);

        end
    
if 0
        % Reporting results
        fprintf('Results for %s stimulus:\n', labels{cc});
        fprintf('Regression model: %f, %f, %f\n',... 
            nlModel.Coefficients.Estimate(1), nlModel.Coefficients.Estimate(2),... 
            nlModel.Coefficients.Estimate(3));
        fprintf('Standard errors: %f, %f, %f\n',... 
            nlModel.Coefficients.SE(1), nlModel.Coefficients.SE(2),... 
            nlModel.Coefficients.SE(3));
        fprintf('t-Statistics: %f, %f, %f\n',... 
            nlModel.Coefficients.tStat(1), nlModel.Coefficients.tStat(2),... 
            nlModel.Coefficients.tStat(3));
        fprintf('P-values: %f, %f, %f\n',... 
            nlModel.Coefficients.pValue(1), nlModel.Coefficients.pValue(2),...
            nlModel.Coefficients.pValue(3) );
        fprintf('R-squared: %f\n\n', nlModel.Rsquared.Ordinary);
        fprintf('Condition Number of the Covariance Matrix: %f\n\n', (condNumber));

        fprintf('Residual normality test (KS test) p-value: %d %f\n\n', h, pValueNormality);
        
end
        if 1 % per-observer parameter fitting
            parsOptimizedObs = [];
            idx = 1;
            for oo = 1:length(obs)
                ds_cc_obs = ds_cc_clean(strcmp(ds_cc_clean.observer, obs{oo}), :);
                if isempty(ds_cc_obs) || (height(ds_cc_obs) <= 2)
                    continue;
                end
                
                warning('off', 'stats:nlinfit:Overparameterized');

                modelFunc = @(pars, X) (Secc_watson2018(X(:,1), X(:,2), pars));
                nlModel = fitnlm([ds_cc_obs.S, ds_cc_obs.eccentricity_deg], ds_cc_obs.threshold_ppd, ...
                    modelFunc, initialPars);

                parsOptimizedObs(idx, :) = nlModel.Coefficients.Estimate';
                idx = idx+1;
            end

            numBootstraps = 1000; % Number of bootstrap samples

            % Sample with replacement from the original data
            bootstrapParams = parsOptimizedObs(randi(height(parsOptimizedObs),... 
                numBootstraps, 1), :);

            for i = 1:numBootstraps
                bootstrapPredictions(:, i) = Secc_watson2018(ds_cc_clean.S(1),... 
                    eccentricity_plot, bootstrapParams(i, :));
            end
            % Preallocate arrays to store the lower and upper bounds of the prediction intervals
            lowerBounds = zeros(length(eccentricity_plot), 1);
            upperBounds = zeros(length(eccentricity_plot), 1);

            % Calculate the prediction intervals for each point in 'eccentricity_plot'
            for i = 1:length(eccentricity_plot)
                % Sort predictions to find percentiles
                sortedPredictions = sort(bootstrapPredictions(i, :), 'ascend');
                
                % 2.5th and 97.5th percentiles for the 95% CI
                lowerBounds(i) = sortedPredictions(max(1, floor(numBootstraps... 
                    * (100-ci_factor)/200)));
                upperBounds(i) = sortedPredictions(min(numBootstraps, ceil(numBootstraps... 
                    * (ci_factor+(100-ci_factor)/2)/100)));
            end
            
            fittedCurve_lower = lowerBounds';
            fittedCurve_upper = upperBounds';
        else

            numBootstraps = 1000; % Number of bootstrap samples
            bootstrapParams = zeros(numBootstraps, length(parsOptimized));
            for i = 1:numBootstraps
                % Sample with replacement from the original data
                bootstrapSample = ds_cc_clean(randi(height(ds_cc_clean), height(ds_cc_clean), 1), :);
                
                % Fit the model to the bootstrap sample
                bootModel = fitnlm([bootstrapSample.S, bootstrapSample.eccentricity_deg],... 
                    bootstrapSample.threshold_ppd, modelFunc, initialPars);
                
                % Store the estimated parameters
                bootstrapParams(i, :) = bootModel.Coefficients.Estimate;

                bootstrapPredictions(:, i) = Secc_watson2018(ds_cc_clean.S(1),... 
                    eccentricity_plot, bootstrapParams(i, :));
            end
            
            % Preallocate arrays to store the lower and upper bounds of the prediction intervals
            lowerBounds = zeros(length(eccentricity_plot), 1);
            upperBounds = zeros(length(eccentricity_plot), 1);

            % Calculate the prediction intervals for each point in 'eccentricity_plot'
            for i = 1:length(eccentricity_plot)
                % Sort predictions to find percentiles
                sortedPredictions = sort(bootstrapPredictions(i, :), 'ascend');
                
                % 2.5th and 97.5th percentiles for the 95% CI
                lowerBounds(i) = sortedPredictions(max(1, floor(numBootstraps... 
                    * (100-ci_factor)/200)));
                upperBounds(i) = sortedPredictions(min(numBootstraps, ceil(numBootstraps... 
                    * (ci_factor+(100-ci_factor)/2)/100)));
            end
            
            fittedCurve_lower = lowerBounds';
            fittedCurve_upper = upperBounds';
        end
            
        warning('on', 'all');

        if log_axis
            plot(eccentricity_plot, log10(fittedCurve), '--k',...
                'DisplayName', 'Exponential fits');
        else
            plot(eccentricity_plot, (fittedCurve), '--k',...
                'DisplayName', 'Exponential fits');
        end

        if log_axis
            text(ds_cc_ee_pi(:, 2)+1, log10(ds_cc_ee_pi(:, 5))+0.07,... 
                num2str(round(ds_cc_ee_pi(:, 3), 0)));
        else
            text(ds_cc_ee_pi(:, 2)+1, ds_cc_ee_pi(:, 5)+7,... 
                num2str(round(ds_cc_ee_pi(:, 3), 0)));
        end

    else
        hh(next_plt_idx) = errorbar(ds_cc_ee_pi(:, 1), ds_cc_ee_pi(:, 2),... 
            ds_cc_ee_pi(:,3)-ds_cc_ee_pi(:,2),... % lower ci
            ds_cc_ee_pi(:, 4)-ds_cc_ee_pi(:,2),... % upper ci
            'Color', 'k', 'Marker', 'o', 'MarkerEdgeColor', 'k',... 
            'MarkerFaceColor', cols(cc, :), 'LineStyle', '-',...
            'MarkerSize', 10, 'CapSize', 20, 'LineWidth', 1,...
            'DisplayName', plot_name);
    end

    title(labels{cc}, 'Fontweight', 'Normal');
    
    if cc == 2
        xlabel('Retinal eccentricity (visual degrees)')
    end

    if cc == 1
        ylabel({'Resolution limit', '(pixels per degree)'});
    end

    if cc == 3
        hh(1) = plot(NaN, NaN, '-', 'Color', [cols(cc, :), 0.3], ...
            'LineWidth', 10,... 
            'DisplayName', 'Probability density estimate');
        hh(2) = plot(NaN, NaN, '-', 'Color', [cols(cc, :), 1], ...
            'LineWidth', 10,... 
            'DisplayName', sprintf('%dth percentile', ci_factor));
        hh(3) = errorbar(50, 1, 1, 'Marker', 'o', 'Color', 'k',... 
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'white',... 
            'LineStyle', 'none', 'MarkerSize', 50, 'Capsize', 50,... 
            'LineWidth', 1,...
            'DisplayName', sprintf('Median (%d%% CI)', ci_factor));
        hh(4) = plot(NaN, NaN, '-', 'Color', 'k', 'LineWidth', 2,...
            'DisplayName', 'Mean');
        hh(5) = plot(NaN, NaN, '--', 'Color', 'k', 'LineWidth', 1,...
            'DisplayName', 'Model fit (Watson 2018)');
        legend boxoff
        legend(hh, 'NumColumns', 3,... 
            'Position', [0.1913    0.8824    0.6279    0.1117]);
    end

    grid on
    if log_axis
        if cc == 1
            xlim([-24, 25]);
        else
            xlim([-4, 25]);
        end
        
        ylim(log10([1, 200]));

        ytick = [0.3 1 3 10 30 100 300];
        set(gca, 'YTick', log10(ytick));
    else
        if cc == 1
            xlim([-24, 25]);
        else
            xlim([-4, 25]);
        end

        ylim(([-5, 130]));
        
        ytick = [0:20:150];
        set(gca, 'YTick', (ytick));
    end
       
    if cc == 1
        xtick = [eccs(3:end)];
    else
        xtick = [eccs];
    end
    set(gca, 'XTick', xtick);
    set(gca, 'XTickLabel', xtick);
    
    
    if cc == 1
        % Define the positions for the custom labels
        custom_xtick_pos = [-20, -10]; % Example positions, adjust as needed
        
        hold on,
        
        if log_axis 
            text_offset = -0.1;
        else
            text_offset = -10;
        end
        % Add the custom labels using text function
        % Custom label: white text on black box
        text(custom_xtick_pos(1), text_offset, ' Text', 'Color', 'white', 'BackgroundColor', 'black', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontWeight', 'bold', ...
            'FontSize', 8, 'Margin', 2);
        
        % Custom label: black text on white box with black outline
        % Position and size of the box
        box_width = 10; % Adjust box width as needed
        box_height = 0.1; % Adjust box height as needed
        box_pos = [custom_xtick_pos(2) - box_width/2, -0.1 - box_height/2, box_width, box_height];
                
        % Add the text inside the box
        text(custom_xtick_pos(2), text_offset, ' Text', 'Color', 'black', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontWeight', 'bold', ...
            'FontSize', 8);  
    end
    
    if cc == 1
        set(gca, 'YTickLabel', ytick);
    end
end

% Find all text objects in the figure
textObjects = findall(gcf, 'Type', 'text');
set(textObjects, 'FontSize', FontSize);

axesObjects = findall(gcf, 'Type', 'axes');
set(axesObjects, 'FontSize', FontSize);

legendHandle = legend;  % If a legend exists, get the handle
if ~isempty(legendHandle)
    set(legendHandle, 'FontSize', FontSize);
end



