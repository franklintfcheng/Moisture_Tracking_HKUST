%%
function [str_scaleavg, scaleavg_all, scaleavg_sig_all, X_bandpass, X_reconstruct, X_reconstruct_s1s2] = plot_mywavelet(varargin)

%/ Author: Fandy
%/ Date of creation: 8 Jun 2021

%/ create a set of valid parameters and their default value
pnames = {'X', 'X_dates', 'str_Xname', 'detrend_mode',  'deseasonalize_mode',  'contf_levels', 'cmap', 'unit', 's1',  's2',...
          'YDir', 'fontsize', 'create_fig', 'draw_fig', 'draw_scaleavg', 'draw_cbar_only', 'fig_width_deplete', 'titlename', 'plotting_folder', 'savefig', 'save_scaleavg_path', 'save_bandpass', 'save_reconstruct'};  
      
dflts  = { [],       [],          [],              1,                     1,              [],     [],     [],   [],    [],...
          'normal', 13,            1,          1,                  1,                 0,                   0,                [],               [],         0,                  [],                     0,                   0};

[X, X_dates, str_Xname, detrend_mode, deseasonalize_mode, contf_levels, cmap, unit, s1, s2,...
 YDir, fontsize, create_fig, draw_fig, draw_scaleavg, draw_cbar_only, fig_width_deplete, titlename, plotting_folder, savefig, save_scaleavg_path, save_bandpass, save_reconstruct]...
            = internal.stats.parseArgs(pnames, dflts, varargin{:});        %/ parse function arguments
        
if size(X, 1) ~= 1 && size(X, 2) ~= 1      error('Input X must be an 1D array!');      end
if savefig && isempty(plotting_folder)     error('Specify the plotting_folder!');      end

%/ Setting based on date format
if numel(num2str(X_dates(1))) == 8                                         %/ assume yyyymmdd
    mod_divider      = 1e4; 
    first_date       = 101;
    CONV_to_yr       = 1/365.25;
    CONV_restore     = 30;        
    CONV_to_sec      = 24*3600;
    slct_period_ind  = 33:92;                                              %/ not showing the very short/long periods
    str_timescale    = 'daily';
    timescale_unit   = 'day';

    
elseif  numel(num2str(X_dates(1))) == 6                                    %/ assume yyyymm

    mod_divider      = 1e2; 
    first_date       = 1;
    CONV_to_yr       = 1/12;
    CONV_restore     = 1;    
    CONV_to_sec      = 30*24*3600;
    slct_period_ind  = [];                                                 %/ not showing the very short/long periods
    str_timescale    = 'mthly';
    timescale_unit   = 'month';
       
else
    error('Check the format of the input date! It has to be either in yyyymmdd or yyyymm!')
end

%/ output the scale avg data with default time scales if s1 s2 are empty
if isempty(s1) || isempty(s2)                                          
    %/ s1 and s2 -> must be in the same time unit of the input data.
    s1 = [1*CONV_restore*CONV_to_yr, 3*CONV_restore*CONV_to_yr,  1,  3,   7]/CONV_to_yr;
    s2 = [3*CONV_restore*CONV_to_yr,                         1,  3,  7,  10]/CONV_to_yr;
end
    

str_detrend = []; str_deseason = [];
if detrend_mode             X = detrend(X);                str_detrend  = ' detrend';      end         %/ remove long-term impact (e.g., global warming)
if deseasonalize_mode       X = deseasonalize(X, X_dates); str_deseason = ' deseason';     end         %/ remove seasonality


% warning('testing bandpass')
% fs    = 1/CONV_to_sec;                     %/ sample rate (Hz)
% fpass = 1./([s2(4) s1(4)]*CONV_to_sec);    %/ passband frequency (Hz)  
% X = bandpass(X, fpass, fs, 'Steepness', 1);


%/ Data processing
% normalize by standard deviation (not necessary, but makes it easier
% to compare with plot on Interactive Wavelet page, at
% "http://paos.colorado.edu/research/wavelets/plot/"
X_variance        = std(X)^2;
X_norm            = (X - mean(X))/sqrt(X_variance); 

X_years           = unique(floor(X_dates/mod_divider));
time              = 1:length(X_dates);
ind_1stDateOfYear = find(mod(X_dates,mod_divider) == first_date);          

FigName            = sprintf('wavelet %s %s%s%s', str_Xname, str_timescale, str_detrend, str_deseason);
FigName_underscore = strrep(FigName, ' ', '_');
        
if draw_cbar_only
    fprintf('*** Drawing cbar only... *** \n');
    figure
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf, 'color','w');

    caxis([min(contf_levels) max(contf_levels)]);
    colormap(gca, cmap)

    cbar_location = 'southoutside';
    cb = colorbar(cbar_location);
    cbar_interval   = 2;
    cbar_YTick      = contf_levels(2:cbar_interval:end-1);
    cbar_YTickLabel = cbar_YTick;

    axis off
    if ismember(cbar_location, {'southoutside'})
        set(cb, 'position', [.1 .3 .4 .05]);    %/ [xposition yposition width height]
        set(cb, 'YAxisLocation','right')
        cb_fontsize = 20;
    elseif ismember(cbar_location, {'eastoutside'})
        set(cb, 'position', [.1 .3 .04 .4]);    %/ [xposition yposition width height]
        set(cb, 'YAxisLocation','right')
        cb_fontsize = 30;
    end
    set(cb, 'YTick', cbar_YTick, 'YTickLabel', cbar_YTickLabel, 'Fontsize', cb_fontsize) %/ cbar Ytick for diverging colormap
    set(get(cb,'Title'),'String', '', 'Fontsize', cb_fontsize)
    drawnow; pause(0.05);
    
    if savefig
%         export_fig(char(strcat(plotting_folder, FigName_underscore,'_cbar.png')),'-r300','-png','-opengl', '-nocrop'); %'-transparent');
        export_fig(char(strcat(plotting_folder, FigName_underscore,'_cbar.pdf')),'-pdf','-painters', '-c[inf, inf, inf, inf]', '-nocrop'); %'-transparent');
    end

    return
end


if create_fig && draw_fig
    figure
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf, 'color','w');
end

axesWidth  = 1.5;
time_Xtick = time(ind_1stDateOfYear);
Xtick_intvl = 10;

%/ Wavelet Parameterization
n = length(time);
xlimit = [1 n]; %plotting range
dt = 1;
pad = 1;      % pad the time series with zeroes (recommended)
dj = 0.125;   % this will do 8 sub-octaves per octave
s0 = 2*dt;    % this says start at a scale of 2dt
j1 = -1;      % use default
%j1 = 7/dj;    % this says do 7 powers-of-two with dj sub-octaves each
lag1 = 0.72;  % lag-1 autocorrelation for red noise background
mother = 'Morlet';

% Wavelet transform:
[wave,period,scale,coi] = wavelet(X_norm,dt,pad,dj,s0,j1,mother);
power = (abs(wave)).^2;          % compute wavelet power spectrum


% Significance levels: (variance=1 for the normalized data)
[signif,fft_theor] = wave_signif(1.0,dt,scale,0,lag1,-1,-1,mother);
sig95 = (signif')*(ones(1,n));  % expand signif --> (J+1)x(N) array
sig95 = power ./ sig95;         % where ratio > 1, power is significant


% Global wavelet spectrum & significance levels:
global_ws = X_variance*(sum(power')/n);   % time-average over all times
dof = n - scale;  % the -scale corrects for padding at edges
global_signif = wave_signif(X_variance,dt,scale,1,lag1,-1,dof,mother);


%--- Plot time series
%subplot('position',[0.1 0.75 0.65 0.2])
%     plot(time,X_daily_norm)
%     set(gca,'XLim',xlim(:))
%     xlabel('Time (years)')
%     ylabel(X_name)
%     title(X_name)
%     set(gca,'XLim',xlim(:), ...
%         'XTick',time_Xtick(1:3:end), ...
%         'XTickLabel',X_years(1:3:end),'FontSize', fontsize);
%     hold off

%--- Contour plot wavelet power spectrum
if draw_fig
    sp_x      = 0.13;
    sp_y      = 0.4;
    sp_width  = 0.6*fig_width_deplete;
    sp_height = 0.44;
    
    gsp_width = 0.1*fig_width_deplete;
    subplot('position',[sp_x sp_y sp_width sp_height])   %/ x  y  width height

    if isempty(slct_period_ind)
        
        [B, ~] = max(coi);
        slct_period_ind = find(period <= B);
%         slct_period_ind = 1:length(period);   
    end

    
    %/ original data along y-axis is in log(period), where period is in month.
    %/ (1) subset the period and convert its unit to years
    %/ (2) create the preferred yticklabel (in years)
    %/ (3) convert the unit of yticklabel back to month and save as Yticks

    period_in_yr = period(slct_period_ind)*CONV_to_yr;                              %/ (1)
    Yticklabel   = 2.^(fix(min(log2(period_in_yr))):fix(max(log2(period_in_yr))));  %/ (2)
    Yticks = Yticklabel/CONV_to_yr;                                                 %/ (3)
    
%     Yticks = 2.^(fix(log2(min(period(slct_period_ind)))):fix(log2(max(period(slct_period_ind))))); % convert log2(period) back to period in y-axis; 
%     Yticklabel = round(Yticks/CONV_from_yr, 1);        %/ convert from days to years

    warning('Note that the both the y-axis and the data are on log2. That''s why we have -ve values in power spectrum.');
    
%     if isempty(contf_levels)
%         contf_range     = prctile(reshape(log2(power), [], 1), [0 100]);             %/ set upperbound to be 90th percentile
%         contf_levels = linspace(contf_range(1)-2, contf_range(end)+2, 11);   
%         cmap = brewermap(length(contf_levels)-1, '*RdBu');
%     end
%     
    contourf(time,log2(period),log2(power), [-99999 contf_levels], 'LineColor','none');
%     xlabel('Time (years)')
    xlabel('')
    ylabel('Period (years)')
    colormap(cmap);
    caxis([min(contf_levels) max(contf_levels)]);
    shading interp;
    % colorbar('southoutside');
    % title(strcat('(a) ',X_name,' Wavelet Power Spectrum'),'FontSize',16)
    set(gca,'XLim',xlimit(:), ...
        'XTick',time_Xtick(5:Xtick_intvl:end), ...
        'XTickLabel',X_years(5:Xtick_intvl:end),'FontSize',fontsize);

    set(gca,'YLim',log2([min(period(slct_period_ind)),max(period(slct_period_ind))]), ...
        'YDir', YDir, ...
        'YTick',log2(Yticks), ...
        'YTickLabel',Yticklabel,'FontSize',fontsize);
    
    % 95% significance contour, levels at -99 (fake) and 1 (95% signif)
    hold on
    contour(time,log2(period),sig95,[-99,1], 'color', [51 204 0]./255, 'Linewidth',2);
    hold on
    % cone-of-influence, anything "below" is dubious
    plot(time,log2(coi),'--','Linewidth',3,'Color',[1 1 1]);
%     plot(time,log2(coi),'-','Linewidth',2,'Color',[51 204 0]./255);
    set(gca,'linewidth',axesWidth)
    hold off

    th = title(titlename);
%     titlePos = get(th, 'position');
%     titlePos(2) = titlePos(2) + 0.07;   %/ change y position of title
%     set(th, 'position', titlePos)
    set(th, 'fontsize', fontsize+15, 'fontweight', 'normal');
    drawnow; pause(0.05);


    %--- Plot global wavelet spectrum
    subplot('position',[sp_x+sp_width+0.01 sp_y gsp_width sp_height])   %/ x  y  width height
    
    ind_nonsig = find(global_ws < global_signif);
    global_sig = global_ws;
    global_sig(ind_nonsig) = NaN;
    clear ind_sig
    clear ind_nonsig

%   %/ my way to plot global spectrum
%     plot(global_ws,log2(period),'--','LineWidth', 2, 'Color', [0.2 0.2 0.2])
%     hold on
%     plot(global_sig,log2(period),'-','LineWidth', 2.5, 'Color', [0.2 0.2 0.2])
    
    %/ traditional way to plot global spectrum
    plot(global_ws, log2(period), '-', 'LineWidth', 2.5, 'Color', [0.2 0.2 0.2])
    hold on
    plot(global_signif, log2(period), '--', 'LineWidth', 1.8, 'Color','r');

    hold off
    xlabel('')
%     xlabel(sprintf('Power %s',unit))

    set(gca,'YLim',log2([min(period(slct_period_ind)),max(period(slct_period_ind))]), ...
        'YDir', YDir, ...
        'YTick', [], ...
        'YTickLabel', {[]},'FontSize',fontsize);

    if draw_scaleavg
        xaxisLocation = 'top';
    else
        xaxisLocation = 'bottom';
    end
    set(gca,'XLim',[0,1.25*max(global_ws)], 'xaxisLocation', xaxisLocation, 'FontSize', fontsize)
    set(gca,'linewidth',axesWidth)

end



%--- Plot aa-bb yr scale-average time series
if draw_fig && draw_scaleavg
    subplot('position',[sp_x 0.2 sp_width 0.17])
end

str_scaleavg = cell(length(s1),1);
for k = 1:length(s1)
    if s1(k) < 1/CONV_to_yr                                        %/ check the lower-end of the time scale only.
        
        if 1/CONV_to_yr == 365.25
            if s2(k) == 1/CONV_to_yr
                str_scaleavg{k} = sprintf('%d-%d-%s', s1(k)/30, 12, 'month'); %/ just for a consistent string name.
            else
                str_scaleavg{k} = sprintf('%d-%d-%s', s1(k)/30, s2(k)/30, 'month'); %/ just for a consistent string name.
            end
            
        elseif 1/CONV_to_yr == 12
            str_scaleavg{k} = sprintf('%d-%d-%s', s1(k), s2(k), timescale_unit);
        end
        
    else
        str_scaleavg{k} = sprintf('%d-%d-%s', s1(k)*CONV_to_yr, s2(k)*CONV_to_yr, 'year');
    end
end

C = brewermap(length(s1)+1,'Set1');
ind_yellow = findismember_loop(C(:,2), [1 1 0.2]);
C(ind_yellow,:) = [];   %/ remove bridght yellow;


%/ Save scale-avg / bandpass data
scaleavg_all     = nan(length(s1), length(time));
scaleavg_sig_all = nan(length(s1), length(time));

cnt = 1; scaleavg_plot_sig95 = []; str_scaleavg_cnt = {}; X_bandpass = []; X_reconstruct = [];  X_reconstruct_s1s2= [];
for k = 1:length(s1)
    avg = find((scale >= s1(k)) & (scale < s2(k)));
    Cdelta = 0.776;                     % this is for the MORLET wavelet
    scale_avg = (scale')*(ones(1,n));   % expand scale --> (J+1)x(N) array
    scale_avg = power ./ scale_avg;     % [Eqn(24)]
    scale_avg = X_variance*dj*dt/Cdelta*sum(scale_avg(avg,:),1);   % [Eqn(24)]

    scaleavg_signif = wave_signif(X_variance,dt,scale,2,lag1,-1,[s1(k), s2(k)],mother);

    ind_nonsig = find(scale_avg < scaleavg_signif);
    scale_avg_sig = scale_avg;
    scale_avg_sig(ind_nonsig) = NaN;
    
    %/ store scale-avg data
    scaleavg_all(k,:) = scale_avg;
    scaleavg_sig_all(k,:) = scale_avg_sig;
   
    if all(isnan(scale_avg_sig))       continue;        end      %/ skip figure plotting and data saving if the scaleavg is entirely insig.
    
    if draw_fig && draw_scaleavg
        str_scaleavg_cnt{cnt} = str_scaleavg{k};                 %/ legend string
    
        %/ take log2 for a better visualization
        plot(time, log2(scale_avg), '--', 'LineWidth', 1.3, 'Color',C(k,:));  %/ BUG: dash lines won't show if linewidth is too thick! ....
        hold on;

        %/ scale avg line plot (sig)
        scaleavg_plot_sig95(cnt)=plot(time, log2(scale_avg_sig), 'LineWidth', 2.5, 'Color',C(k,:));
        % scaleavg_plot_sig95(k) = plot([1 length(time)],scaleavg_signif+[0,0],'--','LineWidth',1.1);
        hold on
    end

    %/ save (sig) scale avg data for later correlation with CIs!
    if ~isempty(save_scaleavg_path) 
        scaleavg_filename = string(strcat(save_scaleavg_path, sprintf('scaleavg_%s_%s.mat', titlename, str_scaleavg{k})));
        scaleavg_filename = strrep(scaleavg_filename, ' ', '_');

        scaleavg_data = [X_dates, scale_avg'];

        fprintf('*** Saving data: %s *** \n', scaleavg_filename)
        save(scaleavg_filename, 'scaleavg_data', '-v7.3');
    end
    
    %/ bandpass filter on the sig. time scale of the input signal
    if save_bandpass
        
        fs    = 1/CONV_to_sec;                         %/ sample rate (Hz)
        fpass = 1./([s2(k) s1(k)]*CONV_to_sec);        %/ passband frequency (Hz)  NOTE: s1 and s2 are in the same time unit of input. 
        
        
        X_norm_bandpass = bandpass(X_norm, fpass, fs); %/ since we do normalization on X before inputing it into wavelet!! So make it consistent.
        bandpass_data = [X_dates, X_norm_bandpass'];

        bandpass_filename = string(strcat(save_scaleavg_path, sprintf('bandpass_%s_%s.mat', titlename, str_scaleavg{k})));
        bandpass_filename = strrep(bandpass_filename, ' ', '_');

        fprintf('*** Saving bandpass data: %s *** \n', bandpass_filename)
        save(bandpass_filename, 'bandpass_data', '-v7.3');
        
        if k == 4
            X_bandpass = X_norm_bandpass * sqrt(X_variance) + mean(X);  %/ only for verification (3-7-yr)
        end
    end
    
    if save_reconstruct
   
        real_W    = real(wave)';
        phi_log_0 = pi^(-1/4);
        C_delta   = 0.776;

        %/ full wavelet reconstruction
        X_norm_reconstruct = dj*dt^(1/2)/C_delta/phi_log_0*sum(real_W.*scale.^(-1/2), 2);     %/ Eq.11 in Torrence and Compo 1998
        X_reconstruct      = X_norm_reconstruct * sqrt(X_variance) + mean(X);      
 
        RMSE = sqrt(mean((X_reconstruct - X').^2));
        fprintf('*** RMSE for the reconstructed data = %.4f ***\n', RMSE);
        
        %/ wavelet-based bandpass filter
        ind = find((scale >= s1(k)) & (scale < s2(k)));
        
        X_norm_reconstruct_s1s2 = dj*dt^(1/2)/C_delta/phi_log_0*sum(real_W(:,ind)./scale(ind).^(1/2), 2);
        X_reconstruct_s1s2      = X_norm_reconstruct_s1s2 * sqrt(X_variance) + mean(X);      
 
        reconstruct_data = [X_dates, X_reconstruct_s1s2];
        
        reconstruct_filename = string(strcat(save_scaleavg_path, sprintf('reconstruct_%s_%s.mat', titlename, str_scaleavg{k})));
        reconstruct_filename = strrep(reconstruct_filename, ' ', '_');

        fprintf('*** Saving reconstruction data: %s *** \n', reconstruct_filename)
        save(reconstruct_filename, 'reconstruct_data', '-v7.3');
    end
    
    cnt = cnt + 1;
end

if draw_fig && draw_scaleavg
    xlabel('Time (years)')
    ylabel('log2(Avg variance)')
%     ylabel(sprintf('Avg variance %s',unit))
    set(gca,'XLim',xlimit(:), ...
            'XTick',time_Xtick(5:Xtick_intvl:end), ...
            'XTickLabel',X_years(5:Xtick_intvl:end),'FontSize',fontsize);
    set(gca,'linewidth',axesWidth)
    lgd = legend(scaleavg_plot_sig95(:), str_scaleavg_cnt{:}, 'location', 'eastoutside');
    lgd.FontSize = 14;

    %------ legend at the bottom right corner ------%
    set(lgd, 'Position', [sp_x+sp_width+0.01 0.2 gsp_width 0.17]);
end

if draw_fig && savefig
    export_fig(char(strcat(plotting_folder, FigName_underscore,'.png')),'-r300','-png','-opengl', '-nocrop'); %'-transparent');
%     export_fig(char(strcat(plotting_folder, FigName_underscore,'.pdf')),'-pdf','-painters', '-c[inf, inf, inf, inf]', '-nocrop'); %'-transparent');
end


end