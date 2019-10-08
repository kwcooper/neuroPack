function view_all_lfp()

% TODO:
%   Abstract this code into a working function!

    num_snaps = 5;     % set the number of lfp snapshots that you'd like
    recLength = 1000;  % set how long you'd like the lfp clip to be
    save_figs = 1;
    
    for offset = linspace(1,size(lfpStruct.data,2)-recLength-1,num_snaps)
        offset = round(offset);  % round to get a whole number index

        % show all of the LFP on one big colored plot.
        figure('Position', [1933, -72, 1005, 1299]);
        num_plots = 2;
        [num_ch, num_tp] = size(lfpStruct.data);
        chansPerPlot = num_ch / num_plots;
        inc = 0;
        start_ch = 1; end_ch = chansPerPlot;

        % Compute the amount of LFP that you'd like to show and where
        %offset = 100000;
        startInd = 1+offset; endInd = recLength+offset;
        timeAx = linspace(startInd/fs, endInd/fs, endInd-startInd+1);
        sLFP = spreadLFP(lfpStruct.data(:, startInd:endInd));

        % edit the colors for the channels
        % initialize a grey matrix
        colors = repmat([.5,.5,.5],64,1);
        % figure; imshow(colors);
        % Color the matrix based on the channel map
        names = fieldnames(ratMetaData.(ratNames{ratInd}).chMap);
        goodChan = false;
        % generate some fancy colors
        chanCol = hsv(length(names));
        for nIdx = 1:length(names)
            cIdx = ratMetaData.(ratNames{ratInd}).chMap.(names{nIdx});
            if goodChan; newCol = [1 0 0];
            else; newCol = chanCol(nIdx,:); end
            colors(cIdx,:) = repmat(newCol,length(cIdx),1);
        end
        % figure; imshow(colors)
        %colors = flipud(colors);

        chanCol = hsv(length(names));
        regionsLeg = '\fontsize{10}';
        for nIdx = 1:length(names)
            cIdx = ratMetaData.(ratNames{ratInd}).chMap.(names{nIdx});
            newCol = chanCol(nIdx,:);
            textBuff = ['\color[rgb]{', num2str(newCol(1)),',',num2str(newCol(2)),',',num2str(newCol(3)),'}   ', names{nIdx}];
            regionsLeg = [regionsLeg, textBuff];
        end

        suptitle({[ratMetaData.(ratNames{ratInd}).ratName, ' Ephys'], regionsLeg});
        for pltInd = 1:num_plots
            subplot(num_plots/2,num_plots,pltInd);
            h = plot(timeAx, sLFP(start_ch:end_ch,:)');
            xlabel({'Seconds',['Chs: ', num2str(start_ch), ' to ', num2str(end_ch)]});
            set(gca, 'YTick', []);
            %colors = hsv(chansPerPlot); % basic rainbox plot
            set(h, {'color'}, num2cell(colors(start_ch:end_ch,:), 2));

            inc = inc + chansPerPlot;
            start_ch = start_ch + inc; end_ch = end_ch + inc;
        end

        if save_figs % Save the figure;
            %imagePath = 'C:\Users\FortinLab\Documents\Analysis\openEphysSandbox\figs\SM04\';
            imagePath = ratMetaData.(ratNames{ratInd}).imgDir;
            imgName = [ratMetaData.(ratNames{ratInd}).ratName, '_Ephys_coloredCh_', num2str(round(timeAx(1))),'-',num2str(round(timeAx(end))),'s.png'];
            saveas(gcf,[imagePath,imgName]);
        end

    end
    if save_figs; fprintf('Images Saved!\n'); end
    
    