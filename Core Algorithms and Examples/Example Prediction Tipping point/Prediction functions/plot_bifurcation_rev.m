%% Bifurcation diagram under this parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Author: Anastasia Ilina / Junheng Li

% Performs plotting of the bifurcation model fit to the state variable data
% as well as the bifurcation state diagram, demonstrating the areas with 1,
% 2 and 3 solutions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SVariableCriticalZone, c_3sol_output]  = plot_bifurcation_rev(dd, xx_smooth, tvec, participantName, params_opt1, ifreturn3sol)

SVariableCriticalZone = [];
c_3sol_output = [];

handle = figure; % Open a new figure and store the handle in h

hold on

plot(tvec,xx_smooth,'LineWidth',2, 'Color','k')
if ~isempty(dd(:,1))
    if length(tvec) == length(dd(:,1))
        plot(tvec,dd(:,1),'LineWidth',4,'Color','blue')
    end
else 
    disp('Fitted Bifurcation is empty!')
end
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
set(gca, 'Box', 'off')
set(gca, 'FontName', 'Helvetica')
ylabel('System state')
xlabel('Time (min)')
ax = gca;
%ylim([0,8])
line([0,0],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
xticks(-30:5:10)
xlim([-30,10])


if ~ isempty(dd(:,2))
    c = dd(:,2);

    r = params_opt1(1);
    K = params_opt1(2);
    h= params_opt1(3);
    m = params_opt1(4);

    x_theo = zeros(length(c));
    c_3sol = [];

    handle = figure;
    hold on

    syms x

    for idxc = 1:length(c)

        % Roots of the theoretical equation when dx/dt = 0;
        %    x = roots([1/10, -1, (c(idxc)+1/10), -1,0]);
        eqn = r*x*(1-x/K) - c(idxc)*x^2 / (x^2+h^2);
        xsol_all = solve(eqn==0);
        xsol_all = double(xsol_all);

        xsol = [];
        for jjj = 1:length(xsol_all)
            if isreal(xsol_all(jjj))
                if abs(xsol_all(jjj)) ~= 0
                    xsol = [xsol,xsol_all(jjj)];
                end
            end
        end
        % Find 3 solution points
        if length(xsol) == 3
            c_3sol = [c_3sol,c(idxc)];
        end

        scatter(c(idxc)*ones(1,length(xsol)),xsol,2,'k')

    end

    xlabel('Control variable')
    ylabel('System state')
    set(gca,'FontSize', 16)
    set(gca,'TickDir','out')
    set(gca,'ticklength',3*get(gca,'ticklength'))
    set(gca,'lineWidth',2)


    %% Line plots for bifurcation theoretical model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % Try to do line plots for easy later processing
    % Find the range where there will be 3 solutions in system

    % c_min = 1.7880;
    % c_max = 2.6040;
    if isempty(c_3sol)
        disp('No bifurcation point identified')
    else


    %% plot the bifurсation diagram and fitted bifurcation in one plot together
    handle = figure;
    set(handle, 'Position', [100, 100, 500, 700]); % [left, bottom, width, height]
    % Create a subplot
    ax1 = subplot(2, 1, 1);
    hold on
    
    plot(tvec,xx_smooth,'LineWidth',2, 'Color','k')
    if ~isempty(dd(:,1))
        if length(tvec) == length(dd(:,1))
            plot(tvec,dd(:,1),'LineWidth',4,'Color','blue')
        end
    else 
        disp('Fitted Bifurcation is empty!')
    end
    set(gca,'FontSize', 16)
    set(gca,'TickDir','out')
    set(gca,'ticklength',2*get(gca,'ticklength'))
    set(gca,'lineWidth',2)
    set(gca, 'Box', 'off')
    set(gca, 'FontName', 'Helvetica')
    ylabel('System state')
    xlabel('Time (min)')
    ax = gca;
    %ylim([0,8])
    line([0,0],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
    xlim([-30,10])

    % Construct the title with line breaks
    graphTitle = {'Fitted Bifurcation mode', ...
        ['Participant ', num2str(participantName)]}; %, ... %' Mon. decr. model ', num2str(is_monotonic)],
        %['R^2 = ', num2str(predictionSummaryTableRow.R_squared), ', RMSE = ', num2str(predictionSummaryTableRow.Min_RMSE)]};

    subplot(2, 1, 2)

     %% Line plots for bifurcation theoretical model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % Try to do line plots for easy later processing
    % Find the range where there will be 3 solutions in system

    % c_min = 1.7880;
    % c_max = 2.6040;
    if isempty(c_3sol)
        disp('No bifurcation point identified')
    else


        idx1 = find(c == c_3sol(1));
        idx2 = find(c == c_3sol(end));

        hold on
        % Plot in different ranges
        cvec = [];
        xvec = [];
        for idxc = 1:idx1-1

            eqn = r*x*(1-x/K) - c(idxc)*x^2 / (x^2+h^2);
            xsol_all = solve(eqn==0);
            xsol_all = double(xsol_all);
            xsol = [];
            for jjj = 1:length(xsol_all)
                if isreal(xsol_all(jjj))
                    if abs(xsol_all(jjj)) ~= 0
                        xsol = [xsol,xsol_all(jjj)];
                    end
                end
            end
            if length(xsol) ~= 1
                warning('Check code')
            end
            cvec = [cvec,c(idxc)];
            xvec = [xvec,xsol];

        end
        %plot(tvec, xvec, 'k', 'lineWidth', 2)
        plot(cvec, xvec,'k','lineWidth',2);

        cvec = [];
        xvec = [];
        for idxc = idx1:idx2

            eqn = r*x*(1-x/K) - c(idxc)*x^2 / (x^2+h^2);
            xsol_all = solve(eqn==0);
            xsol_all = double(xsol_all);
            xsol = [];
            for jjj = 1:length(xsol_all)
                if isreal(xsol_all(jjj))
                    if abs(xsol_all(jjj)) ~= 0
                        xsol = [xsol,xsol_all(jjj)];
                    end
                end
            end
            if length(xsol) ~= 3
                warning('Check code')
            end
            cvec = [cvec,c(idxc)];
            xvec = [xvec;xsol];

        end
        %plot(tvec, xvec(:,1),'k','lineWidth',2);
        %plot(tvec, xvec(:,2),'k--','lineWidth',2);
        %plot(tvec, xvec(:,3),'k','lineWidth',2);

        plot(cvec, xvec(:,1),'k','lineWidth',2);
        plot(cvec, xvec(:,2),'k--','lineWidth',2);
        plot(cvec, xvec(:,3),'k','lineWidth',2);
        SVariableCriticalZone = xvec(:,2);
        criticalSValue = max(SVariableCriticalZone);
        yline(criticalSValue, 'g--', 'LineWidth', 1.5); % Plot criticalSValue as a horizontal line
        cvec = [];
        xvec = [];
        for idxc = idx2+1:length(c)

            eqn = r*x*(1-x/K) - c(idxc)*x^2 / (x^2+h^2);
            xsol_all = solve(eqn==0);
            xsol_all = double(xsol_all);
            xsol = [];
            for jjj = 1:length(xsol_all)
                if isreal(xsol_all(jjj))
                    if abs(xsol_all(jjj)) ~= 0
                        xsol = [xsol,xsol_all(jjj)];
                    end
                end
            end
            if length(xsol) ~= 1
                warning('Check code')
            end
            cvec = [cvec,c(idxc)];
            xvec = [xvec,xsol];

        end
        plot(cvec, xvec,'k','lineWidth',2);

        %ylim([0,14])
        xlabel('Control variable')
        ylabel('System state')
        set(gca,'FontSize', 16)
        set(gca,'TickDir','out')
        set(gca,'ticklength',2*get(gca,'ticklength'))
        set(gca,'lineWidth',2)
        hold off;
        xlim([1, max(cvec)])

    end
    
    axes(ax1)
    yline(criticalSValue, 'g--', 'LineWidth', 1.5)
    


 end
    if ifreturn3sol
        c_3sol_output = c_3sol;
    else
        SVariableCriticalZone = [];
    end 
else
    disp('Control parameter is empty!')
end


