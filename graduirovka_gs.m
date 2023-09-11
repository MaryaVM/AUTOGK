function fitresult = graduirovka_gs(b_gauss, e)
%CREATEFIT(B_GAUSS,E)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : b_gauss
%      Y Output: e
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 06-May-2023 15:41:50


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( b_gauss, e );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
fitresult = fit( xData, yData, ft );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'e vs. b_gauss', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'b_gauss', 'Interpreter', 'none' );
ylabel( 'e', 'Interpreter', 'none' );
grid on

