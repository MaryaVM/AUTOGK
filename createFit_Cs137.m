function fitresult8 = createFit_Cs137(t, isotope_spectrum1)
%CREATEFIT3(T,ISOTOPE_SPECTRUM1)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : t
%      Y Output: isotope_spectrum1
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 18-Apr-2023 20:17:33


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( t, isotope_spectrum1 );

% Set up fittype and options.
ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0];
opts.StartPoint = [28395 203 51.9127483009613];

% Fit model to data.
fitresult8 = fit( xData, yData, ft, opts );

% Plot fit with data.
%figure( 'Name', 'untitled fit 1' );
%h = plot( fitresult, xData, yData );
%legend( h, 'isotope_spectrum1 vs. t', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
%xlabel( 't', 'Interpreter', 'none' );
%ylabel( 'isotope_spectrum1', 'Interpreter', 'none' );
%grid on


