function fitresult7 = createFit_Y88_2(t, isotope_spectrum2)
%CREATEFIT2(T,ISOTOPE_SPECTRUM1)
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

%  Auto-generated by MATLAB on 18-Apr-2023 13:38:59


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( t, isotope_spectrum2 );

% Set up fittype and options.
ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0];
opts.StartPoint = [663 840 42.0409107510239];

% Fit model to data.
fitresult7 = fit( xData, yData, ft, opts );

% Plot fit with data.
%figure( 'Name', 'untitled fit 1' );
%h = plot( fitresult2, xData, yData );
%legend( h, 'isotope_spectrum1 vs. t', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
%xlabel( 't', 'Interpreter', 'none' );
%ylabel( 'isotope_spectrum1', 'Interpreter', 'none' );
%grid on


