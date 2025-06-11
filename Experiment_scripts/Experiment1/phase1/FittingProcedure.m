% Calculate the frequence of "white" response for each constant
data=importdata('Frequence_kevin.txt',' ');
mat=data.data;
Rho=[0.6,0.625,0.65,0.675,0.7];
nr=length(Rho);
subNo=97;

% Fitting with a cumulative gaussian using the Plamedes toolbox
StimLevels=mat(:,2); % data values for the x axis 
NumPos=mat(:,3); % the number of trials in which the observer gave a correct response 
OutOfNum=mat(:,4); %the nb of trials for each stimulus level
PF= @PAL_CumulativeNormal; % a Cumulative gaussian is used to fit the data
paramsValues=[0.52 50 0 0]; % correspond to the paramateres of the fitting function alpha, beta, gamma, and lambda
paramsFree=[1 1 1 1]; % the first two parameteres are free and the two others are fixed

% The curve fitting procedure 
[paramsValues LL exitflag]=PAL_PFML_Fit(StimLevels,NumPos,OutOfNum,paramsValues,paramsFree,PF,'lapseLimits',[0 1],'guessLimits', [0 1]);

% Plot the graph showing the data and the smooth fitted function
StimLevelsFine=[min(StimLevels):(max(StimLevels)-min(StimLevels))/1000:max(StimLevels)];
Fit=PF(paramsValues,StimLevelsFine);
plot(StimLevels,mat(:,5),'k.','markersize',40);
set(gca, 'fontsize',12);
hold on;
plot(StimLevelsFine,Fit,'g-','linewidth',4);
axis([0.4 0.6 0 1]);
TITLE(['The data with the fitted function (mu=0.5 ; sigma=0.04)']);
XLABEL('% of White');
YLABEL('P(White response)');

% Calculate the stimilus corresponding to rho and 1-rho
phase2stimulus=strcat('StimulusForPhase2_',num2str(subNo),'.txt'); % name of data file to write the 2 resulting stimulus corresponding to rho and 1-rho 
phase2stimuluspointer = fopen(phase2stimulus,'wt'); % open ASCII file for writing
fprintf(phase2stimuluspointer,'SubNo rho WhiteState BlackState\n'); % Heading of the frequence file

% parameters of the psychometric function as in the book "Psychophysics" 
alpha=paramsValues(1);
beta=paramsValues(2);
gamma=paramsValues(3);
lamda=paramsValues(4);
% Change of the parameters to allow the use of normcdf function and its inverse
mu0=alpha;
sigma0=1/beta;
%WhiteState=StimLevelsFine(find(abs(Fit-0.6)<0.005, 1 ));
%BlackState=StimLevelsFine(find(abs(Fit-0.4)<0.005, 1 ));
for i=1:nr
WhiteState(i)=norminv((Rho(i)-gamma)/(1-gamma-lamda),mu0,sigma0);
BlackState(i)=norminv((1-Rho(i)-gamma)/(1-gamma-lamda),mu0,sigma0);

% Write results to file: 
fprintf(phase2stimuluspointer,'%i %1.3f %1.4f %1.4f\n', ...
        subNo, ...
        Rho(i), ...
        WhiteState(i), ...
        BlackState(i));
end

% Cleanup at end of experiment - Close window, show mouse cursor, close
% result file, switch Matlab/Octave back to priority 0 -- normal
% priority:
fclose('all');

