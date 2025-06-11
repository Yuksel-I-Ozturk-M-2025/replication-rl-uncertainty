function FitPsychometricFunction(subNo)
%FitPsychometricFunction(1);


% Importing data from the text file where results for subject 'subNo' are recorded 
data=importdata(['FrequencePhase1_Main_',num2str(subNo),'.txt'],' ');
mat=data.data;

% Fitting with a cumulative gaussian using the Plamedes toolbox
StimLevels=mat(:,2); % contrasts (values for the x axis) 
NumPos=mat(:,3); % number of trials in which the observer chose the light response 
OutOfNum=mat(:,4); % number of trials for each stimulus level
PF= @PAL_CumulativeNormal; % a cumulative gaussian is used to fit the data
paramsValues=[0.5 50 0 0]; % correspond to the paramateres of the fitting function alpha, beta, gamma, and lambda
paramsFree=[1 1 1 1]; % All four parameteres are free

% The curve fitting procedure 
[paramsValues LL exitflag]=PAL_PFML_Fit(StimLevels,NumPos,OutOfNum,paramsValues,paramsFree,PF,'lapseLimits',[0 1],'guessLimits', [0 1]);

% Plot the graph showing the collected data and the fitted function
StimLevelsFine=[min(StimLevels):(max(StimLevels)-min(StimLevels))/1000:max(StimLevels)];
Fit=PF(paramsValues,StimLevelsFine);
plot(StimLevels,mat(:,5),'k.','markersize',40);
set(gca, 'fontsize',12);
hold on;
plot(StimLevelsFine,Fit,'g-','linewidth',4);
ylim([0,1]);
xlim([min(StimLevels),max(StimLevels)]);
title('Data with the fitted psychometric function');
xlabel('% of White');
ylabel('P(light response)');
%saveas(gcf,['FitPreTrain_',num2str(subNo),'.png'])
%close;


end