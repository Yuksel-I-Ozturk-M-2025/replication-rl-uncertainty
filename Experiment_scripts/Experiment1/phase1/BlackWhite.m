function BlackWhite(subNo,hand)
%BlackWhite(subNo,hand);

% Skip the screen test
Screen('Preference', 'SkipSyncTests', 2 );

% Clear Matlab window:
clc;

% check for Opengl compatibility, abort otherwise:
AssertOpenGL;

% Reseed the random-number generator for each expt.
rand('state',sum(100*clock));

%ListenChar(2); %disable keyboard output

% Make sure keyboard mapping is the same on all supported operating systems
% Apple MacOS/X, MS-Windows and GNU/Linux:
KbName('UnifyKeyNames');

% Use input variable "hand" to determine response mapping for this session.
% 1 means white and 0 means black
if (hand==0)
    triangle=1; % "white" response via a click inside the triangle 
    rectangle=0; % "black" response via left a click inside the rectangle
else
    triangle=0; % "black" response via a click inside the triangle 
    rectangle=1; % "white" response via left a click inside the rectangle
end



%%%%%%%%%%%%%%%%%%%%%%
% file handling
%%%%%%%%%%%%%%%%%%%%%%

% Define filenames of input files and result file:
datafilename = strcat('BlackWhite_',num2str(subNo),'.txt'); % name of data file to write to
freqfilename = strcat('Frequence_',num2str(subNo),'.txt'); % name of data file to write the frequence of white response to
phase2stimulus=strcat('StimulusForPhase2_',num2str(subNo),'.txt'); % name of data file to write the 2 resulting stimulus corresponding to rho and 1-rho 

% check for existing result file to prevent accidentally overwriting
% files from a previous subject/session (except for subject numbers > 20):
if subNo<99 && fopen(datafilename, 'rt')~=-1
    fclose('all');
    error('Result data file already exists! Choose a different subject number.');
else
    datafilepointer = fopen(datafilename,'wt'); % open ASCII file for writing
    fprintf(datafilepointer,'SubNo TrialNo WhiteProp Resp Ac\n'); % Heading of the data file  
end


%%%%%%%%%%%%%%%%%%%%%%
% experiment
%%%%%%%%%%%%%%%%%%%%%%

% Embed core of code in try ... catch statement. If anything goes wrong
% inside the 'try' block (Matlab error), the 'catch' block is executed to
% clean up, save results, close the onscreen window etc.
try

    % Hide the mouse cursor:
    HideCursor;
    
    % Get screenNumber of stimulation display. We choose the display with
    % the maximum index, which is usually the right one, e.g., the external
    % display on a Laptop:
    screens=Screen('Screens');
    screenNumber=max(screens);
    
    % Returns as default the mean gray value of screen:
    gray=188; %GrayIndex(screenNumber); 
    %gray=135;
    
    % Open a double buffered fullscreen window on the stimulation screen
    % 'screenNumber' and choose/draw a gray background. 'w' is the handle
    % used to direct all drawing commands to that window - the "Name" of
    % the window. 'wRect' is a rectangle defining the size of the window.
    % See "help PsychRects" for help on such rectangles and useful helper
    % functions:
    [w, wRect]=Screen('OpenWindow',screenNumber, gray);
    ifi=Screen('GetFlipInterval',w); %need this for accurate timing
    
    % Set text size (Most Screen functions must be called after
    % opening an onscreen window, as they only take window handles 'w' as
    % input:
    Screen('TextSize', w, 18);
    
    % Do dummy calls to GetSecs, WaitSecs, KbCheck to make sure
    % they are loaded and ready when we need them - without delays
    % in the wrong moment:
    KbCheck;
    WaitSecs(0.1);
    GetSecs;
    
    % Presenting the First message with instructions at the begining of the
    % First phase
    if (hand==0)
        message = 'Estimation phase ...\n\nUsing the mouse, click on the triangle if you think the patch is light\nclick on the rectangle if you think the patch is dark\n\n... press mouse button to begin ...';
    else
        message = 'Estimation phase ...\n\nUsing the mouse, click on the rectangle if you think the patch is light\nclick on the triangle if you think the patch is dark\n\n... press mouse button to begin ...';
    end
    
    % Write instruction message for subject, nicely centered in the
    % middle of the display, in white color. As usual, the special
    % character '\n' introduces a line-break:
    DrawFormattedText(w, message, 'center', 'center', WhiteIndex(w));

    % Update the display to show the instruction text:
    Screen('Flip', w);
        
    % Wait for mouse click:
    GetClicks(w);
    
    % Clear screen to background color (our 'gray' as set at the
    % beginning):
    Screen('Flip', w);
    
    % get the array of structures containing the display parameters
    [rho,mu,sigma,constants,pwhite,nconstants,ntrialsperconstant,ntrials]=DispParams;
    
    % Wait a second before starting trial
    WaitSecs(1.000);
    
    % Present a sequence of (ntrials) stimilus per configuration (a given % of white) 
    for i=1:ntrials
         if (i==round(ntrials/2)+1)
             WaitSecs(1.000);
             DrawFormattedText(w, '...Break...\n\nPress mouse button to continue', 'center', 'center', WhiteIndex(w));
             Screen('Flip', w);
             WaitSecs(5.000);
             GetClicks(w);
         end
         OneTrial(i,pwhite(i),w,ifi,triangle,rectangle,datafilepointer,subNo);
        %OneTrial(i,pwhite(i),w,whiteresp,blackresp,datafilepointer,subNo,hand);
    end
    
    fclose('all');
    % Calculate the frequence of "white" response for each constant
    freqfilepointer = fopen(freqfilename,'wt'); % open ASCII file for writing
    fprintf(freqfilepointer,'SubNo Constant NbWhiteResp NbTrialPerConstant WhiteFreqResp\n'); % Heading of the frequence file
    nbWhite=zeros(nconstants,1);
    nbBlack=zeros(nconstants,1);
    whitefreq=zeros(nconstants,1);
    data=importdata(datafilename,' ');
    mat=data.data;
    nconstants_tmp=unique(mat(:,3));
    for j=1:length(nconstants_tmp)
%         tmp=mat(mat(:,3)==nconstants_tmp(j),:);
%         nbWite(j)=sum(tmp(:,4));
        for i=1:ntrials
            if ( mat(i,3)==nconstants_tmp(j) && mat(i,4)==1 )
                whitefreq(j)=whitefreq(j)+(1/ntrialsperconstant);
                nbWhite(j)=nbWhite(j)+1;
            end
        end
        nbBlack(j)=ntrialsperconstant-nbWhite(j);
    % Write results to file: 
    fprintf(freqfilepointer,'%i %1.2f %i %i %1.2f\n', ...
            subNo, ...
            constants(j), ...
            nbWhite(j),...
            ntrialsperconstant,...
            whitefreq(j));     
    end
    
    fclose('all');

    % Fitting with a cumulative gaussian using the Plamedes toolbox
    StimLevels=constants'; % data values for the x axis 
    NumPos=nbWhite; % the number of trials in which the observer gave a correct response 
    OutOfNum=ntrialsperconstant*ones(size(StimLevels)); %the nb of trials for each stimulus level
    PF= @PAL_CumulativeNormal; % a Cumulative gaussian is used to fit the data
    paramsValues=[0.5 50 0 0]; % correspond to the paramateres of the fitting function alpha, beta, gamma, and lambda
    paramsFree=[1 1 1 1]; % the first two parameteres are free and the two others are fixed
    
    % The curve fitting procedure 
    [paramsValues LL exitflag]=PAL_PFML_Fit(StimLevels,NumPos,OutOfNum,paramsValues,paramsFree,PF,'lapseLimits',[0 1],'guessLimits', [0 1]);
    
    % Plot the graph showing the data and the smooth fitted function
%     StimLevelsFine=[min(StimLevels):(max(StimLevels)-min(StimLevels))/1000:max(StimLevels)];
%     Fit=PF(paramsValues,StimLevelsFine);
%     plot(StimLevels,whitefreq,'k.','markersize',40);
%     set(gca, 'fontsize',12);
%     hold on;
%     plot(StimLevelsFine,Fit,'g-','linewidth',4);
%     TITLE(['The data with the fitted function (mu=',num2str(mu),';sigma=',num2str(sigma),')']);
%     XLABEL('% of White');
%     YLABEL('P(White response)');
    
    % Calculate the stimilus corresponding to rho and 1-rho
    phase2stimulus=strcat('StimulusForPhase2_',num2str(subNo),'.txt'); % name of data file to write the 2 resulting stimulus corresponding to rho and 1-rho 
    phase2stimuluspointer = fopen(phase2stimulus,'wt'); % open ASCII file for writing
    fprintf(phase2stimuluspointer,'SubNo WhiteState BlackState\n'); % Heading of the frequence file

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
    WhiteState=norminv((rho-gamma)/(1-gamma-lamda),mu0,sigma0);
    BlackState=norminv((1-rho-gamma)/(1-gamma-lamda),mu0,sigma0);

    % Write results to file: 
    fprintf(phase2stimuluspointer,'%i %1.4f %1.4f\n', ...
            subNo, ...
            WhiteState, ...
            BlackState);   
    
    % Cleanup at end of experiment - Close window, show mouse cursor, close
    % result file, switch Matlab/Octave back to priority 0 -- normal
    % priority:
    Screen('CloseAll');
    ShowCursor;
    fclose('all');
    Priority(0);
    
    % End of experiment:
    return;
    
    catch
    % catch error: This is executed in case something goes wrong in the
    % 'try' part due to programming error etc.:
    
    % Do same cleanup as at the end of a regular session...
    Screen('CloseAll');
    ShowCursor;
    fclose('all');
    Priority(0);
    %ListenChar(0); %re-enable keyboard output
    
    % Output the error message that describes the error:
    psychrethrow(psychlasterror);
end % try ... catch %