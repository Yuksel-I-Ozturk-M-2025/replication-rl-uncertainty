%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   State uncertainty condition    %%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function StateUncertainty(subNo,hand)
%StateUncertainty(99,0);
%StateUncertainty(99,1);

% Skip the screen test
Screen('Preference', 'SkipSyncTests', 2 );

% Clear Matlab window:
clc;

% check for Opengl compatibility, abort otherwise:
AssertOpenGL;

% Reseed the random-number generator for each expt.
rand('state',sum(100*clock));

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

% For warm-up trials
datafilenamePre = strcat('StateUncertaintyPhase1_PreTrain_',num2str(subNo),'.txt'); % name of data file to write to
freqfilenamePre = strcat('FrequencePhase1_PreTrain_',num2str(subNo),'.txt'); % name of data file to write the frequence of white response to
% check for existing result file to prevent accidentally overwriting
% files from a previous subject/session (except for subject numbers > 20):
if subNo<99 && fopen(datafilenamePre, 'rt')~=-1
    fclose('all');
    error('Result data file already exists! Choose a different subject number.');
else
    datafilepointerPre = fopen(datafilenamePre,'wt'); % open ASCII file for writing
    fprintf(datafilepointerPre,'SubNo TrialNo WhiteProp State Resp Acc ReactionTime\n'); % Heading of the data file  
end

% For main trials
datafilenameMain = strcat('StateUncertaintyPhase1_Main_',num2str(subNo),'.txt'); % name of data file to write to
freqfilenameMain = strcat('FrequencePhase1_Main_',num2str(subNo),'.txt'); % name of data file to write the frequence of white response to
phase2stimulus=strcat('StimulusForPhase2_',num2str(subNo),'.txt'); % name of data file to write the 2 resulting stimulus corresponding to rho and 1-rho 

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
    
    % grey value to use:
    grey=174; % This value needs to be changed depending on the luminence function of the monitor
    
    % Open a double buffered fullscreen window on the stimulation screen
    % 'screenNumber' and choose/draw a gray background. 'w' is the handle
    % used to direct all drawing commands to that window - the "Name" of
    % the window. 'wRect' is a rectangle defining the size of the window.
    % See "help PsychRects" for help on such rectangles and useful helper
    % functions:
    [w, wRect]=Screen('OpenWindow',screenNumber, grey);
    ifi=Screen('GetFlipInterval',w); %need this for accurate timing
    
    % Do dummy calls to GetSecs, WaitSecs, KbCheck to make sure
    % they are loaded and ready when we need them - without delays
    % in the wrong moment:
    KbCheck;
    WaitSecs(0.1);
    GetSecs;
    
    % initialize KbCheck and variables to make sure they're
    % properly initialized/allocted by Matlab - this to avoid time
    % delays in the critical reaction time measurement part of the
    % script:
    [KeyIsDown, endrt, KeyCode]=KbCheck;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%      Phase 1: Estimation      %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Set text size (Most Screen functions must be called after
    % opening an onscreen window, as they only take window handles 'w' as
    % input:
    Screen('TextSize', w, 24);
    
    % Presenting the First message with instructions at the begining of the
    % First phase
    if (hand==0)
        message = 'Phase 1 ...\n\nUsing the mouse, click on the triangle if you think the patch is light\n\nclick on the rectangle if you think the patch is dark\n\n... press mouse button to begin ...';
    else
        message = 'Phase 1 ...\n\nUsing the mouse, click on the rectangle if you think the patch is light\n\nclick on the triangle if you think the patch is dark\n\n... press mouse button to begin ...';
    end
    
    % Write instruction message for subject, nicely centered in the
    % middle of the display, in white color. As usual, the special
    % character '\n' introduces a line-break:
    DrawFormattedText(w, message, 'center', 'center', WhiteIndex(w));

    % Update the display to show the instruction text:
    Screen('Flip', w);
        
    % Wait for mouse click:
    GetClicks(w);
    
    tic; % "tic; ... toc;" to estimate the duration of the experiment
    
    % Clear screen to background color (our 'gray' as set at the beginning):
    Screen('Flip', w);
    
    % Set text size (the size to use for reward feedback)
    Screen('TextSize', w, 50);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pre-training (to estimate sigma) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % upload the display parameters
    [rho,constants_warm,pwhite_warm,nconstants_warm,ntrialsperconstant_warm,ntrials_warm]=DispParamsPreTraining;
       
    % Wait half a second before starting trial
    WaitSecs(0.5);
    
    % Present a sequence of (ntrials) stimilus per configuration (a given % of white) 
    for i=1:ntrials_warm  
        OneTrialPhase1(i,pwhite_warm(i),w,ifi,triangle,rectangle,datafilepointerPre,subNo); % Procedure for each trial
    end
    
    fclose('all');
    
    % Set text size (the size to use for other text presentation than feedback)
    Screen('TextSize', w, 24);
    
    % Break (to allow also for estimating sigma)
    WaitSecs(1.000);
    DrawFormattedText(w, '...Break...\n\nPress mouse button to continue', 'center', 'center', WhiteIndex(w));
    Screen('Flip', w);
    
    % Set text size (the size to use for reward feedback)
    Screen('TextSize', w, 50);
    
    %%%%%%%%%%%%%%%%%%%% Fitting a pyschometric function to the warm-up data to get sigma %%%%%%%%%%%%%%%%
    
    % Calculate the frequence of "white" response for each constant (i.e., contrast)
    freqfilepointerPre = fopen(freqfilenamePre,'wt'); % open ASCII file for writing
    fprintf(freqfilepointerPre,'SubNo Constant NbWhiteResp NbTrialPerConstant WhiteFreqResp\n'); % Heading of the frequence file
    
    % vectors containing the number of light/dark responses for each constant
    nbWhite_warm=zeros(nconstants_warm,1);
    nbBlack_warm=zeros(nconstants_warm,1);
    
    % vector containing light response frequencies for each constant
    whitefreq_warm=zeros(nconstants_warm,1);
    
    %importing necessary data from the result file
    data_warm=importdata(datafilenamePre,' ');
    mat_warm=data_warm.data;
    nconstants_warm_tmp=unique(mat_warm(:,3));
    
    % Computing number of white/black dots for each contrast
    for j=1:length(nconstants_warm_tmp)
        for i=1:ntrials_warm % use only main trials (and discard the warm-up trials)
            if ( mat_warm(i,3)==nconstants_warm_tmp(j) && mat_warm(i,5)==1 )
                whitefreq_warm(j)=whitefreq_warm(j)+(1/ntrialsperconstant_warm);
                nbWhite_warm(j)=nbWhite_warm(j)+1;
            end
        end
        nbBlack_warm(j)=ntrialsperconstant_warm-nbWhite_warm(j);
        
    % Write results to file: 
    fprintf(freqfilepointerPre,'%i %1.4f %i %i %1.4f\n', ...
            subNo, ...
            constants_warm(j), ...
            nbWhite_warm(j),...
            ntrialsperconstant_warm,...
            whitefreq_warm(j));     
    end
    
    fclose('all');

    % Fitting a cumulative gaussian using the Plamedes toolbox
    StimLevels_warm=constants_warm'; % data values for the x axis 
    NumPos_warm=nbWhite_warm; % the number of trials in which the observer gave a correct response 
    OutOfNum_warm=ntrialsperconstant_warm*ones(size(StimLevels_warm)); %the nb of trials for each stimulus level
    PF= @PAL_CumulativeNormal; % a Cumulative gaussian is used to fit the data
    paramsValues=[0.5 50 0 0]; % correspond to the paramateres of the fitting function alpha, beta, gamma, and lambda
    paramsFree=[1 1 1 1]; % the first two parameteres are free and the two others are fixed
    
    % The curve fitting procedure 
    [paramsValues LL exitflag]=PAL_PFML_Fit(StimLevels_warm,NumPos_warm,OutOfNum_warm,paramsValues,paramsFree,PF,'lapseLimits',[0 1],'guessLimits', [0 1]);
      
    % Calculate the stimilus corresponding to rho and 1-rho
    % parameters of the psychometric function as in the book "Psychophysics" 
    alpha_warm=paramsValues(1);
    beta_warm=paramsValues(2);
    gamma_warm=paramsValues(3);
    lamda_warm=paramsValues(4);
    
    % Change of the parameters to allow the use of normcdf function and its inverse
    mu_warm=alpha_warm;
    sigma_warm=1/beta_warm;
    WhiteState_warm=norminv((rho-gamma_warm)/(1-gamma_warm-lamda_warm),mu_warm,sigma_warm);
    BlackState_warm=norminv((1-rho-gamma_warm)/(1-gamma_warm-lamda_warm),mu_warm,sigma_warm);

    % End break
    GetClicks(w);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main trials (to estimate the contrasts to use in the test phase) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Calculate the frequence of "white" response for each constant (i.e., contrast)
    datafilepointerMain = fopen(datafilenameMain,'wt'); % open ASCII file for writing
    fprintf(datafilepointerMain,'SubNo TrialNo WhiteProp State Resp Acc ReactionTime\n'); % Heading of the data file 
    
    % upload the display parameters
    [rho,constants,pwhite,nconstants,ntrialsperconstant,ntrials,sigma_adj]=DispParamsMain(sigma_warm);
    
    % Present a sequence of (ntrials) stimilus per configuration (a given % of white) 
    for i=1:ntrials  
        % Break half way through
        if (i==round(ntrials/2)+1)
            WaitSecs(1.000);
            % Set text size
            Screen('TextSize', w, 24);
            DrawFormattedText(w, '...Break...\n\nPress mouse button to continue', 'center', 'center', WhiteIndex(w));
            Screen('Flip', w);
            WaitSecs(5.000);
            GetClicks(w);
            % Set text size
            Screen('TextSize', w, 50);
        end
        OneTrialPhase1(i,pwhite(i),w,ifi,triangle,rectangle,datafilepointerMain,subNo); % Procedure on each trial
    end
    
    fclose('all');
    
    %%%%%%%%%%%%%%%%%%%% Fitting a pyschometric function to the main data to get sigma %%%%%%%%%%%%%%%%
    
    % Calculate the frequence of "white" response for each constant (i.e., contrast)
    freqfilepointerMain = fopen(freqfilenameMain,'wt'); % open ASCII file for writing
    fprintf(freqfilepointerMain,'SubNo Constant NbWhiteResp NbTrialPerConstant WhiteFreqResp\n'); % Heading of the frequence file
    
    % vectors containing the number of light/dark responses for each constant
    nbWhite=zeros(nconstants,1);
    nbBlack=zeros(nconstants,1);
    
    % vector containing light response frequencies for each constant
    whitefreq=zeros(nconstants,1);
    
    %importing necessary data from the result file
    data=importdata(datafilenameMain,' ');
    mat=data.data;
    nconstants_tmp=unique(mat(:,3));
    
    % Computing number of white/black dots for each contrast
    for j=1:length(nconstants_tmp)
        for i=1:ntrials % use only main trials (and discard the warm-up trials)
            if ( mat(i,3)==nconstants_tmp(j) && mat(i,5)==1 )
                whitefreq(j)=whitefreq(j)+(1/ntrialsperconstant);
                nbWhite(j)=nbWhite(j)+1;
            end
        end
        nbBlack(j)=ntrialsperconstant-nbWhite(j);
        
    % Write results to file: 
    fprintf(freqfilepointerMain,'%i %1.4f %i %i %1.4f\n', ...
            subNo, ...
            constants(j), ...
            nbWhite(j),...
            ntrialsperconstant,...
            whitefreq(j));     
    end
    
    fclose('all');

    % Fitting a cumulative gaussian using the Plamedes toolbox
    StimLevels=constants'; % data values for the x axis 
    NumPos=nbWhite; % the number of trials in which the observer gave a correct response 
    OutOfNum=ntrialsperconstant*ones(size(StimLevels)); %the nb of trials for each stimulus level
    PF= @PAL_CumulativeNormal; % a Cumulative gaussian is used to fit the data
    paramsValues=[0.5 50 0 0]; % correspond to the paramateres of the fitting function alpha, beta, gamma, and lambda
    paramsFree=[1 1 1 1]; % the first two parameteres are free and the two others are fixed
    
    % The curve fitting procedure 
    [paramsValues LL exitflag]=PAL_PFML_Fit(StimLevels,NumPos,OutOfNum,paramsValues,paramsFree,PF,'lapseLimits',[0 1],'guessLimits', [0 1]);
        
    % Calculate the stimilus corresponding to rho and 1-rho
    % parameters of the psychometric function as in the book "Psychophysics" 
    alpha=paramsValues(1);
    beta=paramsValues(2);
    gamma=paramsValues(3);
    lamda=paramsValues(4);
    
    % Change of the parameters to allow the use of normcdf function and its inverse
    mu0=alpha;
    sigma0=1/beta;
    WhiteState=norminv((rho-gamma)/(1-gamma-lamda),mu0,sigma0);
    BlackState=norminv((1-rho-gamma)/(1-gamma-lamda),mu0,sigma0);

    % Creating a file to record the contrast to use in the second phase
    phase2stimulus=strcat('StimulusForPhase2_',num2str(subNo),'.txt'); % name of data file to write the 2 resulting stimulus corresponding to rho and 1-rho 
    phase2stimuluspointer = fopen(phase2stimulus,'wt'); % open ASCII file for writing
    fprintf(phase2stimuluspointer,'SubNo WhiteState BlackState Sigma Sigma_adjusted\n'); % Heading of the frequence file
    % Write results to file:
    fprintf(phase2stimuluspointer,'%i %1.4f %1.4f %1.4f %1.4f\n', ...
            subNo, ...
            WhiteState, ...
            BlackState, ...
            sigma_warm, ...
            sigma_adj);   
        
    % Set text size
    Screen('TextSize', w, 24);
    
    % Presenting the First message with instructions at the begining of the
    % First phase

    message_final = 'First phase completed. Thank you!!!...\n\n... press mouse button to go to the second phase ...';
    
    % Write instruction message for subject, nicely centered in the
    % middle of the display, in white color. As usual, the special
    % character '\n' introduces a line-break:
    DrawFormattedText(w, message_final, 'center', 'center', WhiteIndex(w));
    
    % Update the display to show the instruction text:
    Screen('Flip', w);
    
    % End break
    GetClicks(w);
    
    % Save file
    save(['DataPhase1_',num2str(subNo),'.mat'])  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%      Phase 2: Learning     %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    hand2=round(rand(1));% generate a random mapping to use
    % Use input variable "hand" to determine response mapping for this session.
    % 1 means white and 0 means black
    if (hand2==1)
        square=1; % "white" response via a click inside the square 
        circle=0; % "black" response via left a click inside the square
    else
        square=0; % "black" response via a click inside the square 
        circle=1; % "white" response via left a click inside the circle
    end
    
    % Define filenames of input files and result file:
    datafilename = strcat('StateUncertaintyPhase2_',num2str(subNo),'.txt'); % name of data file to write to

    % check for existing result file to prevent accidentally overwriting
    % files from a previous subject/session (except for subject numbers > 20):
    if subNo<99 && fopen(datafilename, 'rt')~=-1
        fclose('all');
        error('Result data file already exists! Choose a different subject number.');
    else
        datafilepointer = fopen(datafilename,'wt'); % open ASCII file for writing
        fprintf(datafilepointer,'SubNo TrialNo Mapping WhiteProp State Resp Acc Reward Balance ReactionTime\n'); % Heading of the data file  
    end
    
    % Presenting the first message with instructions at the begining of the
    % First phase
    message = 'Phase 2...\n\nLearn the correct mapping between the light and dark states and the two response options\n\n... press mouse button to begin ...';
    
    % Write instruction message for subject, nicely centered in the
    % middle of the display, in white color. As usual, the special
    % character '\n' introduces a line-break:
    DrawFormattedText(w, message, 'center', 'center', WhiteIndex(w));

    % Update the display to show the instruction text:
    Screen('Flip', w);
    
    % Wait half a second before starting trial
    WaitSecs(5.000);
        
    % Wait for mouse click:
    GetClicks(w);
    
    tic; % "tic; ... toc;" to estimate the duration of the experiment
    
    balance=0; % This represents the winning of the subject so far (fixed payment of £4 to start with)
    
    % Clear screen to background color (our 'gray' as set at the
    % beginning):
    Screen('Flip', w);
    
    % upload the display parameters
    [rho,r0,WhiteState,BlackState,pwhite,ntrials,ntrialsBS,ntrialsAS,rate]=DispParamsPhase2(subNo);
    
    % Wait half a second before starting trial
    WaitSecs(0.5); 
    
    % Present a sequence of (ntrials) stimilus per configuration (a given % of white). A switch will occur at the middle 
    % procedure for each of the 100 trials prior to the switch 
    for i=1:ntrialsBS 
        balance=OneTrialPhase2(i,pwhite(i),WhiteState,w,ifi,square,circle,datafilepointer,subNo,hand2,r0,balance); %(the output is the subject's current balance)
    end
    
    hand2=1-hand2; % switch the mapping on the 101th trial
    if (hand2==1)
        square=1; % "white" response via a click inside the square 
        circle=0; % "black" response via left a click inside the square
    else
        square=0; % "black" response via a click inside the square 
        circle=1; % "white" response via left a click inside the circle
    end

    % procedure for each of the 100 trials post switch
    for i=ntrialsBS+1:ntrials  
        balance=OneTrialPhase2(i,pwhite(i),WhiteState,w,ifi,square,circle,datafilepointer,subNo,hand2,r0,balance); %(the output is the subject's current balance)
    end
    
    % Wait a second after the last trial
    WaitSecs(1.000);
    
    % Converting points to dollars (or pounds)
    payment_learn=rate*balance;
    
    % Show the the total winnings of the subject
    
    Screen('TextSize', w, 24); % Set text size 

    message = ['Congratulations! Your bonus earnings are: $' num2str(payment_learn,'%1.1f'),'\n\n- Thank you for your participation -\n\n\n\n\n\n... press mouse button to finish ...'];

    % Write the reward message, nicely centered in the
    % middle of the display, in white color. As usual, the special
    % character '\n' introduces a line-break:
    DrawFormattedText(w, message, 'center', 'center', WhiteIndex(w));

    % Update the display to show the instruction text:
    Screen('Flip', w);
    
    GetClicks(w);
    
    % Cleanup at end of experiment - Close window, show mouse cursor, close
    % result file, switch Matlab/Octave back to priority 0 -- normal
    % priority:
    save(['DataAll_',num2str(subNo),'.mat'])
    Screen('CloseAll');
    ShowCursor;
    fclose('all');
    Priority(0);
    
    toc; % Show elapsed time
    
    % End of experiment:
    return;
    
catch % catch error: This is executed in case something goes wrong in the 'try' part due to programming error etc.:
    
    % Do same cleanup as at the end of a regular session...
    Screen('CloseAll');
    ShowCursor;
    fclose('all');
    Priority(0);
    
    % Output the error message that describes the error:
    psychrethrow(psychlasterror);
end % try ... catch %

end