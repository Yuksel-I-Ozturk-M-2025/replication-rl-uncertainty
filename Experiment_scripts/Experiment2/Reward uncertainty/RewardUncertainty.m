%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Reward uncertainty condition    %%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function RewardUncertainty(subNo)
%RewardUncertainty(99);

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

hand=round(rand(1));% generate a random mapping

% Use input variable "hand" to determine response mapping for this session.
% 1 means white and 0 means black
if (hand==1)
    square=1; % "white" response via a click inside the square 
    circle=0; % "black" response via left a click inside the square
else
    square=0; % "black" response via a click inside the square 
    circle=1; % "white" response via left a click inside the circle
end

%%%%%%%%%%%%%%%%%%%%%%
% file handling
%%%%%%%%%%%%%%%%%%%%%%

% Define filenames of input files and result file:
datafilename = strcat('RewardUncertainty_',num2str(subNo),'.txt'); % name of data file to write to

% check for existing result file to prevent accidentally overwriting
% files from a previous subject/session (except for subject numbers > 20):
if subNo<99 && fopen(datafilename, 'rt')~=-1
    fclose('all');
    error('Result data file already exists! Choose a different subject number.');
else
    datafilepointer = fopen(datafilename,'wt'); % open ASCII file for writing
    fprintf(datafilepointer,'SubNo TrialNo Mapping State Resp Acc Reward Balance ReactionTime\n'); % Heading of the data file  
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
    
    % grey value to use:
    grey=174; % grey = GrayIndex(screenNumber); % Returns as default the mean grey value of screen
    
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
    
    % Set text size (Most Screen functions must be called after
    % opening an onscreen window, as they only take window handles 'w' as
    % input:
    Screen('TextSize', w, 24);
    
    % Presenting the first message with instructions at the begining of the
    % First phase
    message = 'Task starting ...\n\nLearn the correct mapping between the white and black states and the two response options\n\n... press mouse button to begin ...';
    
    % Write instruction message for subject, nicely centered in the
    % middle of the display, in white color. As usual, the special
    % character '\n' introduces a line-break:
    DrawFormattedText(w, message, 'center', 'center', WhiteIndex(w));

    % Update the display to show the instruction text:
    Screen('Flip', w);
        
    % Wait for mouse click:
    GetClicks(w);
    
    tic; % "tic; ... toc;" to estimate the duration of the experiment
    
    % Clear screen to background color (our 'gray' as set at the
    % beginning):
    Screen('Flip', w);
    
    % upload the display parameters
    [rho,r0,ntrials,ntrialsBS,ntrialsAS,stateVect,rate]=DispParams();
    
    % Wait a second before starting trial
    WaitSecs(0.5);
    
    balance=0; % This represents the winning of the subject so far (fixed payment of £4 to start with) 
    
    % Present a sequence of (ntrials) stimilus per configuration (a given % of white). A switch will occur at the middle 
    % procedure for each of the 100 trials prior to the switch 
    for i=1:ntrialsBS 
        balance=OneTrial(i,stateVect(i),rho,r0,w,ifi,square,circle,datafilepointer,subNo,hand,balance); %(the output is the subject's current balance)
    end
    
    hand=1-hand; % switch the mapping on the 101th trial
    if (hand==1)
        square=1; % "white" response via a click inside the square 
        circle=0; % "black" response via left a click inside the square
    else
        square=0; % "black" response via a click inside the square 
        circle=1; % "white" response via left a click inside the circle
    end

    % procedure for each of the 100 trials post switch
    for i=ntrialsBS+1:ntrials  
        balance=OneTrial(i,stateVect(i),rho,r0,w,ifi,square,circle,datafilepointer,subNo,hand,balance); %(the output is the subject's current balance)
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
    save(['Data_',num2str(subNo),'.mat'])
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