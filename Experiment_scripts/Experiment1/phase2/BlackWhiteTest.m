%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Test Phase    %%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function BlackWhiteTest(subNo)
%BlackWhiteTest(subNo);

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


hand=round(rand(1));% generate a random mapping
% Use input variable "hand" to determine response mapping for this session.
% 1 means white and 0 means black
if (hand==1)
    square=1; % "white" response via a click inside the square 
    cercle=0; % "black" response via left a click inside the square
else
    square=0; % "black" response via a click inside the square 
    cercle=1; % "white" response via left a click inside the circle
end

%%%%%%%%%%%%%%%%%%%%%%
% file handling
%%%%%%%%%%%%%%%%%%%%%%

% Define filenames of input files and result file:
datafilename = strcat('BlackWhiteTest_',num2str(subNo),'.txt'); % name of data file to write to

% check for existing result file to prevent accidentally overwriting
% files from a previous subject/session (except for subject numbers > 20):
if subNo<99 && fopen(datafilename, 'rt')~=-1
    fclose('all');
    error('Result data file already exists! Choose a different subject number.');
else
    datafilepointer = fopen(datafilename,'wt'); % open ASCII file for writing
    fprintf(datafilepointer,'SubNo TrialNo Mapping square circle WhiteProp State Resp Ac Reward Balance\n'); % Heading of the data file  
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
    %gray=188; %GrayIndex(screenNumber); 
    gray=135;
    
    % Open a double buffered fullscreen window on the stimulation screen
    % 'screenNumber' and choose/draw a gray background. 'w' is the handle
    % used to direct all drawing commands to that window - the "Name" of
    % the window. 'wRect' is a rectangle defining the size of the window.
    % See "help PsychRects" for help on such rectangles and useful helper
    % functions:
    [w, wRect]=Screen('OpenWindow',screenNumber, gray);
    ifi=Screen('GetFlipInterval',w); %need this for accurate timing

    xCent = wRect(3)/2;
    yCent = wRect(4)/2;
    
    % Set text size (Most Screen functions must be called after
    % opening an onscreen window, as they only take window handles 'w' as
    % input:
    Screen('TextSize', w, 14);
    
    % Do dummy calls to GetSecs, WaitSecs, KbCheck to make sure
    % they are loaded and ready when we need them - without delays
    % in the wrong moment:
    KbCheck;
    WaitSecs(0.1);
    GetSecs;
    
    % Presenting the first message with instructions at the begining of the
    % First phase
    message = 'Test phase ...\n\nLearn the correct mapping between the light and dark states and the two response options\n\n... press mouse button to begin ...';
    
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
    [rho,WhiteState,BlackState,pwhite,ntrials,ntrialsBS,ntrialsAS]=DispParams;
    
    % Wait a second before starting trial
    WaitSecs(1.000);
    
    balance=4; % This represents the winning of the subject so far  
    % Present a sequence of (ntrials) stimilus per configuration (a given %
    % of white). A switch will occur at the middle  
    for i=1:ntrialsBS
%         if (i==round(ntrialsBS/2)+1)
%              WaitSecs(1.000);
%              DrawFormattedText(w, 'Press mouse button to continue', 'center', 'center', WhiteIndex(w));
%              Screen('Flip', w);
%              GetClicks(w);
%          end
        balance2=OneTrial(i,pwhite(i),WhiteState,w,ifi,square,cercle,datafilepointer,subNo,hand,balance);
        balance=balance2;
    end
    hand=1-hand; % switch the mapping
    if (hand==1)
        square=1; % "white" response via a click inside the square 
        cercle=0; % "black" response via left a click inside the square
    else
        square=0; % "black" response via a click inside the square 
        cercle=1; % "white" response via left a click inside the circle
    end
    %hand=abs(round(rand(1))-hand); % the switch may occur or not     
    for i=ntrialsBS+1:ntrials
%         if (i==ntrialsBS+round(ntrialsAS/2)+1)
%              WaitSecs(1.000);
%              DrawFormattedText(w, 'Press mouse button to continue', 'center', 'center', WhiteIndex(w));
%              Screen('Flip', w);
%              GetClicks(w);
%          end
        balance2=OneTrial(i,pwhite(i),WhiteState,w,ifi,square,cercle,datafilepointer,subNo,hand,balance);
        balance=balance2;
    end
    
    % Show the the total winnings of the subject
    WaitSecs(1.000);
    Screen('TextSize', w, 20);
    message = ['Congratulation! Your total earnings is: £' num2str(balance),'\n\n- Thank you for your participation -\n\n\n\n\n\n... press mouse button to finish ...'];

    % Write the reward message, nicely centered in the
    % middle of the display, in white color. As usual, the special
    % character '\n' introduces a line-break:
    DrawFormattedText(w, message, 'center', 'center', WhiteIndex(w));

    % Update the display to show the instruction text:
    Screen('Flip', w);
    
    GetClicks(w);
%     fclose('all');
%     data=importdata(datafilename,' ');
%     mat=data.data;
%     whitemat=mat(mat(:,4)==WhiteState,:); % extract data containing only white states
%     l=length(whitemat(:,1)); % nb of rows of whitemat
%     k=10; % nb of trials used to calculate the average points
%     Rho=rho*ones(l);
%     Ac=whitemat(:,7);
%     AvgAc=zeros(l-k+1,1);
%     for j=0:l-k
%        AvgAc(j+1)=mean(Ac(1+j:k+j)); 
%     end
%     % Plot the graph to monitor learning  
%     clf;
%     plot((1:l),Ac,'bo'); % Plot the points corresponding to the accuracy
%     set(gca, 'fontsize',12);
%     hold on;
%     plot((1:l),Rho,'b--','linewidth',1); % Plot the straight line of equation y=rho
%     hold on;
%     plot((k:l),AvgAc,'r-','linewidth',1); % Plot the Average points to appreciate the learning curve
%     TITLE(['The learning graph (% white=',num2str(WhiteState),')']);
%     XLABEL('trial');
%     YLABEL('accuracy');
    
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