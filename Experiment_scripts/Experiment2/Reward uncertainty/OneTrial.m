function bal2=OneTrial(trial,state,rho,r0,w,ifi,sq,cir,file,No,hd,bal)
%balance=OneTrial(i,stateVect(i),rho,r0,w,ifi,square,circle,datafilepointer,subNo,hand,balance); 

% trial is the nb of the current trial 
% state is the state of the patch in the current trial (0 for black/dark or 1 for white/light)
% rho is the uncertainty level (65%)
% r0 is the value of positive reward 
% w is the name used to call the screen
% sq 0/or 1 depending on the association between the square shape and the light response
% cir 0/or 1 depending on the association between the cercle shape and the light response
% file is the '.txt' file where data will be written 
% No is the nb of the subject
% hd is the mapping between the stimilus and the keyboard answer (takes 1 or any other nb )
% bal is the current balance

% wait a bit between trials
WaitSecs(0.75);
[xc,yc]=WindowCenter(w);

% initialize KbCheck and variables to make sure they're
% properly initialized/allocted by Matlab - this to avoid time
% delays in the critical reaction time measurement part of the
% script:
[KeyIsDown, endrt, KeyCode]=KbCheck;

% when you're ready to start monitoring for the response (e.g., until  clicked)
quit=KbName('q');

%%%%%%%%%%%%%%%%%%%%%%%%
% Necessary parameters
%%%%%%%%%%%%%%%%%%%%%%%%

% Settings and parameters (necessary computations done before displaying the stimuli to increase timing precision)

% Duration of the different stim presentations in secs.
fixperiod=0.750; %750 ms for the fixation display
testperiod=0.250; %250 ms for the stimulus
rewardperiod=0.750; %750 ms for the reward display

CONST=10; % determines the size of the mouse cursor

% draw a cercle and a square randomly situated on a circle (C) of a centre C(xc,yc) and of a radius r  
r=250; % radius of (C)
d=175; % the length of a side of the square/circle 
theta1=0.5*pi*rand(1); % will be used to generate the center of the circle from (C) using the polar equation of a circle
theta2=0.5*pi*rand(1); % will be used to generate the center of the circle from (C) using the polar equation of a circle
pos1=round(rand(1));
pos2=round(rand(1));
if (pos1==0 && pos2==0)
    xo1=xc+r*cos(theta1);
    yo1=yc-r*sin(theta1);
    xo2=xc-r*cos(theta2);
    yo2=yc+r*sin(theta2);
elseif (pos1==0 && pos2==1)
    xo1=xc-r*cos(theta1);
    yo1=yc-r*sin(theta1);
    xo2=xc+r*cos(theta2);
    yo2=yc+r*sin(theta2); 
elseif (pos1==1 && pos2==0)
    xo1=xc-r*cos(theta1);
    yo1=yc+r*sin(theta1);
    xo2=xc+r*cos(theta2);
    yo2=yc-r*sin(theta2); 
elseif (pos1==1 && pos2==1)
    xo1=xc+r*cos(theta1);
    yo1=yc+r*sin(theta1);
    xo2=xc-r*cos(theta2);
    yo2=yc-r*sin(theta2); 
end

circle0 = [xo1-d/2,yo1-d/2,xo1+d/2,yo1+d/2]; % circle response
square0 = [xo2-d/2,yo2-d/2,xo2+d/2,yo2+d/2]; % square response 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displaying the stimulus, the response options and the reward feedback 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Displaying the stimilus as an image got from a matrix of pixels  
mydotmatrix=255*state*ones(400,400); %this will become our patch (400 by 400 pixels)
tex=Screen('MakeTexture', w, mydotmatrix);

%draw fixation display
Screen('DrawLine',w,[255 0 0],xc-15,yc,xc+15,yc,3);
Screen('DrawLine',w,[255 0 0],xc,yc-15,xc,yc+15,3);
vbl=Screen('Flip',w);

% Draw texture image to backbuffer. It will be automatically
% centered in the middle of the display if you don't specify a
% different destination:
Screen('DrawTexture', w, tex);

% Show the fixation cross at next possible display refresh cycle,
% and record stimulus onset time in 'startrt':
[vbl, startrt]=Screen('Flip', w, vbl+fixperiod); % show the fixation cross for 'fixperiod' secs ('fixperiod' would have passed from 'vbl' til 'vbl+fixperiod' time of first display of the fixation cross)

% display the response option (this first display is outside the while loop
% below for a more precise estimation of reaction time)
Screen('FillOval',w,[0 0 0],[xc-CONST,yc-CONST,xc+CONST,yc+CONST]); %draw the cursor yourself: indicate mouse position with a black, filled circle (or choose whatever!)
Screen('FrameOval',w,[255 255 255],circle0,3);
Screen('FrameRect', w, [255 255 255],square0,3);

% Show the stimulus at next possible display refresh cycle,
% and record stimulus onset time in 'startrt':
[vbl, startrt]=Screen('Flip', w, vbl+testperiod); % show the stimulus for 'testperiod' secs 

%centre the mouse
mousex=xc; 
mousey=yc;
SetMouse(mousex,mousey,w);
[mousex,mousey,buttons]=GetMouse(w); %obtain current mouse co-ordinates

startrt = GetSecs; % record time when the response options were first displayed (useful for computing reaction time)

while (KeyCode(quit)==0)
    
    % detect the response of the subject
    % The square reponse
    if ((buttons(1)) && (mousex>=xo2-d/2) && (mousex<=xo2+d/2) && (mousey>=yo2-d/2) && (mousey<=yo2+d/2))
        resp=sq;
        % compute response time
        rt = round(1000*(GetSecs - startrt));
        break;
    end
    % detect the circle reponse
    if (buttons(1) && ((mousex-xo1)^2+(mousey-yo1)^2<=(d/2)^2))   
        resp=cir;
        % compute response time
        rt = round(1000*(GetSecs - startrt));
        break;
    end

    [KeyIsDown, endrt, KeyCode]=KbCheck;
    
    % Display the response options in every video frame until a response is selected
    [mousex,mousey,buttons]=GetMouse(w); %obtain current mouse co-ordinates
    Screen('FillOval',w,[0 0 0],[mousex-CONST,mousey-CONST,mousex+CONST,mousey+CONST]); %draw the cursor yourself: indicate mouse position with a red, filled circle (or choose whatever!)
    Screen('FrameOval',w,[255 255 255],circle0,3);
    Screen('FrameRect', w, [255 255 255],square0,3);
    vbl=Screen('Flip',w,vbl+0.5*ifi); %do this for every video frame

end

% compute reward and accuracy ("correct" if light response selected with a
%white stimulus, or if a dark response selected with a dark stimulus)
if ( (resp==1 && state==1) || (resp==0 && state==0) )
    ac = 1;
    reward = r0*binornd(1,rho); % positive reward with freq rho (65%) if "light" response (for e.g., circle) chosen with "light" state
else
    ac = 0;
    reward = r0*binornd(1,1-rho); % no reward with freq rho (65%) if "dark" response (for e.g., square) chosen with "light" state 
end

% Update the balance
bal=bal+reward;

% Display reward feedback
% read stimulus image into matlab matrix 'imdata':
stimfilename=strcat('Images/',num2str(reward),'.jpg');
imdata=imread(stimfilename);

% make texture image out of image matrix 'imdata'
tex=Screen('MakeTexture', w, imdata);

% Draw texture image to backbuffer. It will be automatically
% centered in the middle of the display if you don't specify a
% different destination:
Screen('DrawTexture', w, tex);

% Clear screen to background color (our 'gray' as set at the
% beginning):
Screen('Flip', w);

% Update the display to show the reward feedback for "duration_reward" seconds
[vbl, startrt]=Screen('Flip', w, vbl+rewardperiod);

% Write trial result to file:
fprintf(file,'%i %i %i %i %i %i %i %i %i\n', ...
        No, ...
        trial, ...
        hd, ...
        state, ...
        resp, ...
        ac, ...
        reward, ...
        bal, ...
        rt);   

bal2=bal;

% Delete all textures
Screen('Close');

end
