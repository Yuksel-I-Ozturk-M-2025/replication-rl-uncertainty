function OneTrialPhase1(trial,pw,w,ifi,tr,rect,file,No)
%OneTrialPhase1(i,pwhite(i),w,ifi,triangle,rectangle,datafilepointer,subNo);

% trial is the index of the current trial 
% pw is the proportion of white dots to be used for our display
% w is the name used to call the screen
% ifi is the monitor flip intervall of the screen window w (1 frame)
% tr 0/or 1 depending on the association between the triangle shape and the light response
% rect 0/or 1 depending on the association between the rectangle shape and the light response
% file is the '.txt' file where data will be written 
% No is the number of the subject

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
rewardperiod=0.750; %750 ms for the accuracy feedback

% to be used after for simplifying the reading of the results
if (pw>=0.5)
    state=1; % light state
else 
    state=0; % dark state
end

CONST=10; % determines the size of the mouse cursor

% draw a circle and a square randomly situated on a circle (C) of a centre C(xc,yc) and of a radius r  
r=250; % radius of (C)
d=125; % the length of a side of the square
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

% Edge points of the triangle shape to use (of center (xo2,yo2))
% First point
x1=xo2;
y1=yo2-d;
% Second point
x2=xo2+sqrt(3)/2*d;
y2=yo2+d/2;
% Third point
x3=xo2-sqrt(3)/2*d;
y3=yo2+d/2;

rectangle0=[xo1-d/2,yo1-d,xo1+d/2,yo1+d];
triangle0= [x1,y1;x2,y2;x3,y3];

% Slope and intercept of each the triangle sides (to use to detect mouse clicks inside the triangle)
% Slopes
a12=(y2-y1)/(x2-x1);
a13=(y3-y1)/(x3-x1);
a23=(y3-y2)/(x3-x2);
% Intercepts
b12=(x1*y2 - x2*y1)/(x1-x2);
b13=(x1*y3 - x3*y1)/(x1-x3);
b23=(x2*y3 - x3*y2)/(x2-x3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displaying the stimulus, the response options  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Displaying the stimilus as an image got from a matrix of pixels
mydotmatrix=zeros(400,400); %this will become our noisy image
imwidth=size(mydotmatrix,2);
imheight=size(mydotmatrix,1);
tmpvector=mydotmatrix(:); %turn into a vector
npixels=length(tmpvector); %number of pixels: 400 x 400
nwhite=round(npixels*pw); %number of white pixels
tmpvector(1:nwhite)=1;
tmpvector=tmpvector(randperm(npixels));
tmpvector=tmpvector*255; %turn zeros and ones into grey scale values
mydotmatrix=reshape(tmpvector,imheight,imwidth); %turn vector back into matrix with original dimensions
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
[vbl, startrt]=Screen('Flip', w, vbl+fixperiod); % show the fixation cross for 'fixperiod' secs 

% display the response option (this first display is outside the while loop
% below for a more precise estimation of reaction time)
Screen('FillOval',w,[0 0 0],[xc-CONST,yc-CONST,xc+CONST,yc+CONST]); %draw the cursor yourself: indicate mouse position with a red, filled circle (or choose whatever!)
Screen('FrameRect',w,[255 255 255],rectangle0, 3);
Screen('FramePoly', w, [255 255 255],triangle0, 3);

% Show the stimulus at next possible display refresh cycle,
% and record stimulus onset time in 'startrt':
[vbl, startrt]=Screen('Flip', w, vbl+testperiod); % show the stimulus for 'testperiod' secs 

%centre the mouse at the start of the trial
mousex=xc; 
mousey=yc;
SetMouse(mousex,mousey,w);
[mousex,mousey,buttons]=GetMouse(w); %obtain current mouse co-ordinates

startrt = GetSecs; % record time when the response options were first displayed (to compute reaction time)

while (KeyCode(quit)==0) % do this unless 'q' is pressed
    
    % detect the response of the subject
    % The rectangle reponse
    if (buttons(1) && (mousex>=xo1-d/2) && (mousex<=xo1+d/2) && (mousey>=yo1-d) && (mousey<=yo1+d))  
        resp=rect;
        % compute response time
        rt = round(1000*(GetSecs - startrt));
        break;
    end
    % The triangle reponse
    if ( (buttons(1)) && (mousey>a12*mousex+b12) && (mousey>a13*mousex+b13) && (mousey<a23*mousex+b23) ) % Triangle can be seen as a convex determined by 3 half planes
        resp=tr;
        % compute response time
        rt = round(1000*(GetSecs - startrt));
        break;
    end

    [KeyIsDown, endrt, KeyCode]=KbCheck;
    
    % Display the response options
    [mousex,mousey,buttons]=GetMouse(w); %obtain current mouse co-ordinates
    Screen('FillOval',w,[0 0 0],[mousex-CONST,mousey-CONST,mousex+CONST,mousey+CONST]); %draw the cursor yourself: indicate mouse position with a red, filled circle (or choose whatever!)
    Screen('FrameRect',w,[255 255 255],rectangle0, 3);
    Screen('FramePoly', w, [255 255 255],triangle0, 3);
    vbl=Screen('Flip',w,vbl+0.5*ifi); %do this for every video frame

end

% compute accuracy
% record "correct" if light response selected with a white stimulus, or if a dark response selected with a dark stimulus
ac=0;
if ( (resp==1 && state==1) || (resp==0 && state==0) )
     ac=1;
end   


% Show feedback
if (ac == 1) % Display the word 'correct' in green

    DrawFormattedText(w, 'Correct', 'center', 'center', [0, 150, 0]);
    
else % Display the word 'incorrect' in red
    
    DrawFormattedText(w, 'Incorrect', 'center', 'center', [255, 0, 0]);

end

% Clear screen to background color (our 'gray' as set at the
% beginning):
vbl=Screen('Flip', w);

% Update the display to show the reward feedback for "duration_reward" seconds
[vbl, startrt]=Screen('Flip', w, vbl+rewardperiod);

% Write trial result to file: 
fprintf(file,'%i %i %1.4f %i %i %i %i\n', ...
        No, ...
        trial, ...
        pw,...
        state, ...
        resp, ...
        ac, ...
        rt);    
    
% Delete all textures
Screen('Close');

end
