function bal2=OneTrial(trial,pw,wstate,w,ifi,sq,cer,file,No,hd,bal)
%function OneTrial(trial,pw,w,wresp,bresp,file,No,hd)

% trial is the nb of the current trial 
% pw is pwhite vector containing the parametres to be used for our display
% w is the name used to call the screen
% wresp is the mapping between the black stimilus and the associated
% keyboard answer
% bresp is the mapping between the black stimilus and the associated keyboard answer
% file is the '.txt' file where data will be written 
% No is the nb of the subject
% hd is the mapping between the stimilus and the keyboard answer (takes 1 or any other nb )

% Duration of image presentation in secs.
duration=0.250;
foreperiod=43; %506 ms for the fixation display
testperiod=21; %247 ms for the stimulus
rewardperiod=129; %506 ms for the reward display
% wait a bit between trials
WaitSecs(0.5);
[xc,yc]=windowcenter(w);

% initialize KbCheck and variables to make sure they're
% properly initialized/allocted by Matlab - this to avoid time
% delays in the critical reaction time measurement part of the
% script:
[KeyIsDown, endrt, KeyCode]=KbCheck;

% to be used after for simplifying the reading of the results
if (pw==wstate)
    state=1;
else 
    state=0;
end

% Displaying the stimilus as an image got from a matrix of pixels
mydotmatrix=zeros(400,400); %this will become our noisy image
imwidth=size(mydotmatrix,2);
imheight=size(mydotmatrix,1);
tmpvector=mydotmatrix(:); %turn into a vector
npixels=length(tmpvector); %how many pixels? 128 x 128
nwhite=round(npixels*pw); %how many white pixels?
tmpvector(1:nwhite)=1;
tmpvector=tmpvector(randperm(npixels));
tmpvector=tmpvector*255; %turn zeros and ones into grey scale values
mydotmatrix=reshape(tmpvector,imheight,imwidth); %turn vector back into matrix with original dimensions
%Screen('PutImage', w, mydotmatrix, [xCent-imwidth/2 yCent-imheight/2 (xCent-1)+imwidth/2 (yCent-1)+imheight/2]); 
tex=Screen('MakeTexture', w, mydotmatrix);

%draw fixation display
%Screen('FillRect',w,135);
Screen('FillRect',w,188);
Screen('DrawLine',w,[255 0 0],xc-7,yc,xc+7,yc,2);
Screen('DrawLine',w,[255 0 0],xc,yc-7,xc,yc+7,2);
vbl=Screen('Flip',w);

Screen('DrawTexture', w, tex);
% Show stimulus on screen at next possible display refresh cycle,
% and record stimulus onset time in 'startrt':
[vbl, startrt]=Screen('Flip', w, vbl+(foreperiod-0.5)*ifi);


%Screen('FillRect',w,135);
Screen('FillRect',w,188);
vbl=Screen('Flip', w, vbl+(testperiod-0.5)*ifi);

CONST=10;
mousex=xc; %centre the mouse at the start of the trial
mousey=yc;
SetMouse(mousex,mousey,w);
%[mousex,mousey,buttons]=GetMouse(w);


% draw a cercle and a square randomly situated on a cercle (C) of a centre C(xc,yc) and of a radius r  
r=250; % radius of (C)
theta1=0.5*pi*rand(1); % will be used to generate the center of the cercle from (C) using the polar equation of a cercle
theta2=0.5*pi*rand(1); % will be used to generate the center of the cercle from (C) using the polar equation of a cercle
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
d=100; % the length of a side of the square

% when you're ready to start monitoring for the response (e.g., until
% clicked)
quit=KbName('q');

while (KeyCode(quit)==0)
    
 [mousex,mousey,buttons]=GetMouse(w); %obtain current mouse co-ordinates
 Screen('FillOval',w,[0 0 0],[mousex-CONST,mousey-CONST,mousex+CONST,mousey+CONST]); %draw the cursor yourself: indicate mouse position with a red, filled circle (or choose whatever!)
 cercle0=[xo1-d/2,yo1-d/2,xo1+d/2,yo1+d/2];
 Screen('FrameOval',w,[255 255 255],cercle0);
 square0= [xo2-d/2,yo2-d/2,xo2+d/2,yo2+d/2];
 screen('FrameRect', w, [255 255 255],square0);
 messageBalance0=['current balance: £' num2str(bal)];
 Screen('DrawText', w, messageBalance0 , 0, 0, [0, 255, 0, 255]);
 vbl=Screen('Flip',w,vbl+0.5*ifi); %do this for every video frame
 
 % detect the response of the subject
 % detect the reponse for the square
 if ((buttons(1)) && (mousex>=xo2-d/2) && (mousex<=xo2+d/2) && (mousey>=yo2-d/2) && (mousey<=yo2+d/2))
     resp=sq;
     break;
 end
 % detect the reponse for the cercle
 if (buttons(1) && ((mousex-xo1)^2+(mousey-yo1)^2<=(d/2)^2))   
     resp=cer;
     break;
 end
 [KeyIsDown, endrt, KeyCode]=KbCheck;

end

% compute accuracy
ac=0;
% code correct if white-response with more-white stimulus,
% Black-response with more-blak stimulus
if ( (resp==1 && state==1) || (resp==0 && state==0) )
     ac=1;
end

% Allocating a reward to each pair (state,action)
reward=0;
if (state==1 && ac==1)
    reward=0.025;
elseif (state==1 && ac==0)
    reward=0;
elseif (state==0 && ac==1)
    reward=0.025;
elseif (state==0 && ac==0)
    reward=0;        
end

if (reward~=0)
bal=bal+reward;    
% Show the reward after subject's response
Screen('TextSize', w, 20);
message = ['You have got the following reward: £' num2str(reward)];

% Write the reward message, nicely centered in the
% middle of the display, in white color. As usual, the special
% character '\n' introduces a line-break:
DrawFormattedText(w, message, 'center', 'center', WhiteIndex(w));

% Update the display to show the instruction text:
[vbl, startrt]=Screen('Flip', w, vbl+(rewardperiod-0.5)*ifi);
end

% Write trial result to file:
fprintf(file,'%i %i %i %1.4f %i %i %i %1.3f %1.3f\n', ...
        No, ...
        trial, ...
        hd, ...
        pw, ...
        state, ...
        resp, ...
        ac, ...
        reward, ...
        bal);   

% Clear screen to background color after fixed 'duration'
% or after subjects response (on test phase)
bal2=bal;
end
