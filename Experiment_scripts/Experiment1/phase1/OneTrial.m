function OneTrial(trial,pw,w,ifi,tr,rect,file,No)
%function OneTrial(trial,pw,w,wresp,bresp,file,No,hd)

% trial is the nb of the current trial 
% pw is pwhite vector containing the parametres to be used for our display
% w is the name used to call the screen
% ifi is the monitor flip intervall of the screen window w (1 frame)
% wresp is the keyboard button associated with the white stimulus 
% bresp is the keyboard button associated the black stimilus
% file is the '.txt' file where data will be written 
% No is the nb of the subject

% Duration of image presentation in secs.
duration=0.250;
foreperiod=43; %506 ms for the fixation display
testperiod=21; %247 ms for the stimulus
% wait a bit between trials
WaitSecs(0.5);
[xc,yc]=windowcenter(w);

% initialize KbCheck and variables to make sure they're
% properly initialized/allocted by Matlab - this to avoid time
% delays in the critical reaction time measurement part of the
% script:
[KeyIsDown, endrt, KeyCode]=KbCheck;

% to be used after for simplifying the reading of the results
if (pw>=0.5)
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

% Draw texture image to backbuffer. It will be automatically
% centered in the middle of the display if you don't specify a
% different destination:
Screen('DrawTexture', w, tex);

% Show the fixation cross at next possible display refresh cycle,
% and record stimulus onset time in 'startrt':
[vbl, startrt]=Screen('Flip', w, vbl+(foreperiod-0.5)*ifi);

%Screen('FillRect',w,135);
Screen('FillRect',w,188);
Screen('Flip', w, vbl+(testperiod-0.5)*ifi);

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
 %Screen('gluDisk', w, [255 0 0], xc+200, yc-200, 100);
 rectangle0=[xo1-d/2,yo1-d,xo1+d/2,yo1+d];
 Screen('FrameRect',w,[255 255 255],rectangle0);
 triangle0= [xo2,yo2-d;xo2+sqrt(3)/2*d,yo2+d/2;xo2-sqrt(3)/2*d,yo2+d/2];
 %poly = [xc, yc; xc-100, yc-100; xc+100, yc-100];
 %Screen('FramePoly', w, [255 255 255], poly);
 screen('FramePoly', w, [255 255 255],triangle0);
 vbl=Screen('Flip',w,vbl+0.5*ifi); %do this for every video frame
 
 % detect the response of the subject
 % detect the reponse for the rectangle
 if (buttons(1) && (mousex>=xo1-d/2) && (mousex<=xo1+d/2) && (mousey>=yo1-d) && (mousey<=yo1+d))  
     resp=rect;
     break;
 end
 % detect the reponse for the triangle
 if ((buttons(1)) && ((mousex-xo2)^2+(mousey-yo2)^2<=d^2)) 
     resp=tr;
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

% Write trial result to file: 
fprintf(file,'%i %i %1.2f %i %i\n', ...
        No, ...
        trial, ...
        pw,...
        resp, ...
        ac);    
    
% Clear screen to background color after fixed 'duration'
% or after subjects response (on test phase)    
