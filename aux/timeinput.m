
function output = timeinput(t,default_string,default_title)
% TIMEINPUT
% Input arguments:  
% t - time delay
% default_string - string which is returned if nothing is entered
%
% Examples:
% If a string is expected
%   x = timeinput(20,'no input')
% If a number is expected
%   x = str2num(timeinput(20,'1'))
%
% ADDOPTED from MATLAB Answers forum, answer by MATLAB support team:
% https://uk.mathworks.com/matlabcentral/answers/95301-how-can-i-implement-a-timer-dependent-input-request-in-matlab

if nargin < 2
   default_string = '';
end
if nargin < 3
   default_title = 'Please insert...';
end

% Creating a figure
h = figure('CloseRequestFcn','','Position',[500 500 450 50],'MenuBar','none',...
    'NumberTitle','off','Name',default_title);
% Creating an Edit field
hedit = uicontrol('style','edit','Units','pixels','Position',[50 15 250 20],'callback','uiresume','string',default_string);
uicontrol(hedit); % focus on the input

% Defining a Timer object
T = timer('Name','mm', ...
   'TimerFcn','uiresume', ...
   'StartDelay',t, ...
   'ExecutionMode','singleShot');

% Starting the timer
start(T)
uiwait(h)
% Defining the return value
output = get(hedit,'String');
% Deleting the figure
delete(h)
% Stopping and Deleting the timer
stop(T)
delete(T)
