function y = rect(x,A,B) 

% x time axis
% A left limit of the rect
% B right limit of the rect
%
% y rect of width B-A and amplitude 1

% init vector of zeros, same length as the time axis
y = zeros(size(x));   

% find the rect window
set = find(x>=A & x<=B); 
% set the rect window to 1
y(set) = 1;

return