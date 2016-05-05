function h = plot2Dframe( e1, e2, orgPt, c1, c2, s )

% Leonel Rozo, 2011
% 
% This function allows to graph a nice 2D frame of reference using
% different colors for each axis.
%   
%   e1:     The vector for the first axis    
%   e2:     The vector for the second axis
%   orgPt:  The origin of the frame
%   c1:     Color for axis e1
%   c2:     Color for axis e2
%   s:      Scale
%

% Handler
h = [];

% e1 axis
h = [ h quiver( orgPt(1), orgPt(2), e1(1), e1(2), s, ...
    'linewidth', 2, 'color', c1)];

% e2 axis
h = [ h quiver( orgPt(1), orgPt(2), e2(1), e2(2), s, ...
    'linewidth', 2, 'color', c2)];
