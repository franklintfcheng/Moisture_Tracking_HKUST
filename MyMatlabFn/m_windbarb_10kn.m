function h=m_windbarb_10kn(long,lat,u,v,varargin)
%  M_WINDBARB Project wind barbs onto map axes
%
%  M_WINDBARB(LONG,LAT,U,V) projects two dimensional wind barbs onto the 
%  current map axes. The vector components (U,V) are assumed to be in units 
%  of knots and are specified at the points (LON,LAT). It handles winds up 
%  to 130 knots. Winds exceeding 130 knots will appear as 130 knots.
%
%  If (u,v) are both positive, barbs appear to the SW of the vector head.
%  (in quiver, the arrowhead would be in the NE)
%
%  M_WINDBARB(...,S) uses the input S to scale the vectors after 
%  they have been automatically scaled to fit within the grid. If omitted, 
%  S = 0.9 is assumed.
%  
%  M_WINDBARB(...,'PropertyName',PropertyValue,...) allows you to specify
%  additional LINE properties.
% 
%  An additional parameter/value pair allows for u/v vectors in units
%  other than knots:
%          'units' : 'knots' | 'm/s' | 'kmh' | 'mph'
%

% Original code:
%  MFILE:   m_windbarb.m
%  MATLAB:  9.0.0 (R2016a)
%  VERSION: 1.0 (19 March 2017)
%  AUTHOR:  Erye
%  CONTACT: tfoterye@gmail.com
%
%  Oct/2017 - code "improved" for m_map style

% Modified code:
%  May/2020 - add a parameter for the reference value of a full barb
%  AUTHOR:  Fandy Cheng

global MAP_PROJECTION MAP_VAR_LIST

% Have to have initialized a map first

if isempty(MAP_PROJECTION)
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end

scale=0.9;
if nargin>4 && ~ischar(varargin{1})
    scale=varargin{1};
    varargin(1)=[];
end

if scale==0
    error(['map:' mfilename ':invalidScale'], ...
            'Invalid scale factor - must be greater than zero.')   
end
%      1 knot = 0.5144 m/s = 1.852 kmh = 1.151 mph

scf=1;  % assume knots
k=1;
while k<=length(varargin)
    switch lower(varargin{k}(1:3))
        case 'uni'
            switch lower(varargin{k+1})
                case 'knots'
                    scf=1;
                case {'m/s','meters/sec','meters/second'}
                    scf=0.5144;
                case 'kmh'
                    scf=1.852;
                case 'mph'
                    scf=1.151;
                otherwise
                    error(['map:' mfilename ':invalidUnits'], ...
                        'Units specified not recognized');
            end
            varargin([k k+1])=[];
        otherwise
            k=k+2;
    end
end


[X,Y]=m_ll2xy(long,lat,'clip','point');


[XN ,YN ]=m_ll2xy([long(:) long(:)]',[lat(:) lat(:)+.001]','clip','off');
[XE ,YE ]=m_ll2xy([long(:) long(:)+(.001)./cos(lat(:)*pi/180)]',[lat(:) lat(:)]','clip','off');
mU=u.*reshape(diff(XE),size(lat))*1000 + v.*reshape(diff(XN),size(lat))*1000;
mV=u.*reshape(diff(YE),size(lat))*1000 + v.*reshape(diff(YN),size(lat))*1000;

umag = sqrt(u.^2+v.^2)/scf ; %wind speed (should be in knots as input)
theta = atan2(mV,mU);

%create 18 logical matrices for 18 possible barbs. Non-zero when the barb
%is called for at that gridpoint.
% g{1} =  umag > 7.5  & umag <= 47.5;
% g{2} =  umag > 17.5 & umag <= 47.5;
% g{3} =  umag > 27.5;
% g{4} = (umag > 37.5 & umag <= 47.5) | (umag > 57.5 & umag <= 97.5);
% g{5} =  umag > 67.5;
% g{6} = (umag > 77.5 & umag <  97.5) | umag > 107.5;
% g{7} =  umag > 87.5 & umag < 97.5   | umag > 117.5;
% g{8} =  umag > 127.5;
% g{9} = (umag > 2.5  & umag <= 7.5 ) | (umag > 12.5 & umag <= 17.5);
% g{10} = umag > 22.5 & umag <= 27.5;
% g{11} =(umag > 32.5 & umag <= 37.5) | (umag > 52.5 & umag <= 57.5);
% g{12} =(umag > 42.5 & umag <= 47.5) | (umag > 62.5 & umag <= 67.5);
% g{13} =(umag > 72.5 & umag <= 77.5) | (umag > 102.5 & umag <= 107.5); 
% g{14} =(umag > 82.5 & umag <= 87.5) | (umag > 112.5 & umag <= 117.5);
% g{15} =(umag > 92.5 & umag <= 97.5) | (umag > 122.5 & umag <= 127.5);
% g{16} = umag > 47.5;
% g{17} = umag > 97.5;
% % g{18} = umag>0;    % All have background line

fullbarb = 10;
flag = fullbarb*5;
%/ full barb
g{1} =  umag >= fullbarb   & umag < flag;
g{2} =  umag >= fullbarb*2 & umag < flag;
g{3} =  umag >= fullbarb*3;
g{4} = (umag >= fullbarb*4 & umag < flag) | (umag >= flag+fullbarb & umag < flag*2);
g{5} =  umag >= fullbarb*7;
g{6} = (umag >= fullbarb*8 & umag < flag*2)  | (umag >= flag*2 + fullbarb);
g{7} = (umag >= fullbarb*9 & umag < flag*2)  | (umag >= flag*2 + fullbarb*2);
g{8} =  umag >= flag*2 + fullbarb*3;
%/ half barb
g{9} = (umag >= fullbarb*0.5 & umag < fullbarb) | (umag >= fullbarb*1.5 & umag < fullbarb*2); %/ h1=h2
g{10} = umag >= fullbarb*2.5 & umag < fullbarb*3; %/h3
g{11} =(umag >= fullbarb*3.5 & umag < fullbarb*4) | (umag >= flag+fullbarb*0.5 & umag < flag+fullbarb);%/h4=h6
g{12} =(umag >= fullbarb*4.5 & umag < flag) | (umag >= flag+fullbarb*1.5 & umag < flag+fullbarb*2); %/ h5=h7

g{13} =(umag >= fullbarb*7.5 & umag < flag+fullbarb*3) | (umag >= flag*2+fullbarb*0.5 & umag < flag*2+fullbarb); %/h8=h11
g{14} =(umag >= fullbarb*8.5 & umag < flag+fullbarb*4) | (umag >= flag*2+fullbarb*1.5 & umag < flag*2+fullbarb*2); %/h9=h12
g{15} =(umag >= fullbarb*9.5 & umag < flag*2)          | (umag >= flag*2+fullbarb*2.5 & umag < flag*2+fullbarb*3); %/h10=h13
g{16} = umag >= flag;
g{17} = umag >= flag*2;
g{18} = umag>0;    % All have background line

%position of each barb relative to grid point: [x0 y0; x1 y1]
c{1} = [-1    0; -1.125 .325];  % Full barb 1
c{2} = [-.875 0; -1     .325];  % Full barb 2
c{3} = [-.75  0; -.875  .325];  % Full barb 3
c{4} = [-.625 0; -.75   .325];  % Full barb 4
c{5} = [-.5   0; -.625  .325];  % Full barb 5
c{6} = [-.375 0; -.5    .325];  % Full barb 6
c{7} = [-.25  0; -.375  .325];  % Full barb 7
c{8} = [-.125 0; -.25   .325];  % Full barb 8
c{9} = [-.875 0; -.9375 .1625]; % Half barb 2
c{10} =[-.75  0; -.8125 .1625]; % Half barb 3
c{11} =[-.625 0; -.6875 .1625]; % Half barb 4
c{12} =[-.5   0; -.5625 .1625]; % Half barb 5
c{13} =[-.375 0; -.4375 .1625]; % Half barb 6
c{14} =[-.25  0; -.3125 .1625]; % Half barb 7
c{15} =[-.125 0; -.1875 .1625]; % Half barb 8
c{16} =[-1    0; -.875  .325];  % Slant other way for triangle 1
c{17} =[-.75  0; -.625  .325];  % Slant other way for taingle 
c{18} =[0     0; -1      0  ];  % Base

%set scale based on average latitude spacing
[m,n]=size(X);
scale2 = scale*(max(max(X))-min(min(X)))/n;

%draw the barbs
Ax=[];Ay=[];
for nn = 1:18
    
   ivals=g{nn}(:);
   if sum(ivals)>0  % number of barbs to draw
     
      %rotation operations
      cthet=cos(theta(ivals));
      sthet=sin(theta(ivals));
      
      x1=c{nn}(1,1)*cthet-c{nn}(1,2)*sthet;
      y1=c{nn}(1,1)*sthet+c{nn}(1,2)*cthet;
      x2=c{nn}(2,1)*cthet-c{nn}(2,2)*sthet;
      y2=c{nn}(2,1)*sthet+c{nn}(2,2)*cthet;
   
      x1 = x1*scale2+X(ivals);
      x2 = x2*scale2+X(ivals);
      y1 = y1*scale2+Y(ivals);
      y2 = y2*scale2+Y(ivals);
      x = [x1 x2 NaN(size(x1))]';   % Speed up by vectorizing this call with
      y = [y1 y2 NaN(size(y1))]';   %    Nans between line segments.
  
      Ax=[Ax;x(:)];
      Ay=[Ay;y(:)];
   end
end
h=line(Ax,Ay,varargin{:});

set(h,'tag','m_windbarb');

if nargout==0
 clear h
end
