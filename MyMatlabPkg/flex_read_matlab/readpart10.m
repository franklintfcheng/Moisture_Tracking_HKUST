function output = readpart10(directory,i,numpart,varargin)
% reads particle files from FLEXPART into structure output
% inputs: 
% directory [string]: path to part files 
% 	e.g. directory='/Users/userid/flexpart/FLEXPART8.1/output/'
% 
% i [integer or string]: order of particle file
% 	e.g. 1 for the firstdate in 'dates' file or  
%        '19800103120000' to read file='partposit_19800103120000'
%
% numpart [integer or string]: number of particles to be read. 
% 	If unknown input any string. This will determine numpart from the size
% 	of the file and return total nr of particles in the file instead of
%	the output given below. 
%
% varargin [integer]: optional number of species nspec. Default is 1. 
%
% output: 
%     outheader: [double] time record in seconds since simulation start
%        npoint: [numpart x 1 double] release point ID-number for each particle
%           xyz: [numpart x 3 double] position coordinates
%       itramem: [numpart x 1 double] relase times of each particle
%          vars: [numpart x 7 double] particle parameters
%         xmass: [numpart x nspec double] particle masses

%xyz (from flexpart fortran partout subroutine):
%1  xlon    =   longitude
%2  ylat    =   latitude
%3  ztra1(i)=   particle heigth in m
%vars:
%1  topo    =   topography heigth in m
%2  pvi     =   potential vorticity
%3  qvi     =   humidity
%4  rhoi    =   air density
%5  hmixi   =   'PBL' heigth in m
%6  tri     =   tropopause heigth in m
%7  tti     =   temperature
%xmass  =   particle mass (for nspec species)

% NB: If numpart < 1, only outheader is returned (output.outheader)

% Written by Ignacio Pisso, January 2008
% Modified
% 5/2011: added output.outheader IP
% 5/2013: modified to take more than one species into account
% aded varargin for more than 1 species. Defalut is 1. 4 arguments required

% Modified by Eivind G Waersted
% 19/12-2014: made npoint (release point ID) be returned for each trajectory
% 12/1-2015: vectorized the file-reading and made the function compute numpart directly from
% the size of the file (NB: it now reads large files ~100 times faster than before)
% 13/1-2015: handle the case that numpart < 1, and added closing of the file.
% 2018-08-09 adapted for distribution with flexpart version 10.3 IP

% t0 = cputime;

if i<1
    disp('i<0 ?')
    disp(i);
end

% determine nspecvi
nspec = 1;
nVarargs = length(varargin);
if nVarargs>0
   nspec = varargin{1}; 
end

% file name
if isnumeric(i) 
    dates=textread(strcat (directory,'/dates'));
    %i=2
    file=['partposit_' num2str(dates(i))];
elseif ischar(i)
    file=['partposit_' i];
end
file = strcat (directory,'/',file);

% t1 = cputime;

% if numpart not given, calculate numpart (total) and return it
if ischar(numpart)
    finfo = dir(file);
    nbytes = finfo.bytes;
    % assuming the file has the exact format written in partoutput.f90 of Flexpart v. 9.02,
    % the size of the file is exactly 4 * [3 + (numpart+1)*(14+nspec)] bytes.
    numpart = round((nbytes/4+3)/(14+nspec)-1);
    output = numpart;
    return
end

% if numpart is given, read the file and return the data
% file

[fid,message] = fopen(file,'r');
if fid < 0
    display(message)
end

% t2 = cputime;

% record time record
fread(fid,1,'int32');
output.outheader=fread(fid,1,'int32');
fread(fid,1,'int32');

% if numpart < 1, return here
if (numpart < 1)
  return
end

% FLEXPART writes the output as:
%           write(unitpartout) npoint(i),xlon,ylat,ztra1(i),
%      +    itramem(i),topo,pvi,qvi,rhoi,hmixi,tri,tti,
%      +    (xmass1(i,j),j=1,nspec)

% as 4byte numbers, with an empty 4byte value before and afterwards

% -> thus, the number of 4byte values read for each trajectory is 14+nspec
nvals = 14 + nspec;

% read data as int32 first, to get itramem and npoint
dataAsInt = fread(fid,[nvals,numpart],'int32');
fclose(fid);

% reread file, this time as float32, to get the rest of the data
fid = fopen(file,'r');
fread(fid,3,'int32'); % skip header
dataAsFloat = fread(fid,[nvals,numpart],'float32');
fclose(fid);

% t3 = cputime;

% create the output structure
output.npoint = transpose(dataAsInt(2,1:numpart));
output.xyz = transpose(dataAsFloat(3:5,1:numpart));
% output.itramem = transpose(dataAsInt(6,1:numpart)); %/ commented out (by FC)

%vars:    True index = relative index + 6. (by FC)
%1  topo    =   topography heigth in m
%2  pvi     =   potential vorticity
%3  qvi     =   humidity
%4  rhoi    =   air density
%5  hmixi   =   'PBL' heigth in m
%6  tri     =   tropopause heigth in m
%7  tti     =   temperature

vars_list = {'topo', 'pvi', 'q', 'rhoi', 'PBL', 'tri', 'T'};

ind_slct_vars = [9,11,13, 7];  %/ put 'topo' at last

output.vars = transpose(dataAsFloat(ind_slct_vars,1:numpart)); %/ only extract the needed columns. (by FC)
output.vars_name = vars_list(ind_slct_vars-6);

output.xmass = transpose(dataAsFloat(14:(13+nspec),1:numpart));



% t4 = cputime;

% Print how long time the function used on each part of the process
% fprintf('Performance:\n %15s: %.2f s\n %15s: %.2f s\n %15s: %.2f s\n %15s: %.2f s\n','Inputs',t1-t0,'Allocations',t2-t1,'Reading',t3-t2,'Structure',t4-t3)

fclose all;    %/ an attempt to resolve the bug "std::exception ... Message Catalog MATLAB:FileIO was not loaded from the file". 
               %/ See https://www.mathworks.com/matlabcentral/answers/124335-caught-std-exception-exception-message-is-message-catalog-matlab-builtins-was-not-loaded-from-th
return

