function s = pcboundary(varargin)
% Draw boundary lines of connected components plotted by pcolor function
% <Syntax>
% pcboundary(C)
% pcboundary(X,Y,C)
% pcboundary(ax,___)
% pcboundary(___,Name,Value)
% s = pcboundary(___)
% 
% <Input/Output>
% The following arguments are the same as that in pcolor function, see 
% pcolor for more details about the way in which they should be specified.
% C - Color matrix
% X - x-coordinates
% Y - y-coordinates
% ax - Target axes
% s - Surface object
% 
% <Optional Input>
% Name-Value pair arguments
% box - outlines and boundary lines between numeric cells and NaN cells
%     specified as 'off' or 'on' (default).
% color - Boundary line color
%     specified as an RGB triplet, a hexadecimal color code, a color name,
%     or a short name. Default is black.
% linestyle - Boundary line style
%     '-' (default) | '--' | ':' | '-.' | 'none'
% linewidth - Boundary line width
%     specified as a positive value in points. Default is 0.5.
% 
% <Examples>
% C = randi(3,20,20);
% C(4:6,4:7)=NaN;
% pcboundary(C,'box','off','color','r','linewidth',2);
% 
% <Copyright>
% Author:   Yang Liu
% Contact:  liuyang-y2003@foxmail.com
% Update:   2022-06-08
% Version:  1.0.4
% 
% See also pcolor

box = 'on'; 
color = 'k';
ls = '-';
lw = 0.5;

n = find(cellfun(@(x) ischar(x),varargin),1);
if isempty(n)
    n = length(varargin)+1;
end

k=n;
while k<=length(varargin)
    switch lower(varargin{k})
        case 'box'
            box = varargin{k+1};
        case 'color'
            color = varargin{k+1};
        case 'linestyle'
            ls = varargin{k+1};
        case 'linewidth'
            lw = varargin{k+1};
    end
    k = k+2;
end

s = pcolor(varargin{1:n-1}); 
s.LineStyle = 'none';
if isvector(s.XData)
    [X,Y] = meshgrid(s.XData,s.YData);
else
    X = s.XData;
    Y = s.YData;
end

hold on

xh = [reshape(X(:,1:end-1),1,[]);reshape(X(:,2:end),1,[])]; 
yh = [reshape(Y(:,1:end-1),1,[]);reshape(Y(:,2:end),1,[])]; 

B = [s.CData(1:end-1,:);nan(1,size(X,2))]-[nan(1,size(X,2));s.CData(1:end-1,:)]; 

B((isnan([s.CData(1:end-1,:);nan(1,size(X,2))])+isnan([nan(1,size(X,2));s.CData(1:end-1,:)]))==1)=strcmp(box,'on');
B(isnan(B))=0;
B = logical(B);

plot(xh(:,reshape(B(:,1:end-1),1,[])),yh(:,reshape(B(:,1:end-1),1,[])),'color',color,'linewidth',lw,'linestyle',ls,...
    'Marker','o','markersize',lw,'markerfacecolor',color,'markeredgecolor','none')

xv = [reshape(X(1:end-1,:),1,[]);reshape(X(2:end,:),1,[])];
yv = [reshape(Y(1:end-1,:),1,[]);reshape(Y(2:end,:),1,[])];
B = [s.CData(:,1:end-1),nan(size(X,1),1)]-[nan(size(X,1),1),s.CData(:,1:end-1)];
B((isnan([s.CData(:,1:end-1),nan(size(X,1),1)]) + isnan([nan(size(X,1),1),s.CData(:,1:end-1)]))==1)=strcmp(box,'on');
B(isnan(B))=0;
B = logical(B);
plot(xv(:,reshape(B(1:end-1,:),1,[])),yv(:,reshape(B(1:end-1,:),1,[])),'color',color,'linewidth',lw,'linestyle',ls,...
    'Marker','o','markersize',lw,'markerfacecolor',color,'markeredgecolor','none');
