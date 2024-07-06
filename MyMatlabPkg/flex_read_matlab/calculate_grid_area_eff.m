function [ar, vo] = calculate_grid_area(header)

% calculates volumn and area for FLEXPART OUTGRID adapted from outgrid_init.f
% June 07 Sabine Eckhardt
%

r_earth=6.371e6;
pi180=pi/180;
ar=zeros(header.numygrid,header.numxgrid);
vo=zeros(header.numygrid,header.numxgrid,header.numzgrid);

hzone=zeros(length(header.latp),1);
ylatp=header.latp+.5*header.dyout; % [-89, -88, ..., ]
ylatm=header.latp-.5*header.dyout; % [-90, -89, ..., ]

ind=find(ylatm<0 & ylatp>0);
hzone(ind)=header.dyout*r_earth*pi180;
cosfactp=cos(ylatp*pi180)*r_earth;
cosfactm=cos(ylatm*pi180)*r_earth;
ind=find(cosfactp<cosfactm);

hzone(ind)=(r_earth.^2-cosfactp(ind).^2).^.5-(r_earth.^2-cosfactm(ind).^2).^.5;
ind=find(cosfactp>=cosfactm);
hzone(ind)=(r_earth.^2-cosfactm(ind).^2).^.5-(r_earth.^2-cosfactp(ind).^2).^.5;
gridarea=2.*pi*r_earth*hzone*header.dxout/360;

for ix=1:header.numxgrid
    ar(:,ix)=gridarea;
end

vo(:,:,1)=ar*header.outheight(1);      
for kz=2:header.numzgrid
   vo(:,:,kz)=ar*(header.outheight(kz)-header.outheight(kz-1));
end  %kz

ar=ar';
