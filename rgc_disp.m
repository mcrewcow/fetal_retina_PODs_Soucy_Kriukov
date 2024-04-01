%close all
clear
clc
warning('off','all')


cd(uigetdir)
filelist = dir('*.csv');

n = size(filelist,1);

mapfile = dir('*.txt');
%mapdata = importdata(mapfile.name);
namelist = [];
for i = 1:length(filelist)
    templist = convertCharsToStrings(filelist(i).name);
    namelist = [namelist; templist];
end

cutoff = [10
37
33
33
31
33];

cell_VXYZ = [];

for w = [3]
    cell_VXYZ = [];
    tempstats = sprintf('JM5.M%d',w);
    tempstatindex = find(contains(namelist,tempstats));
    m = length(tempstatindex);
    
for s = 1:m
   
    data = readtable(sprintf('Statistics for JM5.M%d.OD_resfull_A01_G001_00%02d.oir - C=0.csv',w,s));

    mapdata = importdata(mapfile(1).name);

    tempname = sprintf('00%02d.vsi',s);
    tempindex = find(contains(mapdata,tempname));
    xmap = regexp(mapdata(tempindex+1),'[0-9]+','match');
    ymap = regexp(mapdata(tempindex+2),'[0-9]+','match');
    mapindex = [str2double(xmap{1}), str2double(ymap{1})];
    
    tempVXYZ = [data.Volume_micron_3_, data.X, data.Y, data.Z];
    
    for i = 1:length(tempVXYZ(:,1))
        if tempVXYZ(i,1) > 900
            tempVXYZ(i,5) = floor(tempVXYZ(i,1)/450);
        else
            tempVXYZ(i,5) = 1;
        end
        if tempVXYZ(i,4) > cutoff(w)
            tempVXYZ(i,6) = 1;
        else
            tempVXYZ(i,6) = 0;
        end
        tempVXYZ(i,2) = tempVXYZ(i,2)+512*(mapindex(1,1));
        tempVXYZ(i,3) = tempVXYZ(i,3)+512*(mapindex(1,2));        
    end
    
    cell_VXYZ = [cell_VXYZ;tempVXYZ];
    
end

%% Maximum number of points that can lie inside a circle
% radius = 1000;
% 
% XYvals = [1.243*cell_VXYZ(:,2),1.243*cell_VXYZ(:,3)];
% 
% groupindex = nchoosek(1:length(XYvals),3);
% Centers = zeros(length(groupindex),3);
% 
% for k = 1:size(groupindex, 1)
%     tempgroup = groupindex(k, :);
%     tempXYvals = XYvals(tempgroup, :);
%     [~,tempCenters] = fitcircle3(tempXYvals);
%     
%     x0 = XYvals(:,1) - tempCenters(1);
%     y0 = XYvals(:,2) - tempCenters(2);
%     insideobj = sum(x0.^2+y0.^2-radius^2<0);
%     Centers(k,:) = [tempCenters', insideobj];
%     
% end

% maxindex = find(Centers(:,3)==max(Centers(:,3)));
% x1 = XYvals(:,1) - Centers(maxindex(1),1);
% y1 = XYvals(:,2) - Centers(maxindex(1),2);
% 
% tempradius = 1:10000;
% tempobj = sum(x1.^2+y1.^2<tempradius.^2)';

%% Print to command window
totalcells = sum(cell_VXYZ(:,5))
GCLcells = sum(cell_VXYZ(:,6))/length(cell_VXYZ(:,6))
%injsitecells = max(Centers(:,3))/length(XYvals)
%[~,radius50] = (min(abs(tempobj - length(XYvals)/2)))

%% Figures
% figure()
% hold on
% scatter3(1.243*cell_VXYZ(:,2),1.243*cell_VXYZ(:,3),2*cell_VXYZ(:,4)-cutoff(w));
figure()
hold on
scatter(1.243*cell_VXYZ(:,2),1.243*cell_VXYZ(:,3),30,2*cell_VXYZ(:,4)-cutoff(w),'filled');
xlim([0 8000])
ylim([0 8000])
axis square
%viscircles(Centers(maxindex(1),1:2),radius,'Color','k');

end