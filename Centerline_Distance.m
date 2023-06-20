
%% 1.
root_path = '/Users/amankc/Terminus_Ablation/';
cd(root_path)
glacier_name = 'Daugaard_Jensen_Gletsjer';
glacier_shp = append(glacier_name,'.shp');
centerline_path = fullfile(root_path,glacier_name,'Centerline',['Cenline_',glacier_shp]);
cline = shaperead(centerline_path);

% Run this if you are getting input as a line feature
x = cline.X';
y = cline.Y';
x = rmmissing(x);%removes NaN values
y = rmmissing(y);%removes NaN values

%determine ideal number of points for semi-standard spacing across all glaciers
full_length = sqrt((x(end)-x(1)).^2 + (y(end)-y(1)).^2);
dx = 2; %ideal spacing in meters
inds = round(full_length/dx);

%Apply the function to get evenly sapced points
[pt,dudt,fofthandle] = interparc(inds,x,y,'linear');
points = rmmissing(pt);

%convert the centerline coordinates to along-profile distance from the origin
term(1).center_dist = 0;
for k = 2:length(points)
term(k).center_dist = term(k-1).center_dist+ sqrt((points(k,1)-points(k-1,1)).^2 + (points(k,2)-points(k-1,2)).^2);
end

%instead of loading centerline intersections, load full terminus traces
%terminuspos_path = fullfile(root_path,glacier_name,'Merged_Termini',glacier_shp);
terminuspos_path = fullfile(root_path,glacier_name,'Terminus_Positions',glacier_shp);

S = shaperead(terminuspos_path);

%Assigning quality flag to the authors
names = {'Brough';'Bjork';'Black';'Black_Taryn';'Bunce';'Catania';'Fahrner';'Hill';'Korsgaard';'Moon_Twila';'Sole';'Termpicks';'Wood';'Zhang';
    'Carr';'Murray';'Bevan';'Cheng';'Cheng_D';'ESA';'PROMICE'};
qual_flag = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,3,3,3,3,3];

for j = 1:length(S)
    match = strcmp(S(j).Author,names);
    S(j).flag = qual_flag(match==1);
end
[~,index] = sortrows([S.DecDate].'); 
S = S(index); 
clear index

%% 2. Quality control for multiple traces on a same day
Save = {};
mm=1;
i=1;
while i<=length(S)
    n = 2;
    temp = S(i).DecDate;
    skip = 1;
    index = zeros(0);
    for j = i+1:length(S)
        if temp == S(j).DecDate
            index(1) = S(i).flag;
            index(n) = S(j).flag;
            n = n+1;
        end
    end
    if length(index)>1
        [~,id] = min(index);
        temp2 = S(i+id);
        skip = length(index);

    else
        temp2 = S(i);
    end
    Save {mm} = temp2;
    mm=mm+1;
    i = i + skip;
end
S_filt = cell2mat(Save);
S_filt = unique(S_filt,'rows'); %Gets rid if the two terminus traces are identical
%% 10. Getting the center velocity
S_fil = tePos;
for k = 1:length(tePos)
    decdate(k) = S_fil{1,k}.DecDate;
    termdate(k) = datenum(S_fil{1,k}.Date);
%     S(i).X = rmmissing(S(i).X); S(i).Y = rmmissing(S(i).Y); 
    [xi,yi,ii] = polyxpoly(points(:,1),points(:,2),S_fil{1,k}.X',S_fil{1,k}.Y'); %xi, yi are intersection coordinates
    if ~isempty(ii)
        centerdist(k) = term(ii).center_dist;
    else
        centerdist(k) = NaN;
    end
    clear ii;
end
termdate(isnan(centerdist)) = []; centerdist(isnan(centerdist)) = []; 

%calculate terminus velocity from distance change over time PROBABLY DON'T
%NEED
%then filter according to if dTermdt > 3*speed using Jukes code

%% 11. Filtering using ITS_LIVE Velocity
%rate = glacier_name.velocity;
rate_thresh = 2*15;
Shp_filtered = S_fil;
ind= []; counter = 1;
for k = 1:length(centerdist)-2
    rate = (centerdist(k+1)-centerdist(k))/(termdate(k+1)-termdate(k));
    if rate>0
        if rate>rate_thresh
            ind(counter) = k+1;
            counter = counter+1;
        end
    end
    
    if rate<0
        date_next = termdate(k+2); dist_next = centerdist(k+2);
        rate_next = dist_next/date_next;
        if rate_next >0 && rate_next > rate_thresh
            ind(counter) = k+1;
            counter = counter+1;
        end
    end
end
%%
Shp_filtered(:,ind) = [];  % removing the rows with dips and spikes