%% 1.

NearTerm = readtable([output_path,'Near_Term_Mass_Change_',glacier_name,'.csv']);
NearTermMassChange = NearTerm(:,1);
NearTermMassChange(:,2) = NearTerm(:,5);
NearTermMassChange.Properties.VariableNames{1} = 'Dates';
NearTermMassChange.Properties.VariableNames{2} = 'mass';
first_date = datetime(NearTermMassChange.Dates(1),'ConvertFrom','datenum'); 
len = length(NearTermMassChange.Dates);
last_date = datetime(NearTermMassChange.Dates(len),'ConvertFrom','datenum');
duration = calmonths(between(first_date,last_date));
tdt = datetime(year(first_date),month(first_date):duration+month(first_date)+1,1);%generating array of 
% first day of each month
t = datenum(tdt);
mass_rates = zeros(duration,4);%Storing dates, mass changes and number of termini position 
% used to build the average
varnames = {'dates','datenumb','mass'};

%% 3. Calculation of Terminus Ablation
disch = readtable('/Users/amankc/Terminus_Ablation/gate_D.csv');
date_dis = datetime(table2array(disch(:,1)));
Val_dis = table2array(disch(:,2:end));
disch2 = zeros(size(Val_dis,1)-1,2);
% Manually input Glacier ID
glacierID = 140;
count = 0;
for i = 1:size(Val_dis,2)
    count = count+1;
    if Val_dis(1,i) == glacierID
        disch2(:,1) = datenum(date_dis(2:end,1));
        disch2(:,2) = Val_dis(2:end,count)/365;
    end
end
td = datenum(datetime(year(first_date),month(first_date):duration+month(first_date)+1,16));%generating array of 
% middle day of each month
td = td(2:end-1);
disch2 = unique (disch2,'rows'); %Eliminating duplicated rows
C = td';
years = datetime(mass_rates(:,1),'convertfrom','datenum');
discharge = interp1(disch2(:,1),disch2(:,2),C,'linear','extrap');

%%
discharge = discharge1+discharge2;
%% 
indx = []; coun = 1;
for i = 1:length(mass_rates)
    mass_rates(i,5) = discharge(i);
end
%     if mass_rates(i,2)>discharge(i)
%          indx(coun) = i; % Getting the index/rows of values that are not possible
%          coun = coun+1;
%     end
% end
% 
% mass_rates(indx,:) = []; %Eliminating the rows
    
%%
for i = 1:length(mass_rates)
    if mass_rates (i,2) == 0; % In case of stable terminus ablation
        mass_rates(i,4) = -mass_rates(i,5);
    elseif mass_rates(i,2) < 0; % In case of retreating terminus ablation
        mass_rates(i,4) = -(mass_rates(i,5)-mass_rates(i,2));
    else mass_rates(i,2) > 0; % In case of advancing terminus ablation
        mass_rates(i,4) = -(mass_rates(i,5) - mass_rates(i,2));
    end
end
% plot(datetime(mass_rates(:,1),'convertfrom','datenum'),mass_rates(:,4),'LineWidth',2,'Color','blue');
% hold on
% plot(datetime(mass_rates(:,1),'convertfrom','datenum'),mass_rates(:,5),'LineWidth',2,'Color','red');
% xlabel('Time','FontSize',15); ylabel('Ice Flux (Gt/day)','FontSize',15);
% grid on
% ax = gca; 
% ax.FontSize = 25; 
% legend('Terminus Ablation','Ice Discharge','Location','southwest','FontSize',20)
writematrix(mass_rates,[output_path,strcat('Term_Ablation_',GlacierName,'.csv')]);

%% To get the positive terminus ablation
for i = 1:length(mass_rates)
    if mass_rates (i,2) == 0; % In case of stable terminus ablation
        mass_rates(i,4) = mass_rates(i,5);
    elseif mass_rates(i,2) < 0; % In case of retreating terminus ablation
        mass_rates(i,4) = (mass_rates(i,5)-mass_rates(i,2));
    else mass_rates(i,2) > 0; % In case of advancing terminus ablation
        mass_rates(i,4) = (mass_rates(i,5) - mass_rates(i,2));
    end
end
% plot(datetime(mass_rates(:,1),'convertfrom','datenum'),mass_rates(:,4),'LineWidth',2,'Color','blue');
% hold on
% plot(datetime(mass_rates(:,1),'convertfrom','datenum'),mass_rates(:,5),'LineWidth',2,'Color','red');
% xlabel('Time','FontSize',15); ylabel('Ice Flux (Gt/day)','FontSize',15);
% grid on
% ax = gca; 
% ax.FontSize = 25; 
% legend('Terminus Ablation','Ice Discharge','Location','southwest','FontSize',20)
mass_rates(:,6) = mass_rates(:,4)*365;
writematrix(mass_rates,[output_path,strcat('Term_Ablation_',GlacierName,'.csv')]);

%% find the number of points which have negative terminus ablation
count = 0;
inds = [];
for i = 1:length(mass_rates)
    if mass_rates(i,4)<0
        count = count+1;
        inds(count) = i;
    end
end
%% Plotting

%bar(datetime(mass_rates(300:end-10,1),'convertfrom','datenum'),mass_rates(300:end-10,4)*365,'LineWidth',1,'Color','blue');
bar(datetime(mass_rates(1:end-1,1),'convertfrom','datenum'),mass_rates(1:end-1,4)*365,'blue');
hold on
% low_lim = years(1);
% high_lim = years (length(years)-1);
% xticks(years(1):365*7:high_lim);
plot(datetime(mass_rates(1:end-1,1),'convertfrom','datenum'),mass_rates(1:end-1,5)*365,'LineWidth',1,'Color','red');
xlabel('Time'); ylabel('Ice Flux(Gt/year)');
grid on
ax = gca; 
xticks = get(gca,'ytick');
set(gca,'yticklabel',xticks);
% ax.YAxis.TickLabelFormat = '%d'; 
ax.FontSize = 15; 
legend('Terminus Ablation','Ice Discharge','Location','northwest','FontSize',10)
% add time 
exportgraphics(ax,[outFolderImages,'TA_',GlacierName,'.png'],'Resolution',300);
%% Extra plot
GlacierName = 'Saqqarliup_Sermia';
plot(datetime(TermAblationSaqqarliupSermia(400:end-14,1),'convertfrom','datenum'),TermAblationSaqqarliupSermia(400:end-14,2)*365,'LineWidth',1,'Color','blue');
hold on
plot(datetime(TermAblationSaqqarliupSermia(400:end-14,1),'convertfrom','datenum'),TermAblationSaqqarliupSermia(400:end-14,3)*365,'LineWidth',1,'Color','red');
xlabel('Time'); ylabel('Ice Flux(Gt/year)');
grid on
ax = gca; 
xticks = get(gca,'ytick');
set(gca,'yticklabel',xticks);
% ax.YAxis.TickLabelFormat = '%d'; 
ax.FontSize = 15; 
legend('Terminus Ablation','Ice Discharge','Location','southwest','FontSize',10)
exportgraphics(ax,[outFolderImages,'TA_2',GlacierName,'.png'],'Resolution',300);
%% Finding Error in Gate Meta
error_file = readtable([root_path,'/Mankoffs_Data/gate_err.csv']);
date_dis_err = datetime(table2array(error_file(:,1)));
Val_dis_err = table2array(error_file(:,2:end));
disch_err = zeros(size(Val_dis_err,1)-1,2);
% Manually input Glacier ID
glacierID = 140;
count = 0;
for i = 1:size(Val_dis_err,2)
    count = count+1;
    if Val_dis_err(1,i) == glacierID
        disch_err(:,1) = datenum(date_dis_err(2:end,1));
        disch_err(:,2) = Val_dis_err(2:end,count); %In GT/year
    end
end
td = datenum(datetime(year(first_date),month(first_date):duration+month(first_date)+1,16));%generating array of 
% middle day of each month
td = td(2:end-1);
disch2 = unique (disch_err,'rows'); %Eliminating duplicated rows
C = td';
years = datetime(mass_rates(:,1),'convertfrom','datenum');

for j = 1:length(years)
    err = find(disch2(:,1)>=(mass_rates(j,1)-1) & disch2(:,1)<(mass_rates(j,1)+1));
    if isempty (err)
        %linear error
        ind1 = find(disch2(:,1)<(mass_rates(j,1)-1),1);
        ind2 = find(disch2(:,1)>(mass_rates(j,1)+1),1);
        if isempty(ind1)
            ind1 = ind2;
        end
        %just do the average of two dates (initial and latest)
        discharg_err(j) = sqrt((disch2(ind1,2).^2)/4 + (disch2(ind2,2).^2)/4);
    else
        %weighted error
        ab = 0;
        for i = 1:length(err)
            ab =  ab + (disch2(err(i),2));
        end
        discharg_err(j) = sqrt(ab);
    end

end

%% Error in bed data
bed_file = [root_path,'/BedTopo/BedMachineGreenland.nc'];
x_bed = double(ncread(bed_file,'x'));
y_bed = double(ncread(bed_file,'y'));
bed_err = ncread(bed_file,'errbed');
bed_err = bed_err'; 
[Rast,RA] = readgeoraster([root_path,'/BedTopo/bed_error.tif']);
% Get image info: 
R = geotiffinfo([root_path,'/BedTopo/bed_error.tif']); 
% Get x,y locations of pixels: 
[x_loc,y_loc] = pixcenters(R); 
% Convert x,y arrays to grid: 
[X,Y] = meshgrid(x_loc,y_loc);
for i = 1:length(Poly)
    polygon = Poly(i);
    x_req = Poly{1, i}.Vertices(:,1);
    y_req = Poly{1, i}.Vertices(:,2);
    x_req1 = x_req';
    y_req1 = y_req';
    mask = inpolygon(X,Y,x_req,y_req);
    %plot(X(mask),Y(mask),'r+')
    err_value = Rast(mask);
    bed_error(i) = mean(err_value);
end

%% error in DEMs
arc_DEM_error = 2.28;
korsgaard_error = 5.4;
for i = 1:length(bed_error)
    thick_err(i,1) = (sqrt(bed_error(i)^2+DEM_error(i,1)^2));
    thick_err(i,2) = (sqrt(bed_error(i)^2+DEM_error(i,2)^2));
end
%% Error in Area
leng = length(Shp_filtered);
Satellites = strings(leng,1);
for i = 1:leng
    Satellites(i,1) = S_filt(i).Satellite;
end
Uniq_Sat = unique(Satellites);
for i=1:leng
    x_coord = Shp_filtered{1,i}.X;
    y_coord = Shp_filtered{1,i}.Y;
    d = diff([x_coord(:) y_coord(:)]);
    total_length = sum(sqrt(sum(d.*d,2)));
    Satellites = Shp_filtered{1,i}.Satellite;
    if regexp(Satellites, regexptranslate('wildcard', 'L*')) == 1
        pix_size = 30;
    elseif regexp(Satellites, regexptranslate('wildcard', 'S*')) == 1
        pix_size = 10;
    elseif regexp(Satellites, regexptranslate('wildcard', 'A*')) == 1
        pix_size = 50;
    else regexp(Satellites, regexptranslate('wildcard', 'R*')) == 1;
        pix_size = 8;
    end
    num_pixels = total_length/pix_size;
    Area_err(i) = num_pixels * pix_size^2;%in sq meters
end
for i = 1:leng %VolumeTS(i,3) is in kms
    Vol_err(i,1) = (abs(VolumeTS(i,4)*10^9))*sqrt((Area_err(i)./(AreaTS(i,2)*1e6)).^2+(thick_err(i,1)./(VolumeTS(i,3).*1000)).^2); % In m^3
    Vol_err(i,2) = (abs(VolumeTS(i,4)*10^9))*sqrt((Area_err(i)./(AreaTS(i,2)*1e6)).^2+(thick_err(i,2)./(VolumeTS(i,3).*1000)).^2);
end

NMC_err(:,1) = (abs(Vol_err(:,1)) * 917)/10^12; %In Gt
NMC_err(:,2) = (abs(Vol_err(:,2)) * 917)/10^12; %In gt
ds(1,1) = 1;
% Interpolate between the years before making the comparison
varnames2 = {'dates','datenumb','error_upp','error_down'};
TPa = table(datetime(NearTermMassChange.Dates,'ConvertFrom','Datenum'),datenum(NearTermMassChange.Dates),NMC_err(:,1),NMC_err(:,2),'VariableNames',varnames2);
[C,ia] = unique(TPa.datenumb);
TP1 = TPa(ia,:);
for j = 2:length(t)-2
    terms = find(TP1.datenumb>=t(j) & TP1.datenumb<t(j+1));
    days = t(j+1)-t(j);
    ds(j,1) = days;
    len = length(terms);
    dif = 0;
    if isempty (terms)
        %linear error
        ind1 = find(TP1.datenumb<t(j),1,'last');
        ind2 = find(TP1.datenumb>t(j),1,'first');
        day_diff = TP1.datenumb(ind2) - TP1.datenumb(ind1);
        %just do the average of two dates (initial and latest)
        iu1 = TP1.error_upp(ind1)/day_diff; id1 = TP1.error_down(ind1)/day_diff;
        iu2 = TP1.error_upp(ind2)/day_diff; id2 = TP1.error_down(ind2)/day_diff;
        err_upp(j-1) = sqrt(((iu1^2)/4) + ((iu2^2)/4));
        err_down(j-1) = sqrt(((id1^2)/4) + ((id2^2)/4));

    else
        %weighted error
        a = 0;b = 0;
        id1 = find(TP1.datenumb<t(j),1,'last');
        id2 = find(TP1.datenumb>t(j+1),1,'first');
        dss = t(j)-TP1.datenumb(id1);
        total = 0;
        dif = (TP1.datenumb(terms(1))-(t(j)-1)+dss);
%       ds = dss+dif;
        for i = 1:len
            if i==1
                a1(i) = (TP1.datenumb(terms(1))-(t(j)-1));
                u1 = TP1.error_upp(terms(i)); u2 = TP1.error_upp(id1);
                d1 = TP1.error_down(terms(i)); d2 = TP1.error_down(id1);
                iur(i) = (sqrt(u1^2+u2^2))/dif;
                idn(i) = (sqrt(d1^2+d2^2))/dif;
%                 iur(i) = TP1.error_upp(terms(i))/dif;
%                 idn(i) = TP1.error_down(terms(i))/dif;
            else
                a1(i) = TP1.datenumb(terms(i))-TP1.datenumb(terms(i-1));
                u1 = TP1.error_upp(terms(i)); u2 = TP1.error_upp(terms(i-1));
                d1 = TP1.error_down(terms(i)); d2 = TP1.error_down(terms(i-1));
                iur(i) = (sqrt(u1^2+u2^2))/a1(i);
                idn(i) = (sqrt(d1^2+d2^2))/a1(i);
                if i == len
                    u3 = TP1.error_upp(terms(len)); d3 = TP1.error_down(terms(len));
                    u4 = TP1.error_upp(id2); d4 = TP1.error_down(id2);
                    a1(len+1) = (t(j+1)-1)-(TP1.datenumb(terms(len)));
                    last_diff = TP1.datenumb(id2)-(TP1.datenumb(terms(len)));
                    iur(len+1) = (sqrt(u3^2+u4^2))/last_diff;
                    idn(len+1) = (sqrt(d3^2+d4^2))/last_diff;
                end
%                 iur(i) = TP1.error_upp(terms(i))/a1(i);
%                 idn(i) = TP1.error_down(terms(i))/a1(i);
            end
            %sum of numerator and denominator should be same
            %total dif should be equal to number of days
            total = days;
            
        end
        for i = 1:len
            avg_err_upp(i) = ((a1(i)/total)*iur(i))^2;
            avg_err_low(i) = ((a1(i)/total)*idn(i))^2;
            a = a + avg_err_upp(i);
            b = b + avg_err_low(i);
        end
        if len>1
            avg_err_upp(len+1) = ((a1(len+1)/total)*iur(len+1))^2;
            avg_err_low(len+1) = ((a1(len+1)/total)*idn(len+1))^2;
            a = a + avg_err_upp(len+1);
            b = b + avg_err_low(len+1);
        end

%         dif2 = t(j+1)-TP1.datenumb(terms(len));
%         if dif2>0
%             a1 = TP1.error_upp(terms(len));
%             a2 = TP1.error_down(terms(len));
%             a = a+((dif2/days)*e1)^2;
%             b = b+((dif2/days)*Te2)^2;
%         end
        %calculate the difference between two dates
        err_upp(j-1) = sqrt(a);
        err_down(j-1) = sqrt(b);
        clear a b iur idn a1 avg_err_upp avg_err_low;
     end
end
dis_err = discharg_err'; 
upp_err = (err_upp)*365;
down_err = (err_down)*365;
for i =1:length(mass_rates)-1
    TA_err(i,1) = sqrt(upp_err(i)^2+dis_err(i)^2);
    TA_err(i,2) = sqrt(down_err(i)^2+dis_err(i)^2);
end
writematrix(TA_err,[output_path,strcat('TA_Error_',GlacierName,'.csv')]);
dates_number = datetime(mass_rates(1:end-1,1),'ConvertFrom','datenum');
dates1 = [dates_number;flipud(dates_number)];
TA_uppdown = [mass_rates(1:end-1,4).*365+TA_err(1:end-1,1) ; flipud(mass_rates(1:end-1,4)).*365-flipud(TA_err(1:end-1,2))];
fill(dates1,TA_uppdown,[.7 .7 .7])
hold on;
plot(years(1:end-1),mass_rates(1:end-1,4).*365,'red')
xlabel('Time'); ylabel('Ice Flux(Gt/year)');
grid on
ax = gca; 
xticks = get(gca,'ytick');
set(gca,'yticklabel',xticks);
% ax.YAxis.TickLabelFormat = '%d'; 
ax.FontSize = 15; 
legend('Error bound','Terminus Ablation','Location','northwest','FontSize',10)
exportgraphics(ax,[outFolderImages,'TA_Error_',GlacierName,'.png'],'Resolution',300);
%% 2.
TP = table(datetime(NearTermMassChange.Dates,'ConvertFrom','Datenum'),datenum(NearTermMassChange.Dates),NearTermMassChange.mass,'VariableNames',varnames);
for j = 2:length(t)-2
    terms = find(TP.datenumb>=t(j) & TP.datenumb<t(j+1));
    num_term = length(terms);
    if isempty (terms) %if no data is found between the month, lineraly interpolate 
        % between two earliest and latest dates from the month
        mean_date = round(mean([t(j) t(j+1)])); %getting the middle date of the month
        mass_rates (j-1,1) = mean_date;
        %mass_rates(j,1) = datetime(mean_date,'ConvertFrom','datenum') (To
        %convert into regular dates) Make sure predefined arrays will
        %accept datetime values
        mass_rates (j-1,2) = linear_fun(j,TP,t);
        mass_rates (j-1,3) = 0;
    else
        mass_rates (j-1,1) = round(mean([t(j) t(j+1)]));
        mass_rates(j-1,2) = weighted_fun(terms,TP,t);
        mass_rates(j-1,3) = num_term;
    end
end
%Linearly interp
function rates = linear_fun (j,TP,t)
dm = TP.mass(find(TP.datenumb<t(j),1,'last'))-TP.mass(find(TP.datenumb>t(j),1,'first'));
dt = TP.datenumb(find(TP.datenumb<t(j),1,'last'))-TP.datenumb(find(TP.datenumb>t(j),1,'first'));
rates = dm/dt;
end

%Weighted Interp
function rates = weighted_fun (terms,TP,~)
        l = length(terms);
        first_day = datenum(dateshift(TP.dates(terms(1)),'start','month'));
        last_day = datenum (dateshift(TP.dates(terms(1)),'end','month'));
        first_mass_change = TP.mass(terms(1)) - TP.mass(terms(1)-1);
        first_date_change = TP.datenumb(terms(1)) - TP.datenumb(terms(1)-1);
        if l==1
            mass_change = zeros(l+1);
            date_change = zeros(l+1);
            mass_rate = zeros(l+1);
        else
            mass_change = zeros(l);
            date_change = zeros(l);
            mass_rate = zeros(l);
        end
        mass_rate(1) = first_mass_change/first_date_change;
        if l>1
            for k = 1:l-1
                mass_change(k) = TP.mass(terms(k+1)) - TP.mass(terms(k)-1);
                date_change(k) = TP.datenumb(terms(k+1)) - TP.datenumb(terms(k)-1);
                mass_rate(k+1) = mass_change(k)/date_change(k);
            end
        end
        last_mass_change = TP.mass(terms(l)+1) - TP.mass(terms(l));
        last_date_change = TP.datenumb(terms(l)+1) - TP.datenumb(terms(l));
        mass_rate(l+1) = last_mass_change/last_date_change;
        last_date_diff = last_day-TP.datenumb(terms(l));
%         else
%             mass_rate(l+1) = 0;
%             last_date_diff = 0;
%        end
        first_date_diff = TP.datenumb(terms(1))-first_day;
%       Differencing from the first day of the month
        if first_date_diff == 0
            avg_mass_rate = mass_rate(1)*0;
        else
            avg_mass_rate = mass_rate(1)*first_date_diff;
        end

        for i = 1:length(mass_change)
            avg_mass_rate = avg_mass_rate + mass_change(i);
        end
        if last_date_diff == 0
%           last_date_change = TP.datenumb(terms(l)+1) - TP.datenumb(terms(l));
            sum_avg_mass_rate = avg_mass_rate + (mass_rate(l+1)* 0);
        else
            sum_avg_mass_rate = avg_mass_rate + (mass_rate(l+1)* last_date_diff);
        end
        date_sum = first_date_diff + last_date_diff;
        for j = 1:length(date_change)
            date_sum = date_sum + date_change(j);
        end
        rates = sum_avg_mass_rate/date_sum;
end


