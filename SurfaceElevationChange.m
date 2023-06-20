function [DEMstruct,ChangeStruct,DEM_error]= SurfaceElevationChange(PathToLatestDEM,PathToOldestDEM,PathToBedrockTopo,FrontalAB,Poly)
% Function determines Surface elevation change per year based on two
% DEMS, and calculates ice thickness using surface elvation and bedrock 
% topography (automatically cropped to DEM extent)for each available 
% date of terminus positions. 
%
% arguments: (input)
% PathToLatestDEM - Path to directory of latest DEM. DEM has to be in Polar 
%           Stereographic projection (EPSG:3413).
% PathToOldestDEM - Path to directory of oldest available DEM.DEM has to be 
%           in Polar Stereographic projection (EPSG:3413).
% PathToBedrockTopo - Path the Bedmachine v4 data set.
% FrontalAB - Matrix containing area of ROIs (not used here) and dates of 
%           digitised termini
%
% Input raster data is automatically downsampled to the resolution of the
% bedrock topography data (150 m)
%
% arguments: (output)
% DEM.X - X coordinates of corresponding to caluclated ice thickness values
% DEM.Y - Y coordinates of corresponding to caluclated ice thickness values     
% DEM.H - Ice thickness values calculated assuming linear change in surface
%         elevation. Annual change of surface elevation is determined by
%         differencing the two DEMs and dividing values by number of years
%         between both. Surface elevation is heightend by adding the annual 
%         elevation change multiplied by the elapsed time between terminus 
%         observations for years earlier than the creation date of the
%         oldest DEM. Surface elevation is lowered by subtracting the annual 
%         elevation change multiplied by the elapsed time between terminus 
%         observations for years later than the creation date of the
%         oldest DEM. 
%
%
% NOTE: These calculations assume linear changes in surface elevation and
% therefore do not account for annual dynamic variability
%         

% PathToLatestDEM = 'Data/DEMs/Narsap_Sermia/ArcticDEM_2017_Narsap.tif';
% PathToOldestDEM = 'Data/DEMs/Narsap_Sermia/1985_DEM_NS_Clipped.tif';
% PathToBedrockTopo = 'Data/BedTopo/BedmachineV4_BedTopography.tif';

    [LatestDEM,Rlatest] =geotiffread(PathToLatestDEM);
    [OldDem,Rold] = geotiffread(PathToOldestDEM);
    [Bed,Rbed] = geotiffread(PathToBedrockTopo);



    OldDEMinfo = geotiffinfo(PathToOldestDEM);
    Bedinfo = geotiffinfo(PathToBedrockTopo);
    arc_DEM_error = 2.28;
    korsgaard_error = 5.4;
%%
    
    [ResampNewDEM,RnewResamp] = mapresize(LatestDEM,Rlatest,Rlatest.CellExtentInWorldX/Bedinfo.PixelScale(1),'nearest');
    [ResampOldDEM,RoldResamp] = mapresize(OldDem,Rold,Rold.CellExtentInWorldX/Bedinfo.PixelScale(1));
%     ResampOldDEM(ResampOldDEM>20000)=NaN;
%     ResampOldDEM(ResampOldDEM<-20000)=NaN;
%     xlimits = [RoldResamp.XWorldLimits(1),RoldResamp.XWorldLimits(2)];
%     ylimits = [RoldResamp.YWorldLimits(1),RoldResamp.YWorldLimits(2)];
%   
%     [CropNewDEM,RnewCrop] = mapcrop(ResampNewDEM,RnewResamp,xlimits,ylimits);
%     [CropBed,RbedCrop] = mapcrop(Bed,Rbed,xlimits,ylimits);
%     
%     if sum(size(CropNewDEM)~=size(ResampOldDEM))>0
%         SizeDifference = size(ResampOldDEM)-size(CropNewDEM);
%         CropNewDEM(:,end+SizeDifference(2))=zeros;
%         CropNewDEM(end+SizeDifference(1),:)=zeros;
%     end
%     if sum(size(CropBed)~=size(ResampOldDEM))>0
%         SizeDifference = size(CropBed)-size(ResampOldDEM);
%         CropBed = CropBed(:,1:end-SizeDifference(2));
%         CropBed = CropBed(1:end-SizeDifference(1),:);
%     end
    ResampNewDEM(ResampNewDEM>20000)=NaN;
    ResampNewDEM(ResampNewDEM<-20000)=NaN;
    sizeNewDEM = size(ResampNewDEM);
    sizeOldDEM = size(ResampOldDEM);

    if (sizeNewDEM<sizeOldDEM)==1
    % Cropping Older DEM & Bedrock Topography to extent of newer DEM
    xlimits = [RnewResamp.XWorldLimits(1),RnewResamp.XWorldLimits(2)];
    ylimits = [RnewResamp.YWorldLimits(1),RnewResamp.YWorldLimits(2)];
    
    [CropOldDEM,ROldCrop] = mapcrop(ResampOldDEM,RoldResamp,xlimits,ylimits);
    [CropBed,RbedCrop] = mapcrop(Bed,Rbed,xlimits,ylimits);
else
    xlimits = [RoldResamp.XWorldLimits(1),RoldResamp.XWorldLimits(2)];
    ylimits = [RoldResamp.YWorldLimits(1),RoldResamp.YWorldLimits(2)];
    
    [CropNewDEM,RCropNewDEM] = mapcrop(ResampNewDEM,RnewResamp,xlimits,ylimits);
    [CropBed,RbedCrop] = mapcrop(Bed,Rbed,xlimits,ylimits);
end
    if sum(size(CropOldDEM)~=size(ResampNewDEM))>0
        SizeDifference = size(CropOldDEM) - size(ResampNewDEM);
        CropOldDEM(:,end+SizeDifference(2))=zeros;
        CropOldDEM(end+SizeDifference(1),:)=zeros;
    end
    if sum(size(CropBed)~=size(ResampOldDEM))>0
        SizeDifference = size(CropBed)-size(ResampNewDEM);
        CropBed = CropBed(:,1:end-SizeDifference(2));
        CropBed = CropBed(1:end-SizeDifference(1),:);
    end
    %%% Check if all Images have the same x/y size
if exist('CropOldDEM','var')
    if isequal(CropBed,CropOldDEM,ResampNewDEM)==0 % Check if all files have the same size and if not: 
        sizeDifferenceDEMs = size(ResampNewDEM) - size(CropOldDEM); % Get size difference for DEMs
        if (sum(sizeDifferenceDEMs~=0))>0% if they are the same size
         % if there is a size difference between the DEMs, make them the same size and Bed Topo
                SizeMinDEM = min(size(CropOldDEM), size(ResampNewDEM));
                CropOldDEM = CropOldDEM(1:SizeMinDEM(1),1:SizeMinDEM(2));
                ResampNewDEM = ResampNewDEM(1:SizeMinDEM(1),1:SizeMinDEM(2));
        end
    end
        sizeDifferenceDEM_Bed = size(CropOldDEM) - size(CropBed); % Check size difference between Bed and DEM
        if (sizeDifferenceDEM_Bed~=0) ==1 % if there is a size difference
            SizeMin = min(size(CropOldDEM), size(CropBed));
            CropBed = CropBed(1:SizeMin(1),1:SizeMin(2));
        end
        RbedCrop.RasterSize = size(CropBed);
        RnewResamp.RasterSize = size(ResampNewDEM);
        ROldCrop.RasterSize = size(CropOldDEM);

end
if exist('CropNewDEM','var')
    if isequal(CropBed,ResampOldDEM,CropNewDEM)==0 % Check if all files have the same size and if not: 
        sizeDifferenceDEMs = size(ResampOldDEM) - size(CropNewDEM); % Get size difference for DEMs
        if (sizeDifferenceDEMs~=0)>0% if they are the same size
         % if there is a size difference between the DEMs, make them the same size and Bed Topo
                SizeMinDEM = min(size(CropNewDEM), size(ResampOldDEM));
                CropNewDEM = CropNewDEM(1:SizeMinDEM(1),1:SizeMinDEM(2));
                ResampOldDEM = ResampOldDEM(1:SizeMinDEM(1),1:SizeMinDEM(2));
        end
    end
    sizeDifferenceDEM_Bed = size(CropNewDEM) - size(CropBed); % Check size difference between Bed and DEM
    if (sizeDifferenceDEM_Bed~=0) ==1 % if there is a size difference
        SizeMin = min(size(CropNewDEM), size(CropBed));
        CropBed = CropBed(1:SizeMin(1),1:SizeMin(2));
    end
    RbedCrop.RasterSize = size(CropBed);
    RCropNewDEM.RasterSize = size(CropNewDEM);
    RoldResamp.RasterSize = size(ResampOldDEM);
end

    %% Find acquisition year of earliest and latest DEM and calculate yearly change
  
   
    k = strfind(PathToOldestDEM,'19');
    if ~isempty(k)
        EarliestDEM = str2num(PathToOldestDEM(k:k+3));
    else
        k= strfind(PathToOldestDEM,'20');
        EarliestDEM = str2num(PathToOldestDEM(k:k+3));
    end
    l = strfind(PathToLatestDEM,'20');
    LatestDEMYear = str2num(PathToLatestDEM(l:l+3));


    DEMYearDifference = LatestDEMYear-EarliestDEM;
    [x,y] = pixcenters(RbedCrop,size(CropBed));
    [xDEMfinal,yDEMfinal] = meshgrid(x,y);
    DEM{1}.X = xDEMfinal;
    DEM{1}.Y = yDEMfinal;
%    l = length(FrontalAB);
    ROImask{1} = inpolygon(DEM{1}.X,DEM{1}.Y,Poly{1}.Vertices(:,1),Poly{1}.Vertices(:,2));
    ROImask{243} = inpolygon(DEM{1}.X,DEM{1}.Y,Poly{l}.Vertices(:,1),Poly{l}.Vertices(:,2));
    DiffDEM = nanmean(CropOldDEM(ROImask{1}),'all') - nanmean(ResampNewDEM(ROImask{243}),"all");
    DiffDEM(DiffDEM<-1e+38)=NaN;
    YearlyDiff = DiffDEM/DEMYearDifference;
%     yearlyChange = nanmean(YearlyDiff(:));
    yearlyChange = YearlyDiff;
    dailyChange = yearlyChange./365;

    %% Get Dates of observations of terminus positions
    Tdates = datetime(single(FrontalAB(:,1)),'ConvertFrom','datenum');

    TdateLatestDEM = datetime(PathToLatestDEM(l:l+7),'Format','yyyyMMdd');

    DifferenceToLatestDEM.Hours=TdateLatestDEM-Tdates

    %%
    % Subtract/add surface elevation change form 1985 DEM based ontimedelta
    % between Terminus position observations and get Ice Thickness (H) by
    % subtracting surface elvation from bedrock elevation
    for i=1:length(FrontalAB)
         if exist('CropOldDEM','var')
             [x,y] = pixcenters(RbedCrop,size(CropBed));
             [xDEMfinal,yDEMfinal] = meshgrid(x,y);
             
                if year(Tdates(i))<LatestDEMYear
                    DEM{i}.X = xDEMfinal;
                    DEM{i}.Y = yDEMfinal;
                    
                    timeDelta = LatestDEMYear-year(Tdates(i)); 
                    fprintf('Date before Latest DEM \n')
                    fprintf('Adding surface height \n')
                    Change{i} = timeDelta*abs(yearlyChange);
                    NewSurfaceHeight{i} = ResampNewDEM + Change{i}; %
                    DEM{i}.H = NewSurfaceHeight{i} - CropBed; 
                    ROImask{i} = inpolygon(DEM{i}.X,DEM{i}.Y,Poly{i}.Vertices(:,1),Poly{i}.Vertices(:,2));
                    DEM_err(i,1) = nanmean((CropOldDEM(ROImask{i})-CropBed(ROImask{i})),'all') - nanmean((NewSurfaceHeight{i}(ROImask{i})-CropBed(ROImask{i})),"all");
                    DEM_err(i,2) = nanmean((NewSurfaceHeight{i}(ROImask{i})-CropBed(ROImask{i})),"all") - nanmean((ResampNewDEM(ROImask{i})-CropBed(ROImask{i})),"all");
                end
                if year(Tdates(i))>=LatestDEMYear
                    DEM{i}.X = xDEMfinal;
                    DEM{i}.Y = yDEMfinal;
                    fprintf('Date is the same as Latest DEM \n')
                    fprintf('No change to surface height \n')
                    Change{i} = [];
                    NewSurfaceHeight{i} = ResampNewDEM; %
                    DEM{i}.H = NewSurfaceHeight{i} - CropBed; 
                    DEM_err(i,1) = arc_DEM_error;
                    DEM_err(i,2) = arc_DEM_error;
                end
%                 if year(Tdates(i))>LatestDEMYear
%                     DEM{i}.X = xDEMfinal;
%                     DEM{i}.Y = yDEMfinal;
%                     timeDelta = year(Tdates(i))-LatestDEMYear;
%                     fprintf('Date after latest DEM \n')
%                     fprintf('subtracting surface height \n')
%                     Change{i} = (timeDelta*abs(yearlyChange));
%                     NewSurfaceHeight{i} = ResampNewDEM -(timeDelta*abs(yearlyChange));
%                     DEM{i}.H = NewSurfaceHeight{i} - CropBed;
%                 end
        else
            [x,y] = pixcenters(ResampNewDEM,size(ResampNewDEM));    
            [xDEMfinal,yDEMfinal] = meshgrid(x,y);
                if year(Tdates(i))<LatestDEMYear
                    DEM{i}.X = xDEMfinal;
                    DEM{i}.Y = yDEMfinal;
                    timeDelta = LatestDEMYear - year(Tdates(i)); 
                    fprintf('Date before Latest DEM \n')
                    fprintf('Adding surface height \n')
                    Change{i} = (timeDelta*abs(yearlyChange));
                    NewSurfaceHeight{i} = CropNewDEM +(timeDelta*abs(yearlyChange)) ;
                    DEM{i}.H = NewSurfaceHeight{i} - CropBed;
                    ROImask{i} = inpolygon(DEM{i}.X,DEM{i}.Y,Poly{round(i/2)}.Vertices(:,1),Poly{round(i/2)}.Vertices(:,2));
                    DEM_err(i,1) = nanmean(ReSampOldDEM) - nanmean(DEM{i}.H,"all");
                    DEM_err(i,2) = nanmean(DEM{i}.H,"all") - nanmean(CropNewDEM,"all");
                end
                if year(Tdates(i))>=LatestDEMYear
                    DEM{i}.X = xDEMfinal;
                    DEM{i}.Y = yDEMfinal;
                    fprintf('Date is the same as Latest DEM \n')
                    fprintf('No change to surface height \n')
                    Change{i} = [];
                    NewSurfaceHeight{i} = CropNewDEM;
                    DEM{i}.H = NewSurfaceHeight{i} - CropBed; 
                    DEM_err(i,1) = arc_DEM_error;
                    DEM_err(i,2) = arc_DEM_error;
                end
%                 if year(Tdates(i))>LatestDEMYear
%                     DEM{i}.X = xDEMfinal;
%                     DEM{i}.Y = yDEMfinal;
%                     timeDelta = year(Tdates(i))-LatestDEMYear; 
%                     fprintf('Date after latest DEM \n')
%                     fprintf('subtracting surface height \n')
%                     Change{i} = (timeDelta*abs(yearlyChange));
%                     NewSurfaceHeight{i} = CropNewDEM -(timeDelta*abs(yearlyChange));
%                     DEM{i}.H = NewSurfaceHeight{i} - CropBed;
%                 end
         end
     DEMstruct = DEM;
     ChangeStruct = Change;
     DEM_error = DEM_err;
end
%     for i=1:length(FrontalAB)
%         [x,y] = pixcenters(RoldResamp,size(ResampOldDEM));
%         [xDEM,yDEM] = meshgrid(x,y);
%         
%         if year(Tdates(i))<LatestDEMYear
%             DEM{i}.X = xDEM;
%             DEM{i}.Y = yDEM;
%             timeDelta = years(DifferenceToLatestDEM.Hours(i));
%             disp('Date before Latest DEM, adding surface height')
%             Change{i} = timeDelta*yearlyChange;
%             NewSurfaceHeight{i} = CropOldDEM - (timeDelta*yearlyChange) ;
%             DEM{i}.H = NewSurfaceHeight{i} - CropBed; 
%             ROImask{i} = inpolygon(DEM{i}.X,DEM{i}.Y,Poly{round(i/2)}.Vertices(:,1),Poly{round(i/2)}.Vertices(:,2));
%             DEM_err(i,1) = nanmean(ResampOldDEM(ROImask{i}) - DEM{i}.H(ROImask{i}),"all");
%             DEM_err(i,2) = nanmean(DEM{i}.H(ROImask{i}),"all") - nanmean(CropNewDEM(ROImask{i}),"all");
%             DEM_err(i,3) = nanmean( DEM{i}.H,"all");
%             DEM_err(i,4) = nanmean(ResampOldDEM(ROImask{i}));
%             DEM_err(i,5) = nanmean(CropNewDEM,"all");
%             DEM_err(i,6) = nanmean(CropBed,'all');
%             DEM_err(i,7) = nanmean(NewSurfaceHeight{i},'all');
%         end
% %         if year(Tdates(i))>LatestDEMYear
% %             DEM{i}.X = xDEM;
% %             DEM{i}.Y = yDEM;
% %             timeDelta = years(DifferenceToLatestDEM.Hours(i));
% %             disp('Date after latest DEM, subtracting surface height')
% %             Change{i} = (timeDelta*yearlyChange)*-1;
% %             NewSurfaceHeight{i} = CropNewDEM -(timeDelta*yearlyChange);
% %             DEM{i}.H = NewSurfaceHeight{i} - CropBed;
% %         end
%          if year(Tdates(i))>=LatestDEMYear
%             DEM{i}.X = xDEM;
%             DEM{i}.Y = yDEM;
%             DEM{i}.H = CropNewDEM - CropBed;
%             DEM_err(i,1) = arc_DEM_error;
%             DEM_err(i,2) = arc_DEM_error;
%          end
%     end
%    DEM_error = DEM_err; 
%    DEMstruct = DEM;
%    ChangeStruct = Change;
% % 
% end