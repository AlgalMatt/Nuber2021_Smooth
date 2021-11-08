%This script imports data from the named file and performs
%a gaussian smooth to produce a single mean record for each variable.

%The data a first linearly interpolated to mitigate sampling bias. 

clear all
close all

%% Data import
filename='Matlab_ready_SSS.xlsx'; %import file
opts = detectImportOptions(filename);
sheets = sheetnames(filename);
xend=1200; %xaxis limits in kyr (shorten to speed up script)
window=8; %smoothing window size in kyr
sigma=window/5;
%extract data
  for i=1:numel(sheets)
  dat=readtable(filename, 'Sheet',sheets(i));
  %dat.Age(isnan(dat.SST),:)=[]; 
  data{i}=dat; 
  end


%% interpolation

for i=1:numel(sheets)
    %remove NaNs from data
    SST{i}=[data{i}.Age(~isnan(data{i}.SST)),data{i}.SST(~isnan(data{i}.SST))];
    d18O{i}=[data{i}.d18Oage(~isnan(data{i}.d18O)),data{i}.d18O(~isnan(data{i}.d18O))];
    SSS{i}=[data{i}.SSSage(~isnan(data{i}.SSS)),data{i}.SSS(~isnan(data{i}.SSS))];
    %z scores
    SST{i}(:,3)=(SST{i}(:,2)-mean(SST{i}(:,2), 'omitnan'))./std(SST{i}(:,2), 'omitnan');
    d18O{i}(:,3)=(d18O{i}(:,2)-mean(d18O{i}(:,2), 'omitnan'))./std(d18O{i}(:,2), 'omitnan');
    SSS{i}(:,3)=(SSS{i}(:,2)-mean(SSS{i}(:,2), 'omitnan'))./std(SSS{i}(:,2), 'omitnan');
end


 xTinterp=[]; Tinterp=[]; zTinterp=[];
 xOinterp=[];  Ointerp=[]; zOinterp=[];
 xSinterp=[];  Sinterp=[]; zSinterp=[];
 for i=1:numel(sheets)
     xTinterp = [xTinterp; [SST{i}(:,1):0.01:min([SST{i}(end,1),xend])]'];
     Tinterp=[Tinterp; interp1(SST{i}(:,1),SST{i}(:,2), [SST{i}(:,1):0.01:min([SST{i}(end,1),xend])]')];
     zTinterp=[zTinterp; interp1(SST{i}(:,1),SST{i}(:,3), [SST{i}(:,1):0.01:min([SST{i}(end,1),xend])]')];
     
     xOinterp = [xOinterp; [d18O{i}(:,1):0.01:min([d18O{i}(end,1),xend])]'];
     Ointerp=[Ointerp; interp1(d18O{i}(:,1),d18O{i}(:,2), [d18O{i}(:,1):0.01:min([d18O{i}(end,1),xend])]')];
     zOinterp=[zOinterp; interp1(d18O{i}(:,1),d18O{i}(:,3), [d18O{i}(:,1):0.01:min([d18O{i}(end,1),xend])]')];
     
     if ~isempty(SSS{i})
         xSinterp = [xSinterp; [SSS{i}(:,1):0.01:min([SSS{i}(end,1),xend])]'];
         Sinterp=[Sinterp; interp1(SSS{i}(:,1),SSS{i}(:,2), [SSS{i}(:,1):0.01:min([SSS{i}(end,1),xend])]')];
         zSinterp=[zSinterp; interp1(SSS{i}(:,1),SSS{i}(:,3), [SSS{i}(:,1):0.01:min([SSS{i}(end,1),xend])]')];
     end
 end

 allinterp={sortrows([xTinterp, Tinterp]), sortrows([xOinterp, Ointerp]), sortrows([xSinterp, Sinterp]), ...
    sortrows([xTinterp, zTinterp]), sortrows([xOinterp, zOinterp]), sortrows([xSinterp, zSinterp])};
 ysmooth=cell(1,numel(allinterp));
  xsmooth=0:0.1:xend;

%% Gaussian smoothing

for j = 1:numel(xsmooth)
      minx=xsmooth(j)-window/2;
    maxx=xsmooth(j)+window/2;    
    
for i =1:numel(allinterp)
    %Finds the limits of the window around the data point
    windowindx=find(allinterp{i}(:,1)>=minx & allinterp{i}(:,1)<=maxx);
    %normalize the age so that the data point is at zero (centre of
    %gaussian)
     relativex=allinterp{i}(windowindx, 1)-xsmooth(j);   
    weights=[];
    for k=1:numel(relativex)
        %gaussian smoothing weights
        weights(k)=(1/(2*pi*sigma^2)^0.5)*exp((-relativex(k)^2)/(2*sigma^2));
    end
    %normalise weights so they add up to 1
    norm_wts=weights'./sum(weights);
    %Apply weights to the d13C data to smooth
    ywindow=allinterp{i}(windowindx, 2);
    ysmooth{i}(j)=sum(norm_wts.*ywindow, 'omitnan');
end  
end

%% Data export

 count =1;
    for i=1:numel(sheets)
writematrix([SST{i}(:,1), SST{i}(:,3)],['allzdata_' date '.xlsx'],'Sheet',sheets(i),'Range','A2')
writematrix(["Age", "zSST"],['allzdata_' date '.xlsx'],'Sheet',sheets(i),'Range','A1')
writematrix([d18O{i}(:,1), d18O{i}(:,3)],['allzdata_' date '.xlsx'],'Sheet',sheets(i),'Range','C2')
writematrix(["Age", "zd18O"],['allzdata_' date '.xlsx'],'Sheet',sheets(i),'Range','C1')
if ~isempty(SSS{i})
writematrix([SSS{i}(:,1), SSS{i}(:,3)],['allzdata_' date '.xlsx'],'Sheet',sheets(i),'Range','E2')
writematrix(["Age", "zSSS"],['allzdata_' date '.xlsx'],'Sheet',sheets(i),'Range','E1')
end
count=count+1;
    end

writematrix(["Age", "SST", "d18O", "SSS"],['smooth_SSTd18OSSS_mean_' date '.xlsx'],'Range','A1')
writematrix([xsmooth', ysmooth{1}', ysmooth{2}', ysmooth{3}'],['smooth_SSTd18OSSS_mean_' date '.xlsx'],'Range','A2')
writematrix(["Age", "zSST", "zd18O", "zSSS"],['smooth_zSSTd18OSSS_mean_' date '.xlsx'],'Range','A1')
writematrix([xsmooth', ysmooth{4}', ysmooth{5}', ysmooth{6}'],['smooth_zSSTd18OSSS_mean_' date '.xlsx'],'Range','A2')
    

