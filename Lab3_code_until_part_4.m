clear
close all
clc

[data txt raw] = xlsread('lab3_data.xlsx');

date = data(:,1);
ubctanom = data(:,2); % Celsius, UBC temp anomaly
gtanom = data(:,3); % Celsiu, global temp anomaly
TSI = data(:,4); % W/m^2, solar irradiance at top of atmosphere
AOD = data(:,5); % Aerosol Optical Depth
CO2 = data(:,6); % ppm, atmospheric CO2 concentration
SO2 = data(:,7); % Tg/year, anthropogenic SO2 emissions
MEI = data(:,8); % Multivariate El Nino Index

%% Part 1: Plot each time series and linear regression of Temp Series
figure('pos',[30 30 900 600])
subplot(421)
    plot(date,ubctanom,'k'); hold on
        xlabel('Date'); ylabel({'Temp ' ,'Anomaly ({\circ}C)'})
    title({'UBC Temperature Anomaly' ,'from 1959-2016'})
    axis([1955,2018,-5,6])
    mubc = ~isnan(ubctanom);
    [coef,bint,r,rint,stats] = ...
        regress(ubctanom,[ones(size(ubctanom)) date]);
    ubctanomlinfit = coef(1)+coef(2).*date;
    ubctanomlinfitbotlim = coef(1)+bint(2,1).*date;
    ubctanomlinfittoplim = coef(1)+bint(2,2).*date;
    plot(date(mubc),ubctanomlinfit(mubc),'r'); hold on
%     plot(date(mubc),ubctanomlinfitbotlim(mubc),'r--'); hold on
%     plot(date(mubc),ubctanomlinfittoplim(mubc),'r--'); hold on
    text(1956,5.2,['Slope = ',num2str(round(coef(2),3)),...
        '; R^2=',num2str(round(stats(1),3))...
        '; conf = [' num2str(round(bint(2,1),3)) ','...
            num2str(round(bint(2,2),3)) ']'])
        
subplot(422)
    plot(date,gtanom,'k'); hold on
    xlabel('Date'); ylabel({'Temp ' ,'Anomaly (C^{\circ})'})
    title({'Global Temperature Anomaly' ,' from 1950-2016'})
    axis([1950,2018,-1,1.5])
      mgt = ~isnan(gtanom);
      clear bint
    [coef,bint,r,rint,stats] = ...
        regress(gtanom,[ones(size(gtanom)),date]);
    gtanomlinfit = coef(1)+coef(2).*date;
    gtanomlinfitbotlim = coef(1)+bint(2,1).*date;
    gtanomlinfittoplim = coef(1)+bint(2,2).*date;
    plot(date,gtanomlinfit,'r--'); hold on
%     plot(date,gtanomlinfitbotlim,'r--'); hold on
%     plot(date,gtanomlinfittoplim,'r--'); hold on
    text(1951,1.2,['Slope = ',num2str(round(coef(2),3)),...
        '; R^2=',num2str(round(stats(1),3))...
        '; conf = [' num2str(round(bint(2,1),3)) ','...
            num2str(round(bint(2,2),3)) ']']) 
subplot(423)
    plot(date,TSI,'r'); xlabel('Date');
    ylabel({'Total Solar ' ,'Irradiance (W/m^2)'})
    title({'Total Solar ' ,'Irradiance from 1950-2016'})
subplot(424)
    plot(date,AOD,'b'); xlabel('Date');
    ylabel('AOD')
    title({'Aerosol Optical' ,' Depth from 1950-2016'})
subplot(425)
    plot(date,CO2,'r'); xlabel('Date');
    ylabel('CO_2 (ppm)')
    title({'Atmospheric CO2 ' ,'from 1950-2016'})
subplot(426)
    plot(date,SO2,'b'); xlabel('Date');
    ylabel('SO_2 (Tg/year)')
    title({'Anthropogenic Atmospheric ' ,'SO_2 from 1950-2016'})
subplot(427)
    plot(date,MEI,'g'); xlabel('Date');
    ylabel('MEI')
    title({'Multivariate El Niño' ,' Index from 1950-2016'})
    
    
% Histogram Plots

m1 = date <= 1985;
m2 = date > 1985;
numbins = 40;
figure(2)
subplot(211)
histogram(ubctanom(m1),linspace(min(ubctanom),max(ubctanom),numbins),...
    'normalization','probability'); hold on
    xlabel('Temp Anomaly ({\circ}C)'); ylabel('Count per bin')
    title(['Histogram of UBC Temperature Anomalies, binsize = ',...
        num2str(numbins)])
histogram(ubctanom(m2),linspace(min(ubctanom),max(ubctanom),numbins),...
    'normalization','probability');
legend('<1985','>1985')
subplot(212)
histogram(gtanom(m1),linspace(min(gtanom),max(gtanom),numbins),...
    'normalization','probability'); hold on
    xlabel('Temp Anomaly ({\circ}C)'); ylabel('Count per bin')
    title(['Histogram of Global Temperature Anomalies, binsize = ',...
        num2str(numbins)])
histogram(gtanom(m2),linspace(min(gtanom),max(gtanom),numbins),...
    'normalization','probability');
legend('<=1985','>1985')

%% Part 2: Decadal Timescale Trends

Table1Global=zeros(7,4);
%columns 1,2,3,4,5 are decade (start year), slope,CI(min),CI(max) respectively
figure(3)
n = 1;
for i = 1960:10:2020
    subplot(4,2,n)
    mMEIandganom = date>=i-10 & date<i;
    [coef,bint,r,rint,stats] = ...
            regress(gtanom(mMEIandganom),[ones(size(gtanom(mMEIandganom))),date(mMEIandganom)]);
        gtanomlinfit = coef(1)+coef(2).*date(mMEIandganom);
        plot(date(mMEIandganom),gtanom(mMEIandganom),'k-'); hold on
        plot(date(mMEIandganom),gtanomlinfit,'r-'); hold on
        text(i-10,max(gtanom(mMEIandganom)),['Slope = ',num2str(round(coef(2),3)),...
            ' {\circ}C/yr; R^2=',num2str(round(stats(1),3))...
            '; 95% conf = [' num2str(round(bint(2,1),3)) ','...
            num2str(round(bint(2,2),3)) ']']) 
        xlabel('Date'); ylabel('Temp Anomaly ({\circ}C)')
        title(['Global Temp Anomaly for ' num2str(i-10) '-' num2str(i)])
        Table1Global(n,1)=i-10;
        Table1Global(n,2)=coef(2);
        Table1Global(n,3)=bint(2,1);
        Table1Global(n,4)=bint(2,2);
        
        n = n+1;
end


Table1UBC=zeros(7,4);
%columns 1,2,3,4,5 are decade (start year), slope,CI(min),CI(max) respectively
figure(4)
    n = 1;
for i = 1960:10:2020
    subplot(4,2,n)
    mMEIandganom = date>=i-10 & date<i;
    [coef,bint,r,rint,stats] = ...
            regress(ubctanom(mMEIandganom),[ones(size(ubctanom(mMEIandganom))),date(mMEIandganom)]);
        ubctanomlinfit = coef(1)+coef(2).*date(mMEIandganom);
        plot(date(mMEIandganom),ubctanom(mMEIandganom),'k-'); hold on
        plot(date(mMEIandganom),ubctanomlinfit,'r-'); hold on
        text(i-10,max(ubctanom(mMEIandganom)),['Slope = ',num2str(round(coef(2),3)),...
            ' {\circ}C/yr; R^2=',num2str(round(stats(1),3))...
            '; 95% conf = [' num2str(round(bint(2,1),3)) ','...
            num2str(round(bint(2,2),3)) ']']) 
        xlabel('Date'); ylabel('Temp Anomaly ({\circ}C)')
        title(['UBC Temp Anomaly for ' num2str(i-10) '-' num2str(i)])
        Table1UBC(n,1)=i-10;
        Table1UBC(n,2)=coef(2);
        Table1UBC(n,3)=bint(2,1);
        Table1UBC(n,4)=bint(2,2);
        n = n+1;
end


%% Part 3: Local vs. Global Temp

figure(5)
plot(gtanom,ubctanom,'k.'); hold on
xlabel('UBC Temp ({\circ}C)'); ylabel('Global Temp ({\circ}C)')
mMEIandganom = isnan(ubctanom) | isnan(gtanom);
[R,p] = corrcoef(ubctanom(~mMEIandganom),gtanom(~mMEIandganom));
[coef,bint,r,rint,stats] = regress(ubctanom(~mMEIandganom),[ones(size(gtanom(~mMEIandganom))) gtanom(~mMEIandganom)]);
R1 = stats(1)^2; p1 = stats(3);
text(-0.4,-5,['Correlation = ' num2str(R(2,1)) ' and p-value = ' num2str(p(2,1))]);


%% Part 4: Impact of specific forcingon global temperature anomaly

AllForcingsMatrix=[TSI,AOD,CO2,SO2,MEI];
Table2SimpleLinear = zeros(5,4);
%rows 1,2,3,4,5 are TSI,AOD,CO2,SO2,MEI respectively
%columns 1,2,3,4 are slope, CI,CI,coeffecient of determination respectively
n = 1;
for forcing = AllForcingsMatrix
    mindependentForcing = isnan(forcing) | isnan(gtanom);
    [coef,bint,r,rint,stats] = regress(gtanom(~mindependentForcing),[ones(size(forcing(~mindependentForcing))) forcing(~mindependentForcing)]);
    Table2SimpleLinear(n,1)= coef(2);
    Table2SimpleLinear(n,2)=bint(2,1);
    Table2SimpleLinear(n,3)=bint(2,2);
    Table2SimpleLinear(n,4)=stats(1);

    n=n+1;
end

figure(6)
subplot(2,1,1)
mMEIandganom = ~(isnan(MEI) | isnan(gtanom));
hold on
plot(MEI(mMEIandganom),gtanom(mMEIandganom),'k.'); 
xlabel('MEI');
ylabel('Global Temp Anomaly({\circ}C)');
title('Global Temp Anomaly VS MEI');
clear bint
[coef,bint,r,rint,stats] = regress(gtanom(mMEIandganom),[ones(size(gtanom(mMEIandganom))),MEI(mMEIandganom)]);
glMEIlinfit = coef(1)+coef(2).*MEI(mMEIandganom);
glMEIlinfitbotlim = coef(1)+bint(2,1).*MEI(mMEIandganom);
glMEIlinfittoplim = coef(1)+bint(2,2).*MEI(mMEIandganom);
plot(MEI(mMEIandganom),glMEIlinfit,'r--'); 
%plot(MEI(mMEIandganom),glMEIlinfitbotlim,'b--'); hold on
%plot(MEI(mMEIandganom),glMEIlinfittoplim,'r--'); hold on
axis([-2.5,3.5,-0.7,1.3])
text(-2.4,1,['Slope = ',num2str(round(coef(2),3))]) 
hold off


subplot(2,1,2)
mMEIandubcanom = ~(isnan(MEI) | isnan(ubctanom));
hold on
plot(MEI(mMEIandubcanom),ubctanom(mMEIandubcanom),'k.'); 
xlabel('MEI');
ylabel('UBC Temp Anomaly({\circ}C)');
title('UBC Temp Anomaly VS MEI');
clear bint
clear coef
[coef,bint,r,rint,stats] = regress(ubctanom(mMEIandubcanom),[ones(size(ubctanom(mMEIandubcanom))),MEI(mMEIandubcanom)]);
UBCMEIlinfit = coef(1)+coef(2).*MEI(mMEIandubcanom);
UBCMEIlinfitbotlim = coef(1)+bint(2,1).*MEI(mMEIandubcanom);
UBCMEIlinfittoplim = coef(1)+bint(2,2).*MEI(mMEIandubcanom);
plot(MEI(mMEIandubcanom),UBCMEIlinfit,'r--'); 
%plot(MEI(mMEIandubcanom),UBCMEIlinfitbotlim,'b--'); hold on
%plot(MEI(mMEIandubcanom),UBCMEIlinfittoplim,'r--'); hold on
axis([-2.2,3.2,-6,6])
text(-2.1,5,['Slope = ',num2str(round(coef(2),3))]) 
hold off

%MULTILINEAR REGRESSION
Table2Multi = zeros(7,3);
clear bint
clear coef
%rows 2,3,4,5,6 are TSI,AOD,CO2,SO2,MEI respectively
%row 1 is y intercept
%row 7 is R^2
%columns 1,2,3 are slope, CI,CI respectively


all_forcings_mask = ~(isnan(AOD) |isnan(CO2) |isnan(SO2) |isnan(MEI) |isnan(TSI) | isnan(gtanom));
[coef,bint,r,rint,stats] = regress(gtanom(all_forcings_mask),[ones(size(TSI(all_forcings_mask))) TSI(all_forcings_mask) AOD(all_forcings_mask) CO2(all_forcings_mask) SO2(all_forcings_mask) MEI(all_forcings_mask) ]); 
Table2Multi(1,1)=coef(1);
Table2Multi(2,1)=coef(2);
Table2Multi(3,1)=coef(3);
Table2Multi(4,1)=coef(4);
Table2Multi(5,1)=coef(5);
Table2Multi(6,1)=coef(6);

Table2Multi(1,2)=bint(1,1);
Table2Multi(1,3)=bint(1,2);
Table2Multi(2,2)=bint(2,1);
Table2Multi(2,3)=bint(2,2);
Table2Multi(3,2)=bint(3,1);
Table2Multi(3,3)=bint(3,2);
Table2Multi(4,2)=bint(4,1);
Table2Multi(4,3)=bint(4,2);
Table2Multi(5,2)=bint(5,1);
Table2Multi(5,3)=bint(5,2);
Table2Multi(6,2)=bint(6,1);
Table2Multi(6,3)=bint(6,2);

Table2Multi(7,1)=stats(1);
