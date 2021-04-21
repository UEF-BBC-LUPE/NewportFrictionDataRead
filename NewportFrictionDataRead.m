%% New Newport lubrication / friction data read 
clear all;
close all;
polku = uigetdir; % Pick the file containing measurement data
cd(polku)
list1=dir([pwd filesep '*.txt'])
% list1=dir([pwd filesep '*.xlsx'])
visc=0.89*1e-3; % water viscosity
% visc=81.42*1e-3; % olive oil measured with Rheosense
% visc=397.3571*1e-3; % Belray V-twin 20w50 motor oil viscosity.
%%
for II=1:length(list1)
%     measurement(II)=importdata([path filesep list1(II).name]);
filename=list1(II).name;
measurement(II)=importdata(filename);
    
    % measurement(II).textdata contains the GUI parameters (first row of
    % .txt file
    % initial z-velocity, adjustment velocity, speed decr rate,data read
    % interval (ms), target force, limit + (upper), limit � (lower), nmb of steps
    % step time (s), lift off (um), lift off period (s), sample rate (Hz),
    % file write averaging, sample averaging, ADC maximum value, ADC min value
% variables in the measuremnt(II).data:
% time program; time step; force; torque; position �m; rotation speed rps;
% adjustment velocity mm/s 
end

%%
close all
tt=[-5 360];
for ii=1:length(list1)
    % Read measurement
%     parameters from text file
    GUI=str2num(cell2mat(measurement(ii).textdata)); 
    F_trgt=GUI(5);
    init_vel=GUI(1); % initial velocity before contact
    adj_vel=GUI(2); % adjustment velocity during creep
    steptime=measurement(ii).data(:,1); % step time
    time=measurement(ii).data(:,2); % time for the whole measurement
    fs(ii)=length(time)/(time(end)-time(1)); % calculate sampling frequency
    Force=measurement(ii).data(:,3); % force N
    Force=Force-Force(1,1); % zero the force
    torque=measurement(ii).data(:,4); % torque Nm
    torque=torque-torque(1,1); % zero the torque 
    position=measurement(ii).data(:,5); % position in �m
    v=measurement(ii).data(:,6); % rotation velocity rot/s
    vcomp=measurement(ii).data(:,7); % compression velocity mm/s 
    
    % Two ways to define the time points for analysis
    % 1) Find time points were we have contact based on compression velocity
%     [idx1,idx2]=find(vcomp<init_vel);
%     [idx1,idx2]=find(vcomp<0.5*init_vel);
%     [idx1,idx2]=find(vcomp<5*adj_vel);
%     [idx1,idx2]=find(vcomp<2*adj_vel);

    % 2) Find time points where we have reached 95% of the target force
    [idx1,idx2]=find(Force>0.95*F_trgt); 
    
    figure % plot raw data and the chosen time point for analysis
    subplot(2,4,1)
    plot(time,Force,'--')
    hold on
    plot(time(idx1),Force(idx1),'o')
    xlabel('time (s)')
    ylabel('force (N)')
    xlim(tt)
    subplot(2,4,2)
    plot(time,position,'--')
        hold on
    plot(time(idx1),position(idx1),'o')
    xlabel('time (s)')
    ylabel('position (mm)')
    xlim(tt)
    d=5.2;
    ro=d/2;
    ri=1.0;
    Reff=2/3*ro*1.0e-3;
%     Reff=2/3*((ro*1e-3)^3-(ri*1.0e-3)^3)/((ro*1e-3)^2-(ri*1.0e-3)^2); %
%     For annulus geometry
    mu=abs(torque)./(abs(Force)*Reff);
    
    H=(visc*v*2*pi*Reff)./Force;
    subplot(2,4,3)
    plot(time,torque,'--')
        hold on
    plot(time(idx1),torque(idx1),'o')
    xlabel('time(s)')
    ylabel('torque (Nm)')
    xlim(tt)
    subplot(2,4,4)
    plot(time,mu,'--')
    hold on
    plot(time(idx1),mu(idx1),'o')
    xlabel('time(s)')
    ylabel('COF ')
    xlim(tt)
    subplot(2,4,5)
    plot(time,H,'--')
    hold on
    plot(time(idx1),H(idx1),'o')
    xlabel('time(s)')
    ylabel('Hersey')
    xlim(tt)
    subplot(2,4,6)
    plot(time,vcomp,'--')
    hold on
    plot(time(idx1),vcomp(idx1),'o')
    xlabel('time(s)')
    ylabel('comp vel')
    xlim(tt)
    subplot(2,4,7)
    plot(time(idx1),v(idx1))
    xlabel('time(s)')
    ylabel('rot vel')
    xlim(tt)
    Hersey(ii).data=H(idx1);
    COF(ii).data=mu(idx1);
end
%% plot semilogaritmic Stribeck curves

figure
for jj=[1:10] % edit these indexes according to the data to correspond to velocity etc.
    semilogx(Hersey(jj).data(:),COF(jj).data(:),'ob')
    hold on
end
%
for jj=[11:20]
    semilogx(Hersey(jj).data(:),COF(jj).data(:),'or')
end
%%
for jj=[33:48]
    semilogx(Hersey(jj).data(:),COF(jj).data(:),'ok')
end
xlabel('Hersey')
ylabel('COF')

%% plot logaritmic Stribeck curves
figure
for jj=[1:16]
loglog(Hersey(jj).data(:),COF(jj).data(:),'ob')
hold on
end
for jj=[17:32]
loglog(Hersey(jj).data(:),COF(jj).data(:),'or')
end
for jj=[33:48]
loglog(Hersey(jj).data(:),COF(jj).data(:),'ok')
end
xlabel('Hersey')
ylabel('COF')

%     figure
%     histogram(fs)
%% Next part implements filtering to the data - NOT READY, NEED REVISION/EDIT
% vm=mean(v);
%     wp=(vm)/(2*pi);
%     wpass=[0.05*wp 0.9*wp];
%     ws = 100; 
%     b = (1/ws)*ones(1,ws); % Numerator coefficients
%     a = 1; % Denominator coefficients


% for ii=1:length(list1)
%     figure
% %     GUI=str2num(cell2mat(measurement(ii).textdata));
% %     F_trgt=GUI(5);
%     v=str2double(measurement(ii).data(2:end,6)); % implement band-stop filter to remove noise due to rotation
%     vm=mean(v);
%     wp=(vm)/(2*pi);
%     wpass=[0.05*wp 0.9*wp];
%     ws = 100; 
%     b = (1/ws)*ones(1,ws); % Numerator coefficients
%     a = 1; % Denominator coefficients
% 
%     time=str2double(measurement(ii).data(2:end,2));
%     Force=str2double(measurement(ii).data(2:end,3));
%     Force=Force+abs(Force(1,1)); % Force is in Newtons, but the zero load reading has to be added 
%     Force2=filtfilt(b,a,Force);
% %     Force3=Force-bandstop(Force,wpass);
%     lbintoNm=0.11298482933333;
% %     lbintoNm=1;
%     torque=lbintoNm*str2double(measurement(ii).data(2:end,4)); %Torque is in Newton meters(?)
%         torque=torque-abs(torque(1,1));
%     torque2=filtfilt(b,a,torque);
%     position=str2double(measurement(ii).data(2:end,4));
%     position=position-abs(position(1,1));
%     position2=filtfilt(b,a,position);
%     
%     subplot(2,3,1)
%     plot(time,Force,'--')
%     hold on
%     plot(time,Force2,'-r')
%     xlabel('time (s)')
%     ylabel('force (F)')
%     xlim([5 300])
%     subplot(2,3,2)
%     plot(time,abs(position),'--')
%     hold on
%     plot(time,abs(position2),'-r')
%     xlabel('time (s)')
%     ylabel('position (mm)')
%     xlim([5 300])
%     
% %     ind1=find(Force>0.02,1,'first');
% %     ind2=find(Force>F_trgt,1,'first');
% %     t_comp(ii)=time(ind2);
% %     disp(ii)=(position(ind2)-position(ind1))*1e-3;
% %     F_step(ii)=F_trgt;
%     stress=Force/(3.3^2*pi);
%     strain=abs(position)/10;
%     E=(stress)./(strain);
%     ind=find(strain==max(strain),1,'last');
%     stress_max(ii)=stress(ind);
%     strain_max(ii)=strain(ind);
% %     Reff=2/3*3.3e-3;
%     Ro=3.3; % in millimeters
%     Ri=0;
%     Reff=2/3*((Ro*1.0e-3)^3-(Ri*1.0e-3)^3)/((Ro*1.0e-3)^2-(Ri*1.0e-3)^2) % dimensions now in meters!
%     mu=abs(torque)./(abs(Force)*Reff);
%     mu2=abs(torque2)./(abs(Force2)*Reff);
%     visc=0.00089;
%     H=(visc*v*2*pi*Reff)./Force;
%     H2=(visc*v*2*pi*Reff)./Force2;
%     subplot(2,3,3)
%     plot(time,abs(torque),'--')
%     hold on
%     plot(time,abs(torque2),'-r')
%     xlabel('time(s)')
%     ylabel('torque (Nm)')
%     xlim([5 300])
%     subplot(2,3,4)
%     plot(time,mu,'--')
%     hold on 
%     plot(time,mu2,'--r')
%     xlabel('time(s)')
%     ylabel('COF ')
%     xlim([5 300])
%     subplot(2,3,5)
%     plot(time,H,'--')
%     hold on
%     plot(time,H2,'-r')
%     xlabel('time(s)')
%     ylabel('Hersey')
%     xlim([5 300])
%     ind1=find(time<=250,1,'last');
%     ind2=find(time<=290,1,'last');
%     COF(ii)=mean(mu2(ind1:ind2));
%     Hersey(ii)=mean(H2(ind1:ind2));
% end
% 
% figure
% loglog(Hersey,COF,'o')
% xlabel('Hersey')
% ylabel('COF')