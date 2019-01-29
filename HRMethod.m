
% load the data
close all
clear all
clc

load sig4

%vcg2 = signals;
vcg2 = signals(:,200000:260000);
fs = 1000;
ycomb = zeros(1,length(vcg2));
ycomb = sqrt(vcg2(1,:).^2 + vcg2(2,:).^2 + vcg2(3,:).^2); % sqrt(xˆ2 + yˆ2 + zˆ2)
figure;plot(ycomb);
[qrs_amp_raw,qrs_i_raw,delay]=pan_tompkin(double(ycomb), 1000, 1);
i_R_peaks = qrs_i_raw;
plot(ycomb)
hold on
plot(qrs_i_raw,ycomb(1,qrs_i_raw),'*r')
hold off

plot(ycomb)
hold on
plot(qrs_i_raw,ycomb(1,qrs_i_raw),'*r')
hold off

% Let's select the first template 
iBefore= 183;  %round(180/1000*128); 
iAfter= 400; %round(600/1000*128); 
Ref_delta = 20;
if i_R_peaks(1)-iBefore-Ref_delta < 0
    template=ycomb(i_R_peaks(2)-iBefore-Ref_delta:i_R_peaks(2)+iAfter+Ref_delta); 
    K_starter = 3;
else
    template=ycomb(i_R_peaks(1)-iBefore-Ref_delta:i_R_peaks(1)+iAfter+Ref_delta); 
    K_starter = 2;
end
plot(template) 

% Position of the second piece of signal which is very similar to the
% template 
startp = [];
endp = [];
timecounter = 1;
for k=K_starter:(K_starter+8)

    iRange=20; 
    maxCorr=0; 
    iMaxCorr=NaN; 

    for i=-iRange:iRange 
        thisCorr=sum(template.*ycomb(i_R_peaks(k)-iBefore-Ref_delta-i:i_R_peaks(k)+iAfter+Ref_delta-i)); 
        %thisCorr = xcorr(template,ycomb(i_R_peaks(2)-iBefore-i:i_R_peaks(2)+iAfter-i));
        if thisCorr>maxCorr 
            maxCorr=thisCorr; 
            iMaxCorr=i; 
        end 
    end 
    startp(timecounter) = i_R_peaks(k)-iBefore-iMaxCorr;
    endp(timecounter) = i_R_peaks(k)+iAfter-iMaxCorr;
    timecounter = timecounter + 1;
   % template = (template + ycomb(i_R_peaks(2)-iBefore-iMaxCorr:i_R_peaks(2)+iAfter-iMaxCorr))/2;
end

startTemplatePoint= i_R_peaks(K_starter-1)-iBefore-Ref_delta;
endTemplatePoint = i_R_peaks(K_starter-1)+iAfter+Ref_delta;
Yref = [zeros(3,endTemplatePoint-startTemplatePoint + 1)];
%find Yref x
sumOfP = vcg2(1,startTemplatePoint:endTemplatePoint) + vcg2(1,startp(1)-Ref_delta:endp(1)+Ref_delta);
for p=2:9
  sumOfP  = sumOfP + vcg2(1,startp(p)-Ref_delta:endp(p)+Ref_delta);
end

Yref(1,:) = sumOfP/10;

%find Yref y
sumOfP = vcg2(2,startTemplatePoint:endTemplatePoint) + vcg2(2,startp(1)-Ref_delta:endp(1)+Ref_delta);

for p=2:9
  sumOfP  = sumOfP + vcg2(2,startp(p)-Ref_delta:endp(p)+Ref_delta);
end
Yref(2,:) = sumOfP/10;

%find Yref z
sumOfP = vcg2(3,startTemplatePoint:endTemplatePoint) + vcg2(3,startp(1)-Ref_delta:endp(1)+Ref_delta);

for p=2:9
  sumOfP  = sumOfP + vcg2(3,startp(p)-Ref_delta:endp(p)+Ref_delta);
end

Yref(3,:) = sumOfP/10;

% plot Yref
figure;plot3(Yref(1,:), Yref(2,:), Yref(3,:));


delta = 20; %based on paper
E_Yref = Yref;
angleX = [];
angleY = [];
angleZ = [];
i_ref_peak= [];
for p = 1:size(i_R_peaks,2)-1
    [Final_Q, Ts] = OperatorSplitting(vcg2(1:3,:), E_Yref, delta, i_R_peaks,iBefore,iAfter,p, size(Yref,2), K_starter);
    i_ref_peak(p) = i_R_peaks(p) - Ts;
    fi_y = asin(Final_Q(1,3));
    fi_x = asin(Final_Q(2,3)/cos(fi_y));
    fi_z = asin(Final_Q(1,2)/cos(fi_y));
    angleX(p) = fi_x;
    angleY(p) = fi_y;
    angleZ(p) = fi_z;
  
end
%%
%plot Angles
figure;plot(angleX);
figure;plot(angleY);
figure;plot(angleZ);
%% 
q = i_ref_peak*(1/fs); %should 
w = angleX; 
xq1 = q(1):0.2:q(end); %1000 is sampling rate
angleX_resampled = spline(q,w,xq1);
figure;plot(q,w,'o',xq1,angleX_resampled);

q = i_ref_peak*(1/fs); 
w = angleY; 
xq1 = q(1):0.2:q(end);
angleY_resampled = spline(q,w,xq1);
figure;plot(q,w,'o',xq1,angleY_resampled,'-.');

q = i_ref_peak*(1/fs); 
w = angleZ; 
xq1 = q(1):0.2:q(end);
angleZ_resampled = spline(q,w,xq1);
figure;plot(q,w,'o',xq1,angleZ_resampled,'-.');
% %%
% % PSD
% 
%%
L = length(angleX_resampled) ;

X_FourierT=fft(angleX_resampled);
X_FourierT = abs(X_FourierT);
X_FourierT = X_FourierT.^2;
N=length(angleX_resampled);
Sx = X_FourierT/N;

% figure;plot(fft(angleX_resampled));
% 
Y_FourierT=fft(angleY_resampled);
Y_FourierT = abs(Y_FourierT);
Y_FourierT = Y_FourierT.^2;
N=length(angleY_resampled);
Sy = Y_FourierT/N;
% figure;plot(Sy);
% 
Z_FourierT=fft(angleZ_resampled);
Z_FourierT = abs(Z_FourierT);
Z_FourierT = Z_FourierT.^2;
N=length(angleZ_resampled);
Sz = Z_FourierT/N;
% figure;plot(Sz);%hold on; plot(0.07,0.5,'bo');
% 
% %%
%averaging PSD
averageOfThreePSD = (Sx + Sy + Sz)/3;

%% Cubic Spline Interpolation of heart rate Signal 

% q = i_ref_peak*(1/fs); %should 
% w = angleX; 
% xq1 = q(1):0.2:q(end); %1000 is sampling rate
% angleX_resampled = spline(q,w,xq1);
% figure;plot(q,w,'o',xq1,angleX_resampled);


hr=diff(i_ref_peak);              %no of samples in peaks 
figure;plot(hr);
q = i_ref_peak *(1/fs);  

w = hr;            %peaks 

xq1 = q(1):0.2:q(end);        %1000 is sampling rate 

splHr = spline(q(1:end-1),w,xq1); 

figure;plot(q(1:end-1),w,'.',xq1,splHr); 

 

%% Power Spectrum of Heart Rate 

hrFFT=fft(splHr); 
figure;plot(hrFFT);
lensplHR= length(splHr); 

hrPSD =  (abs(hrFFT).^2)/lensplHR; 
figure;plot(hrPSD)
%%
conj_FFT_X = conj(fft(angleX_resampled));
conj_FFT_Y = conj(fft(angleY_resampled));
conj_FFT_Z = conj(fft(angleY_resampled));

NN = length(conj_FFT_Z)
new_conj_x = (hrFFT.*conj_FFT_X)/NN;
new_conj_y = (hrFFT.*conj_FFT_Y)/NN;
new_conj_z = (hrFFT.*conj_FFT_Z)/NN;
%%
ave = (abs(new_conj_x) + abs(new_conj_y) + abs(new_conj_z))/3;
EE = length(ave);
f = 5*(0:EE)/EE;
% 
figure;semilogy(f(f>0.07 & f< 0.5), ave(f>0.07 & f< 0.5))

[ttt Imasx] = max(ave(f>0.07 & f< 0.5));
result = find(ave == ttt);
our_frequency_third = f(result(1,1));
Respiratory_frequency_by_VCG = round(our_frequency_third*60)

Resp_signal = vcg2(4,:);

Resp_signal_FourierT=fft(Resp_signal);
Resp_signal_FourierT = Resp_signal_FourierT.^2;
N_resp=length(Resp_signal);
Resp_signal_FourierT = abs(Resp_signal_FourierT);
S_resp = Resp_signal_FourierT/N_resp;

f_resp = 1000*(0:N_resp)/N_resp;
figure;plot(f_resp(1:end-1), S_resp);
[Resp_Max I_resp] = max(S_resp);
Resp_frequecy = f_resp(I_resp);
Respiratory_frequency_by_Resp_Signal = round(Resp_frequecy*60)


Perf_EDR = (abs(Resp_frequecy - our_frequency_third)/Resp_frequecy)*100


