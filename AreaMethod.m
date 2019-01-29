
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
plot(ycomb)
hold on
plot(qrs_i_raw,ycomb(1,qrs_i_raw),'*r')
hold off

plot(ycomb)
hold on
plot(qrs_i_raw,ycomb(1,qrs_i_raw),'*r')
hold off
i_R_peaks = qrs_i_raw;
% Let's select the first template 
iBefore= 183;  %round(180/1000*128); 
iAfter= 400; %round(600/1000*128); 
YYY=[];

for o=1:size(i_R_peaks,2)
    if i_R_peaks(o) - iBefore < 0
        continue;
    end
    if i_R_peaks(o) + iAfter > size(vcg2,2)
        signal_A = vcg2(1,i_R_peaks(o) - iBefore : size(vcg2,2));
        signal_B = vcg2(2,i_R_peaks(o) - iBefore : size(vcg2,2));
    else
        signal_A = vcg2(1,i_R_peaks(o) - iBefore : i_R_peaks(o)+iAfter);
        signal_B = vcg2(2,i_R_peaks(o) - iBefore : i_R_peaks(o)+iAfter);
    end
    


area_A = trapz(signal_A);
area_B = trapz(signal_B);

YYY(o) = atan(area_A/area_B)

end
L = length(YYY);
YYY_FourierT=fft(YYY); 

YYY_FourierT = abs(YYY_FourierT); 

YYY_FourierT = YYY_FourierT.^2; 

N=length(YYY_FourierT); 


another_Method = YYY_FourierT/N; 

f = 5*(0:L)/L;
% 
figure;plot(f(1:end-1), another_Method)
[The_Max I] = max(another_Method(f<0.5 & f>0.07));

result_another = find(another_Method == The_Max);
our_frequency = f(result_another(1,1));
Respiratory_frequency_by_VCG = round(our_frequency*60)

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

Perf_EDR = (abs(Resp_frequecy - our_frequency)/Resp_frequecy)*100

