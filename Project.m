clc;
clear all;

%%Reading First Signal%%
[FirstSignal,F1]= audioread('Short_BBCArabic2.wav'); %reading first signal

%%Reading Second Signal%%
[SecondSignal,F2]=audioread('Short_FM9090.wav'); %reading second signal

%%Changing signals from stereo to monophonic%%
First_Signal_Mono=mean(FirstSignal,2);
Second_Signal_Mono=mean(SecondSignal,2);



%%Upsampling the signals%%

First_Signal_Sampling=interp(First_Signal_Mono,12);
Second_Signal_Sampling=interp(Second_Signal_Mono,12);

%%New Sampling%%
F1=F1*12; %new sampling frequancy for first signal
TS1=1/F1; %Time sampling for first signal
N1=[-length(First_Signal_Sampling)/2:(length(First_Signal_Sampling)/2-1)]; %vector for the first signal
L1=length(First_Signal_Sampling); %length of first signal
f1=[-L1/2:L1/2-1]*(F1/L1); %vector for frequancy sampling for the first signal

F2=F2*12; %new sampling frequancy for Second signal
TS2=1/F2; %Time sampling for Second signal
N2=[-length(Second_Signal_Sampling)/2:(length(Second_Signal_Sampling)/2-1)];
L2=length(Second_Signal_Sampling); %length of second signal
f2=[-L2/2:L2/2-1]*(F2/L2);  %vector for frequancy sampling for the second signal

%%Modulation%%
FC=100000; %carrier frequancy
df=50000; %delta carrier frequancy 

Carrier1=cos(2*pi*FC*N1*TS1); %carrier for first signal
Carrier2=cos(2*pi*(FC+df)*N2*TS2); %carrier for second signal


First_Signal_Modulated=First_Signal_Sampling'.*Carrier1; %first signal after modulation
figure(1)
plot(f1,abs(fftshift(fft(First_Signal_Modulated)))) % Plotting first signal after modulation
title('BBC ARABIA after Modulation')
xlabel('Frequency in HZ' )
ylabel ('Amplitude')

Second_Signal_Modulated=Second_Signal_Sampling'.*Carrier2; %second signal after modulation
figure(2)
plot(f2,abs(fftshift(fft(Second_Signal_Modulated)))) % Plotting second signal after modulation
title('FM9090 after Modulation')
xlabel('Frequency in HZ' )
ylabel ('Amplitude')


%check which signal is the longest and padding the short one with zeros
if(L1 > L2)
Second_Signal_Modulated=[Second_Signal_Modulated, zeros(1,L1-L2)];
elseif(L1 < L2 )
First_Signal_Modulated=[First_Signal_Modulated, zeros(1,L2-L1)];
end

Modulated_Signal=First_Signal_Modulated+Second_Signal_Modulated; %transmited signal
figure(3)
plot(f1,abs(fftshift(fft(Modulated_Signal)))) %plotting transmited signal
title('The 2 signals together')
xlabel('Frequency in HZ' )
ylabel ('Amplitude')

%%%%%%%%Getting_First_Signal%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%_RF STAGE_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BW1=10000; %bandwight of first signal
%bandpass signal with
% F_stop1= 90000 khz
% F_pass1= 95000 khz
% F_pass2= 105000 khz
% F_stop2= 110000 khz
RF=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',(FC-BW1),(FC-BW1/2),(FC+BW1/2),(FC+BW1),60,1,60,F1);
RF1=design(RF,'butter'); 
RF_OUT=filter(RF1,Modulated_Signal);
figure(4)
plot(f1,abs(fftshift(fft(RF_OUT))));     
title('BBC ARABIA after RF')
xlabel('Frequency in HZ' )
ylabel ('Amplitude')

%%%%%%%%%%%%%%%%%%%%%%%%%%_Mixer Stage_%%%%%%%%%%%%%%%%%%
F_IF=25000; % IF frequancy
F_OS=FC+F_IF; %oscillator frequancy
Mixer_Carrier=cos(2*pi*F_OS*N1*TS1);
Mixed_Signal=Mixer_Carrier.*RF_OUT;
figure(5)
plot(f1,abs(fftshift(fft(Mixed_Signal))))
title('BBC ARABIA after Mixer')
xlabel('Frequency in HZ' )
ylabel ('Amplitude')

%%%%%%%%%%%%%%%%%%%%_IF STAGE_%%%%%%%%%%%%%%%%%%%%%%

BW_IF=10000; 
%bandpass signal with
% F_stop1= 15000 khz
% F_pass1= 20000 khz
% F_pass2= 30000 khz
% F_stop2= 35000 khz
IF=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',(F_IF-BW_IF),(F_IF-BW_IF/2),(F_IF+BW_IF/2),(F_IF+BW_IF),60,1,60,F1);
IF1=design(IF,'butter');
IF_OUT=filter(IF1,Mixed_Signal);
figure(6)
plot(f1,abs(fftshift(fft(IF_OUT))));    
title('BBC ARABIA after IF_Stage')
xlabel('Frequency in HZ' )
ylabel ('Amplitude')

%%%%%%%%%%%%%%%%%%%%_Base Band Detection_%%%%%%%%%%%%%%%%%%%%%%%
IF_Carrier=cos(2*pi*F_IF*N1*TS1);
IF_OUT1=IF_OUT.*IF_Carrier;
figure(7)
plot(f1,abs(fftshift(fft(IF_OUT1))));
title('BBC ARABIA Base band detection')
xlabel('Frequency in HZ' )
ylabel ('Amplitude')

%%%%%%%%%%%%%%%%%_LOW PASS FILTER_%%%%%%%%%%%%%%%%%%%
LP=fdesign.lowpass('Fp,Fst,Ap,Ast',F_IF,F_IF+1000,0.1,50,F1); %low pass filter
LP1=design(LP,'butter');
LP_OUT=filter(LP1,IF_OUT1);
figure(8)
plot(f1,abs(fftshift(fft(LP_OUT))));
title('BBC ARABIA after Low Pass filter')
xlabel('Frequency in HZ' )
ylabel ('Amplitude')


%%%%%%%%%%%%%%%%%_Down Sample_%%%%%%%%%%%%%%%%
OUTPUT1=downsample(LP_OUT,12);
sound(OUTPUT1,F1/12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



