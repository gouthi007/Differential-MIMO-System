clc;
clear all;
close all;

%Given data
M=8;
N=34;
%Carrier frequency in Hz
fc = 800e6;
% number of samples
ns=1e5; 
%Speed of vehicle in m/s
l = 22.352;
%BandWidth in Hz
BW = 30e3;
%Speed of light in m/s
c = 3e8; 
%Wavelength
n1 = c/fc; 
%Maximum doppler frequency in Hz
maxfd = l/n1; 
Ts = 1/BW;

%Part A
%Parameter Initialization
m1=0;
m2=0;
m3=0;
m4=0;
final1=0;
final2=0;
final3=0;
final4=0;
sum1=0;
sum2=0;
sum3=0;
sum4=0;
 
%Using the Modified Jake's Simulator to generate the Channel Coefficients
%Hadamard matrix
h=hadamard(8);
%Number of Channel samples
for k=1:ns 
      m1=0;m2=0;m3=0;m4=0;
%Number of Oscillators
  for n=1:M 
    on=(2*pi*n)/N;
    bn=(pi*n)/(M+1);
    m1 = m1 + (h(1,n)*exp(1i*bn)*cos(2*pi*k*maxfd*Ts*cos(on)+(2*2)*bn));% for rayleigh channel 1 
    m2 = m2 + (h(2,n)*exp(1i*bn)*cos(2*pi*k*maxfd*Ts*cos(on)+(2*3)*bn));% for rayleigh channel 2
    m3 = m3 + (h(3,n)*exp(1i*bn)*cos(2*pi*k*maxfd*Ts*cos(on)+(2*4)*bn));% for rayleigh channel 3
    m4 = m4 + (h(4,n)*exp(1i*bn)*cos(2*pi*k*maxfd*Ts*cos(on)+(2*5)*bn));% for rayleigh channel 4
  end
  
  h1(k)= abs(m1);
  h2(k)= abs(m2);
  h3(k)= abs(m3);
  h4(k)= abs(m4);
  %Power of channel 1
  sum1=sum1+(h1(k).*h1(k));  
  %Power of channel 2
  sum2=sum2+(h2(k).*h2(k));   
  %Power of channel 3
  sum3=sum3+(h3(k).*h3(k));  
  %Power of channel 4
  sum4=sum4+(h4(k).*h4(k));   
end
  
%Average Power of channel 1
final1 = sum1/ns; 
%Average Power of channel 2
final2 = sum2/ns;
%Average Power of channel 3
final3 = sum3/ns; 
%Average Power of channel 4
final4 = sum4/ns; 

%Normlized Channel Coefficients of channel 1
hp1= h1/sqrt(final1); 
%Normlized Channel Coefficients of channel 2
hp2= h2/sqrt(final2); 
%Normlized Channel Coefficients of channel 3
hp3= h3/sqrt(final3); 
%Normlized Channel Coefficients of channel 4
hp4= h4/sqrt(final4); 
 
% plotting the channels 
t=1:500;
figure(1)
% rayleigh channel 1
plot(t, 10*log(hp1(t)),'g');
hold on;
% rayleigh channel 2
plot(t, 10*log(hp2(t)),'r');
hold on;
% rayleigh channel 3
plot(t, 10*log(hp3(t)),'blue');
hold on;
% rayleigh channel 4
plot(t, 10*log(hp4(t)),'black');
title('Initial 500 Channel Samples for the four Rayleigh Channels');
xlabel('Channel Samples');
ylabel('|H(t)| in dB'); 
legend('Channel 1','Channel 2','Channel 3','Channel 4');

%Part B
SNR = 1:1:15;
Mt = 2;
Mr = 1;
iterations=1e5;
for l=1:4
codeword(:,:,l) = sqrt(Mt)*[exp(1i*(l-1)*pi/2) 0;0 exp(1i*(l-1)*pi/2)];% codeword generation
end

s(:,:,1) = sqrt(Mt)*[1 0 ;0 1];
h1t(:,:,1) = [hp1(1) ; hp2(1)];
N1 = sqrt(1/2)*randn(Mt,Mr)+1i*sqrt(1/2)*randn(Mt,Mr);
 
for i= 1:length(SNR)
y1t(:,:,1) = sqrt(SNR(i)/Mt)*s(:,:,1)*h1t(:,:,1)+N1;
for l = 2:iterations
r1 = randi(4);
C(:,:,l) = codeword(:,:,r1);
s(:,:,l) = (1/sqrt(Mt))*(C(:,:,l)*s(:,:,l-1));
h1t(:,:,l) = [hp1(l);hp2(l)];
%Transmitted signal
y1t(:,:,l) = sqrt(SNR(i)/Mt)*(s(:,:,l)*h1t(:,:,l))+N1;
end
err1=0;
for v = 2:iterations
n=zeros(1,4);

%Decoding 
for u=1:4
n(u)=norm(y1t(:,:,v)-((1/sqrt(Mt))*codeword(:,:,u)*y1t(:,:,v-1)),'fro')^2;
end
[q,r]=min(n);
ct(:,:,v)=codeword(:,:,r);
tf=isequal(ct(:,:,v),C(:,:,v));
% if(tf==0)
if((ct(1,1,v)~=C(1,1,v))||(ct(2,2,v)~=C(2,2,v))||(ct(1,2,v)~=C(1,2,v))||(ct(2,1,v)~=C(2,1,v)))
err1=err1+1;
end
end
% SER, when Mr=1
ser1(i)=err1/(iterations-1); 
end 
 
figure(2)
semilogy(SNR,ser1,'r');
title('SER vs SNR for Mt=2 and Mr=1')
xlabel('SNR in db')
ylabel('SER');
grid on;
hold on;
legend('1 Receiver');

%Part C
Mt = 2;
Mr = 2;
%Codeword generation
for l=1:4
codeword(:,:,l) = sqrt(Mt)*[exp(1i*(l-1)*pi/2) 0;0 exp(1i*(l-1)*pi/2)];
end
s(:,:,1) = sqrt(Mt)*[1 0 ;0 1];
h2t(:,:,1) = [hp1(1) hp2(1);hp3(1) hp4(1)];
N2 = sqrt(1/2)*randn(Mt,Mr)+1i*sqrt(1/2)*randn(Mt,Mr);
for i= 1:length(SNR)
y2t(:,:,1) = sqrt(SNR(i)/Mt)*s(:,:,1)*h2t(:,:,1)+N2;
for l = 2:iterations
r1 = randi(4);
C(:,:,l) = codeword(:,:,r1);
s(:,:,l) = (1/sqrt(Mt))*(C(:,:,l)*s(:,:,l-1));
h2t(:,:,l) = [hp1(l) hp2(l);hp3(1) hp4(l) ];                   
%Transmitted signal 
y2t(:,:,l) = sqrt(SNR(i)/Mt)*(s(:,:,l)*h2t(:,:,l))+N2;
end
err2=0;
for v = 2:iterations
n=zeros(1,4);
%Decoding
for u=1:4
n(u)=norm(y2t(:,:,v)-((1/sqrt(Mt))*codeword(:,:,u)*y2t(:,:,v-1)),'fro')^2;
end
[q,r]=min(n);
ct(:,:,v)=codeword(:,:,r);
tf=isequal(ct(:,:,v),C(:,:,v));
% if(tf==0)
if((ct(1,1,v)~=C(1,1,v))||(ct(2,2,v)~=C(2,2,v))||(ct(1,2,v)~=C(1,2,v))||(ct(2,1,v)~=C(2,1,v)))
err2=err2+1;
end
end
%SER, when Mr=2
ser2(i)=err2/(iterations-1);
end

figure(3)
semilogy(SNR,ser2);
title('SER vs SNR for Mt=2 and Mr=2')
xlabel('SNR in db')
ylabel('SER');
grid on;
hold on;      
legend('2 Receivers');

figure(4)
semilogy(SNR,ser1,'r',SNR,ser2,'g');
legend('2 Receivers','1 Receiver');
title('SER vs SNR for Mt=2 and Mr=1,2')
xlabel('SNR in db')
ylabel('SER');
grid on;
hold on;

