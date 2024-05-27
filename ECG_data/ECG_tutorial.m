%% ȭ�� �� ���� �ʱ�ȭ
clear all
close all
clc
%% �⺻ ���� ���� �� �ҷ�����

x1 = load('ECG3.dat'); % ECG ��ȣ �ҷ�����
fs = 200;              % Sampling rate
N = length (x1);       % ��ȣ�� �� ����
t = [0:N-1]/fs;        % �ð� �� ���
%% �׸� �׸���
figure(1)
subplot(2,1,1)
plot(t,x1)
xlabel('second');ylabel('Volts');title('RAW ECG ��ȣ')
subplot(2,1,2)
plot(t(200:600),x1(200:600))
xlabel('second');ylabel('Volts');title('1��~3�� ���� ��ȣ')
xlim([1 3])

%% DC component ���� �� ����ȭ

x1 = x1 - mean(x1);    % DC ���۳�Ʈ ����
x1 = x1/ max(abs(x1)); % �ִ밪�� 1�� ����ȭ

figure(2)
subplot(2,1,1)
plot(t,x1)
xlabel('second');ylabel('Volts');title('DC���۳�Ʈ ���� �� ����ȭ�� ECG��ȣ')
subplot(2,1,2)
plot(t(200:600),x1(200:600))
xlabel('second');ylabel('Volts');title('1~3�� ���� ��ȣ')
xlim([1 3])

%% �����������(LPF)
load('lpf.mat')

h_LP=filter(lpf.b,lpf.a,[1 zeros(1,12)]); % ������ impulse response

x2 = conv(x1 ,h_LP); % �������
x2 = x2 (6+[1:N]); % ������ ����
x2 = x2/ max(abs(x2)); %����ȭ

figure(3)
subplot(2,1,1)
plot(t,x2)
xlabel('second');ylabel('Volts');title('LPF�� ECG')
xlim([0 max(t)])
subplot(2,1,2)
plot(t(200:600),x2(200:600))
xlabel('second');ylabel('Volts');title('1~3�� ���� ��ȣ')
xlim([1 3])


%% �����������(HPF)
load('hpf.mat')

h_HP=filter(hpf.b,hpf.a,[1 zeros(1,32)]); % ������ impulse response

x3 = conv (x2 ,h_HP); % �������
x3 = x3 (16+[1:N]); % ������ ����
x3 = x3/ max(abs(x3)); % ����ȭ

figure(4)
subplot(2,1,1)
plot(t,x3)
xlabel('second');ylabel('Volts');title('HPF�� ECG')
xlim([0 max(t)])
subplot(2,1,2)
plot(t(200:600),x3(200:600))
xlabel('second');ylabel('Volts');title('1~3�� ���� ��ȣ')
xlim([1 3])


%% Derivative Filter
load('df.mat')

x4 = conv (x3 ,df.h);
x4 = x4 (2+[1:N]);
x4 = x4/ max(abs(x4));

figure(5)
subplot(2,1,1)
plot(t,x4)
xlabel('second');ylabel('Volts');title('DF�� ECG ��ȣ')
subplot(2,1,2)
plot(t(200:600),x4(200:600))
xlabel('second');ylabel('Volts');title('1~3�� ���� ��ȣ')
xlim([1 3])

%% Squaring

x5 = x4 .^2; % ��ȣ�� ����
x5 = x5/ max(abs(x5)); % ����ȭ

figure(6)
subplot(2,1,1)
plot(t,x5)
xlabel('second');ylabel('Volts');title('Squaring ���� ��ȣ')
subplot(2,1,2)
plot(t(200:600),x5(200:600))
xlabel('second');ylabel('Volts');title('1~3�� ���� ��ȣ')
xlim([1 3])


%% Moving Window Integration

load('mw.mat')

x6 = conv (x5 ,mw.h);
x6 = x6 (15+[1: N]);
x6 = x6/ max( abs(x6 ));

figure(7)
subplot(2,1,1)
plot([0:length(x6)-1]/fs,x6)
xlabel('second');ylabel('Volts');title('�̵���� ���� ��ȣ')
subplot(2,1,2)
plot(t(200:600),x6(200:600))
xlabel('second');ylabel('Volts');title('1~3�� ���� ��ȣ')
xlim([1 3])

%% Thresholding

max_h = max(x6);
thresh = mean (x6);
poss_reg =(x6>thresh*max_h)'; % threshold for QRS

figure (8)
subplot(2,1,1)
plot (t(200:600),x6(200:600)) % Integrated signal
hold on; line(t, 0.2, 'LineStyle', '-', 'Color', 'red'); % Threshold
xlabel('second');ylabel('Integrated')
xlim([1 3])

subplot(2,1,2)
plot(t, x1) % basline corrected signal
hold on;
plot (t, poss_reg, 'r') % ������ QRS ����
xlabel('second');ylabel('Integrated')
xlim([1 3])

%% Find QRS peaks
left = find(diff([0 poss_reg])==1);
right = find(diff([poss_reg 0])==-1);

for i=1:length(left)
    [R_value(i) R_loc(i)] = max( x1(left(i):right(i)) );
    R_loc(i) = R_loc(i)-1+left(i); % add offset

    [Q_value(i) Q_loc(i)] = min( x1(left(i):R_loc(i)) );
    Q_loc(i) = Q_loc(i)-1+left(i); % add offset

    [S_value(i) S_loc(i)] = min( x1(left(i):right(i)) );
    S_loc(i) = S_loc(i)-1+left(i); % add offset
end

figure(9)
subplot(2,1,1)
title('ECG Signal with R points');
plot (t,x1, t(R_loc) ,R_value , 'r^', t(S_loc) ,S_value, '*',t(Q_loc) , Q_value, 'o');
legend('ECG','R','S','Q');
subplot(2,1,2)
plot (t,x1, t(R_loc) ,R_value , 'r^', t(S_loc) ,S_value, '*',t(Q_loc) , Q_value, 'o');
xlim([1 3])

%% RR interval and bpm
mRR = mean(diff(R_loc))/fs % mean RR (s);

bpm = 60./mRR % 1min / meanRR = bpm






