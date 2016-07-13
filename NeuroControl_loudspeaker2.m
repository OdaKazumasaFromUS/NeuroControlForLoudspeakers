
close all
clear all

%%
%�B��w�̃��C���[��
h = 100;
%�e���C���[�̃��j�b�g��
%h = [size(x,1)+1;h(:)+1;size(xref,1)];
SVDn = 6;
h = [SVDn*2+1;h(:)+1;SVDn];
%���C���[��
L = numel(h);
W1initial = randn(h(1),h(2));
W2initial = randn(h(2),h(3));

%
yOld = 0;

%%
%%�����̃p�����[�^�����낢��ς��Ă��������B%%
P = SVDn;%% ���莟���y��FIR�t�B���^�̃^�b�v��
m = 1;%�󉹓_�̐�(1~8)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%SD(Standard Deviation)�̌v���ݒ�%%
F1=50;  %SD�v���͈͂̉������g��
F2=500;  %SD�v���͈͂̏�����g��

%% �X�s�[�J�[Simulink���f������IR��ǂݍ���
%�X�s�[�J�̃p�����[�^�i��������MFB�C���m�Ȃ̂��E�E�E�j
Md = 13.8e-3;
Rmd = 2.34;
Re0 = 0;%�d���̓d�C�C���s�[�_���X��������
L0 = 0;
Red = 7.4;%�R�C���̓d�C�C���s�[�_���X��������
Ld = 0.97e-3;
Kd = 4.538e+3;
%Kd = 1e+3;
Cd = 21820;
%V = 1.6e-2;
Ss = 0.018;
%Ka = 1.4e+5*Ss^2./V;
K = Kd;
A = 5.62;
OutGain = 10e+5;

%�X�s�[�J�[�̏�ԋ�ԕ�����
AaLinear = [0,1,0,0; -(K)/Md,-Rmd/Md,0,A/Md; 0,0,0,1; 0,-A/(L0+Ld),-(1/Cd)/(L0+Ld),-(Re0+Red)/(L0+Ld)];%LQR�p
Ba = [0,0,0,1/(L0+Ld)]';
Ca = [OutGain,0,0,0];
Da = 0;
sysPlant = ss(AaLinear,Ba,Ca,Da);

%����`���̓���
gammaK = 1000;%�X�e�B�t�l�X����`�W��
%gammaK = 0;
Anlr = zeros(4,4);%�X�e�B�t�l�X�̔���`��
Anlr(2,1) = -K/Md*gammaK;
gammaA = -1000;%�͌W������`�W��
%gammaA = 0;
Anlr3 = zeros(4,4);
Anlr3(2,4) = A/Md*gammaA;
Anlr2 = zeros(4,4);
Anlr2(4,2) = -A/(L0+Ld)*gammaA;


%goalAmp=1e-5;%�ڕW�̐U��
AmpInput = 1;
fInput =60;
IMD = 0;
fIMD = 10;
SignalChange = 0;%0�ŃC���p���X�C1�Ő����g�B�X�e�b�v�͖����C�����ɒǏ]�������Ȃ����C���ł͂���K�v���Ȃ�
LowPassOnOff = 0;%0�ŃI�t�C1�ŃI��
THDno = 5;%THD�ł����̍����g���[�h���Ƃ邩�B1�ňꎟ���[�h�����B2�œ�ڂ̃��[�h�܂�
alpha = 0.98;
k1 = 10;%�؂�ւ��Q�C��
pertub = 1;%1�Őۓ��I��
noise = 0.3;
noise = 0;
SimTime = 2;
fs = 48000;
fsSim = 48000;
feedback = zeros(1,4);
if SignalChange == 0
    PlantInputGain = 10;
    inputGainForPushingNonlinearity = 10;
    inverseFilterGainCompensation = 1/0.0054;
    plantInputGainCompensation = 1/1.15;
    plantInputGainCompensationInv = 1/1.15;
    plantInputGainCompensationFb = 1/1.15;
    inputGainForPushingNonlinearityForSin = 1;
else
    PlantInputGain = 10;
    inputGainForPushingNonlinearity = 5;
    inverseFilterGainCompensation = 1/0.0054;
    plantInputGainCompensation = 1/1.15;
    plantInputGainCompensationInv = 1/1.15;
    plantInputGainCompensationFb = 1/1.15;
    inputGainForPushingNonlinearityForSin = 1;
end

temp = SignalChange;
SignalChange = 0;
pulseStart = 0;
open_system('Loudspeaker_Nonlinear_SVD_OutPress','loadonly')
sim('Loudspeaker_Nonlinear_SVD_OutPress')

SignalChange = temp;
clear temp

%�ψʂ������o��
irSim = Press;
irSim = irSim(1:end-1,1);
%irSim = (irSim./max(irSim))';%�Ƃ肠�������K�����Ă���C(1,:)�`����
irSim = irSim';

%irSim = filter(b,a,irSim);
%�_�E���T���v�����O
downrate = 1;
t = linspace(0,2,2*fs);
fs = fs/downrate;%sampling rate after decreased
%irSimDown = decimate(irSim,downrate,'fir');
tDown = t(1:2000);
irSimDown = irSim(1,1:2000);



%�X�s�[�J�[���U���ԏ�ԋ�ԃ��f��
sysSim = ss(AaLinear,Ba,Ca,Da,1/fs);

figure
plot(irSimDown)

%%
%���͐M���ւ�LPF�Ɏg�p����`�F�r�V�F�t�t�B���^
nfilt = 4;
rip = .05;	% passband ripple in dB
%[b,a] = cheby1(nfilt, rip, 150/fs);
[b,a] = cheby1(nfilt, rip, [70/fs,20000/fs],'bandpass');

%% �F�X�Ȓl�̐ݒ�

len=length(irSimDown(1,:)); %�v�����g��IR��

%fft�|�C���g�̐ݒ�
point = 2;
while point < length(irSimDown(1,:,1));
    point = point*2;
end

%���Ԏ��̐ݒ�
Time = 0:1/fs:(length(irSim(1,:))-1)/fs;
time = (length(irSim(1,:)-1))/fs;
TimeDown = 0:1/fs:(length(irSimDown(1,:))-1)/fs;
timeDown = (length(irSimDown-1))/fs;

%% SVD���f����
%[A,B,C,D] = era(irSimDown,P); %���ْl����@
Atemp = A;
[A,B,C,D] = era(irSimDown,P);
Dtemp = D;
sys_era = ss(A,B,C,D,1/fs);
sys_era_c = d2c(sys_era,'zoh');

[svd_imp,tSvd] = impulse(sys_era);
%era_model = impulse(sys_era,time);

figure
plot(tDown,irSimDown,tSvd,svd_imp)
title('Comparison between Plant and SVDmodel')
legend('Plant','SVDmodel')

%%
i = 1;
for noise = [0]
    %% �����g���͂̂Ƃ��v�����g�̐����g����
    open_system('Loudspeaker_Nonlinear_SVD2_Press','loadonly')
    sim('Loudspeaker_Nonlinear_SVD2_Press')
    %�ψʂ������o��
    SinSim = Velocity;
    SinSim = SinSim(1:end-1,1);
    %SinSim = (SinSim./max(SinSim));%�Ƃ肠�������K�����Ă���C(1,:)�`����
    SinSim = SinSim';

    %�_�E���T���v�����O
    SinSimDown = decimate(SinSim,downrate);

    %% �����Ńt�B�[�h�o�b�N��������i�ɂ̈ʒu���ω�����C������܂߂Č�ŋt�t�B���^��������B���ʂƂ��ăt�B�[�h�o�b�N����݂��Ă��邩��O���}���⃍�o�X�g�������܂��B�j
    Q = sys_era_c.c'*sys_era_c.c;
    %Q=eye(6);
    R=10;
    FbSVD = lqr2(sys_era_c.a,sys_era_c.b,Q,R);
    %FbSVD = zeros(1,6);
    sys_fbC = ss(sys_era_c.a-sys_era_c.b*FbSVD,sys_era_c.b,sys_era_c.c-sys_era_c.d*FbSVD,sys_era_c.d);
    sys_fb = c2d(sys_fbC,1/fs);

    figure
    impulse(sys_fbC,sys_era_c)
    
    %% zpk�ŕ���

    %SVD�ECAPZ���f���Ƃ��ɋɗ냂�f���ɕϊ�����
    zpkerasys=zpk(sys_era);
    %zpkFbsys=zpk(sys_fb);
    %zpkcapzsys=zpk(sys_capz);

    %% �i�[�s��̗p��

    rec=zeros(m,5);%SD�̋L�^
    ti=struct('s_m',{},'s_i',{},'c_m',{},'c_i',{},'fir',{});
    %% FIR�t�B���^�̂��߂̐ݒ�

    H=zeros(m,len); %IR���i�[����s��
    D=zeros(m,len); %�ڕW�M�����i�[����s��
    Ng=P*2; % �t�B���^�̃^�b�v��


    %% �󉹓_�̐������J��Ԃ�

    %fft_sd( ap(ir(1,:)),fs,F1,F2)

    for k=1:1:m
        %% ��_�̈��艻
        %�e���f���̗�_(�t�B���^�ł͋ɂɂȂ�)�ɂ��āC�s����Ȃ��̂����艻������

        %SVD
        zeroBefore = zpkerasys(k,1).z{1,1};
        zeroAfter = zeroBefore;
        z_svd=zeros(length(zeroBefore)+1,1); %�u��������ꂽ��_���i�[����s��B�t�t�B���^���v���p�[�ɂ��邽�߁C��_��������B�j

        %1�ȏ�̗�_�̐�Βl��u�������āA���f���̌`�ɕ�������ߒ�
        %for kk=1:1:P-1%���_�l�d�l
        for kk=1:1:length(zeroBefore)
            if abs(zeroBefore(kk)) >= 1
                zeroAfter(kk)=1/conj(zeroBefore(kk));
            end

            %     if rz_capz(kk) >= 1
            %         rz_capz(kk)=1/rz_capz(kk);
            %     end

            %     if rz_fb(kk) >= 1
            %         rz_fb(kk)=1/rz_fb(kk);
            %     end

            %�ȉ������艻���ꂽ��_
            z_svd(kk,1)=zeroAfter(kk);
            %z_capz(kk,1)=rz_capz(kk)*exp(1i*theta_capz_z(kk));
            %z_svdFb(kk,1)=rz_fb(kk)*exp(1i*theta_fb_z(kk));
        end

        %% �����O���z���ʃv���b�g
        figure
        subplot(211)
        zplane([cell2mat(zpkerasys.z)],[cell2mat(zpkerasys.p)])
        legend('before stabilization')
        subplot(212)
        zplane(z_svd,cell2mat(zpkerasys.p))
        legend('after stabilization')

        %% �ŏ����@�Ɋ�Â�FIR�t�B���^
%         %���Ƃ̔�r�̂��߂ɍ\������
% 
%         H(k,:)=irSimDown(k,:);
% 
%         [waste, dd]=max(abs(irSimDown(k,:)));
%         delay1=Ng/2+dd;
%         d1=zeros(len,1);
%         d1(delay1+1,1)=1;
% 
%         D(k,:)=d1';
% 
%         Ehd=ccor_vec(H,D,Ng);
%         Ehh=ccor_mat(H,H,Ng);
% 
%         c1=linsolve(Ehh,Ehd);
%         y_fir=conv(c1,irSimDown(k,:));
%         y_fir=ap(y_fir);

        %% �t�B���^�̎���MPAP_IIR(SVD)
        zero = zpkerasys(1,1).p{1,1};
        pole = z_svd;
        gain = 1/zpkerasys(1,1).k;
        sys_svd_iird = zpk(zero,pole,gain,1/fs);
        %sys_svd_iirc = d2c(sys_svd_iird);
        InputGain = 1;
        open_system('SVD_IIR_Sim_Cont_Press','loadonly')
        sim('SVD_IIR_Sim_Cont_Press')
        PlantY = Velocity;
        %PlantY = PlantY./max(abs(PlantY));
        PlantY = PlantY(1:end-1,1)';
        svd_iir = PlantY;

        Fb = FbSVD;
        %Fb = zeros(1,6);

        %�I�u�U�[�o�Q�C������
        Qobs = 1000*eye(6);
        %Qobs(2,2) = 0;
        %Qobs(1,1) = 0;
        %Qobs(3,3) = 0;
        %Qobs(4,4) = 0;
        ObserverNorm = 1;
        Robs = 1e-6;
        Gt = lqr2(sys_era_c.a',sys_era_c.c',Qobs,Robs);
        G = Gt';
        %G = zeros(6,1);
        %feedback= zeros(1,6);
        open_system('NeuroControl_Loudspeaker_offline2','loadonly')
        sim('NeuroControl_Loudspeaker_offline2')

        figure
        plot(time,y,time,yModel,time,ybar)
        legend('y','yModel','ybar')
        xlim([0,0.05])

        PlantY = Velocity;
        %PlantY = PlantY./max(abs(PlantY));
        PlantY = PlantY';
        PlantY = PlantY(1,1:end-1);
        svd_iir_fb = PlantY;

        %close all
        %%
        figure
        plot(time(1:end-1,1),irSim,time(1:end-1,1),svd_iir_fb)
        figure
        plot(time,ybar./max(ybar),time,yModel./max(yModel),time(1:end-1,1),SinSim./max(SinSim))

       
    end
    IrPlant(i,:) = SinSim;
    IrSvd(i,:) = svd_iir;
    IrSvdFb(i,:) = svd_iir_fb;
    i = i + 1;
end