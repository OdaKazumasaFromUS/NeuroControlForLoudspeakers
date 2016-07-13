
close all
clear all

%%
%隠れ層のレイヤー数
h = 100;
%各レイヤーのユニット数
%h = [size(x,1)+1;h(:)+1;size(xref,1)];
SVDn = 6;
h = [SVDn*2+1;h(:)+1;SVDn];
%レイヤー数
L = numel(h);
W1initial = randn(h(1),h(2));
W2initial = randn(h(2),h(3));

%
yOld = 0;

%%
%%ここのパラメータをいろいろ変えてください。%%
P = SVDn;%% 同定次数及びFIRフィルタのタップ数
m = 1;%受音点の数(1~8)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%SD(Standard Deviation)の計測設定%%
F1=50;  %SD計測範囲の下限周波数
F2=500;  %SD計測範囲の上限周波数

%% スピーカーSimulinkモデルからIRを読み込み
%スピーカのパラメータ（研究室のMFB，正確なのか・・・）
Md = 13.8e-3;
Rmd = 2.34;
Re0 = 0;%電源の電気インピーダンス直流成分
L0 = 0;
Red = 7.4;%コイルの電気インピーダンス直流成分
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

%スピーカーの状態空間方程式
AaLinear = [0,1,0,0; -(K)/Md,-Rmd/Md,0,A/Md; 0,0,0,1; 0,-A/(L0+Ld),-(1/Cd)/(L0+Ld),-(Re0+Red)/(L0+Ld)];%LQR用
Ba = [0,0,0,1/(L0+Ld)]';
Ca = [OutGain,0,0,0];
Da = 0;
sysPlant = ss(AaLinear,Ba,Ca,Da);

%非線形性の導入
gammaK = 1000;%スティフネス非線形係数
%gammaK = 0;
Anlr = zeros(4,4);%スティフネスの非線形性
Anlr(2,1) = -K/Md*gammaK;
gammaA = -1000;%力係数非線形係数
%gammaA = 0;
Anlr3 = zeros(4,4);
Anlr3(2,4) = A/Md*gammaA;
Anlr2 = zeros(4,4);
Anlr2(4,2) = -A/(L0+Ld)*gammaA;


%goalAmp=1e-5;%目標の振幅
AmpInput = 1;
fInput =60;
IMD = 0;
fIMD = 10;
SignalChange = 0;%0でインパルス，1で正弦波。ステップは無理，直流に追従させられないし，音ではする必要もない
LowPassOnOff = 0;%0でオフ，1でオン
THDno = 5;%THDでいくつの高調波モードをとるか。1で一次モードだけ。2で二つ目のモードまで
alpha = 0.98;
k1 = 10;%切り替えゲイン
pertub = 1;%1で摂動オン
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

%変位だけ取り出す
irSim = Press;
irSim = irSim(1:end-1,1);
%irSim = (irSim./max(irSim))';%とりあえず正規化してから，(1,:)形式に
irSim = irSim';

%irSim = filter(b,a,irSim);
%ダウンサンプリング
downrate = 1;
t = linspace(0,2,2*fs);
fs = fs/downrate;%sampling rate after decreased
%irSimDown = decimate(irSim,downrate,'fir');
tDown = t(1:2000);
irSimDown = irSim(1,1:2000);



%スピーカー離散時間状態空間モデル
sysSim = ss(AaLinear,Ba,Ca,Da,1/fs);

figure
plot(irSimDown)

%%
%入力信号へのLPFに使用するチェビシェフフィルタ
nfilt = 4;
rip = .05;	% passband ripple in dB
%[b,a] = cheby1(nfilt, rip, 150/fs);
[b,a] = cheby1(nfilt, rip, [70/fs,20000/fs],'bandpass');

%% 色々な値の設定

len=length(irSimDown(1,:)); %プラントのIR長

%fftポイントの設定
point = 2;
while point < length(irSimDown(1,:,1));
    point = point*2;
end

%時間軸の設定
Time = 0:1/fs:(length(irSim(1,:))-1)/fs;
time = (length(irSim(1,:)-1))/fs;
TimeDown = 0:1/fs:(length(irSimDown(1,:))-1)/fs;
timeDown = (length(irSimDown-1))/fs;

%% SVDモデル化
%[A,B,C,D] = era(irSimDown,P); %特異値分解法
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
    %% 正弦波入力のときプラントの正弦波応答
    open_system('Loudspeaker_Nonlinear_SVD2_Press','loadonly')
    sim('Loudspeaker_Nonlinear_SVD2_Press')
    %変位だけ取り出す
    SinSim = Velocity;
    SinSim = SinSim(1:end-1,1);
    %SinSim = (SinSim./max(SinSim));%とりあえず正規化してから，(1,:)形式に
    SinSim = SinSim';

    %ダウンサンプリング
    SinSimDown = decimate(SinSim,downrate);

    %% ここでフィードバック導入する（極の位置が変化する，それを含めて後で逆フィルタ生成する。結果としてフィードバックを内在しているから外乱抑圧やロバスト性が生まれる。）
    Q = sys_era_c.c'*sys_era_c.c;
    %Q=eye(6);
    R=10;
    FbSVD = lqr2(sys_era_c.a,sys_era_c.b,Q,R);
    %FbSVD = zeros(1,6);
    sys_fbC = ss(sys_era_c.a-sys_era_c.b*FbSVD,sys_era_c.b,sys_era_c.c-sys_era_c.d*FbSVD,sys_era_c.d);
    sys_fb = c2d(sys_fbC,1/fs);

    figure
    impulse(sys_fbC,sys_era_c)
    
    %% zpkで分解

    %SVD・CAPZモデルともに極零モデルに変換する
    zpkerasys=zpk(sys_era);
    %zpkFbsys=zpk(sys_fb);
    %zpkcapzsys=zpk(sys_capz);

    %% 格納行列の用意

    rec=zeros(m,5);%SDの記録
    ti=struct('s_m',{},'s_i',{},'c_m',{},'c_i',{},'fir',{});
    %% FIRフィルタのための設定

    H=zeros(m,len); %IRを格納する行列
    D=zeros(m,len); %目標信号を格納する行列
    Ng=P*2; % フィルタのタップ数


    %% 受音点の数だけ繰り返し

    %fft_sd( ap(ir(1,:)),fs,F1,F2)

    for k=1:1:m
        %% 零点の安定化
        %各モデルの零点(フィルタでは極になる)について，不安定なものを安定化させる

        %SVD
        zeroBefore = zpkerasys(k,1).z{1,1};
        zeroAfter = zeroBefore;
        z_svd=zeros(length(zeroBefore)+1,1); %置き換えられた零点を格納する行列。逆フィルタをプロパーにするため，零点一個足した。）

        %1以上の零点の絶対値を置き換えて、複素数の形に復元する過程
        %for kk=1:1:P-1%理論値仕様
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

            %以下が安定化された零点
            z_svd(kk,1)=zeroAfter(kk);
            %z_capz(kk,1)=rz_capz(kk)*exp(1i*theta_capz_z(kk));
            %z_svdFb(kk,1)=rz_fb(kk)*exp(1i*theta_fb_z(kk));
        end

        %% 処理前後のz平面プロット
        figure
        subplot(211)
        zplane([cell2mat(zpkerasys.z)],[cell2mat(zpkerasys.p)])
        legend('before stabilization')
        subplot(212)
        zplane(z_svd,cell2mat(zpkerasys.p))
        legend('after stabilization')

        %% 最小二乗法に基づくFIRフィルタ
%         %他との比較のために構成する
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

        %% フィルタの実装MPAP_IIR(SVD)
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

        %オブザーバゲイン決定
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