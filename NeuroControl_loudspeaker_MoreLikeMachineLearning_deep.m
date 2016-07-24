
%close all

%% SVDやスピーカーの?ﾝ定
%%ここのパラ?ータをいろいろ変えてください?B%%
P = 6;%% 同定次?萩yびFIRフィルタのタップ??
m = 1;%受音点の??1~8)

%% スピーカーSimulinkモデルからIRを読み?桙ﾝ
%スピーカのパラ?ータ?i研究室のMFB?C?ｳ確なのか?E?E?E?j
Md = 13.8e-3;
Rmd = 2.34;
Re0 = 0;%電源の電気インピーダンス直流?ｬ分
L0 = 0;
Red = 7.4;%コイルの電気インピーダンス直流?ｬ分
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

%スピーカーの?ﾔ空間方程式
AaLinear = [0,1,0,0; -(K)/Md,-Rmd/Md,0,A/Md; 0,0,0,1; 0,-A/(L0+Ld),-(1/Cd)/(L0+Ld),-(Re0+Red)/(L0+Ld)];%LQR用
Ba = [0,0,0,1/(L0+Ld)]';
Ca = [OutGain,0,0,0];
Da = 0;
sysPlant = ss(AaLinear,Ba,Ca,Da);

%非??`?ｫの導入
gammaK = 1000;%スティフネス非??`係??
%gammaK = 0;
Anlr = zeros(4,4);%スティフネスの非??`?ｫ
Anlr(2,1) = -K/Md*gammaK;
gammaA = -1000;%力係?粕??`係??
%gammaA = 0;
Anlr3 = zeros(4,4);
Anlr3(2,4) = A/Md*gammaA;
Anlr2 = zeros(4,4);
Anlr2(4,2) = -A/(L0+Ld)*gammaA;


%goalAmp=1e-5;%目標の?U?
AmpInput = 1;
fInput =60;
IMD = 0;
fIMD = 10;
SignalChange = 0;%0でインパルス?C1で?ｳ弦波?Bステップは無??C直流に追?]させられないし?C音ではする必要もない
LowPassOnOff = 0;%0でオフ?C1でオン
THDno = 5;%THDでいくつの?ｒｲ波モードをとるか?B1で一次モードだけ?B2で二つ目のモードまで
alpha = 0.98;
k1 = 10;%?ﾘり替えゲイン
pertub = 1;%1で?ﾛ動オン
noise = 0;
SimTime = 0.1;
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
    inputGainForPushingNonlinearity = 10;
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

%変位だけ取り?oす
irSim = Press;
irSim = irSim(1:end-1,1);
%irSim = (irSim./max(irSim))';%とりあえず?ｳ規化してから?C(1,:)形式に
irSim = irSim';

%irSim = filter(b,a,irSim);
%ダウンサンプリング
downrate = 1;
t = linspace(0,2,2*fs);
fs = fs/downrate;%sampling rate after decreased
%irSimDown = decimate(irSim,downrate,'fir');
tDown = t(1:2000);
irSimDown = irSim(1,1:2000);

%スピーカー離散時間?ﾔ空間モデル
sysSim = ss(AaLinear,Ba,Ca,Da,1/fs);

%%
%入力?M??ﾖのLPFに使用するチェビシェフフィルタ
nfilt = 4;
rip = .05;	% passband ripple in dB
%[b,a] = cheby1(nfilt, rip, 150/fs);
[b,a] = cheby1(nfilt, rip, [70/fs,20000/fs],'bandpass');

%% ?F?Xな値の?ﾝ定

len=length(irSimDown(1,:)); %プラントのIR長

%fftポイントの?ﾝ定
point = 2;
while point < length(irSimDown(1,:,1));
    point = point*2;
end

%時間軸の?ﾝ定
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

%%
i = 1;
k = 1;

%% ?ｳ弦波入力のときプラントの?ｳ弦波応答
open_system('Loudspeaker_Nonlinear_SVD2_Press','loadonly')
sim('Loudspeaker_Nonlinear_SVD2_Press')
%変位だけ取り?oす
SinSim = Velocity;
SinSim = SinSim(1:end-1,1);
%SinSim = (SinSim./max(SinSim));%とりあえず?ｳ規化してから?C(1,:)形式に
SinSim = SinSim';

%ダウンサンプリング
SinSimDown = decimate(SinSim,downrate);

%% ここでフィードバック導入する?i極の位置が変化する?Cそれを含めて後で逆フィルタ?ｶ?ｬする?B結果としてフィードバックを内?ﾝしているから外?抑圧や?バスト?ｫが?ｶまれる?B?j
Q = sys_era_c.c'*sys_era_c.c;
%Q=eye(6);
R=10;
FbSVD = lqr2(sys_era_c.a,sys_era_c.b,Q,R);
%FbSVD = zeros(1,6);
sys_fbC = ss(sys_era_c.a-sys_era_c.b*FbSVD,sys_era_c.b,sys_era_c.c-sys_era_c.d*FbSVD,sys_era_c.d);
sys_fb = c2d(sys_fbC,1/fs);

%% zpkで分解

%SVD?ECAPZモデルともに極零モデルに変換する
zpkerasys=zpk(sys_era);
%zpkFbsys=zpk(sys_fb);
%zpkcapzsys=zpk(sys_capz);

%% 格納?s列の用意

rec=zeros(m,5);%SDの記録
ti=struct('s_m',{},'s_i',{},'c_m',{},'c_i',{},'fir',{});
%% FIRフィルタのための?ﾝ定

H=zeros(m,len); %IRを格納する?s列
D=zeros(m,len); %目標?M??i納する?s列
Ng=P*2; % フィルタのタップ??

%% 零点の安定化
%各モデルの零点(フィルタでは極になる)について?C不安定なものを安定化させる

%SVD
zeroBefore = zpkerasys(k,1).z{1,1};
zeroAfter = zeroBefore;
z_svd=zeros(length(zeroBefore)+1,1); %置き換えられた零点を格納する?s列?B逆フィルタをプ?パーにするため?C零点一個足した?B?j

%1以?繧ﾌ零点の?竭ﾎ値を置き換えて?A複素?狽ﾌ形に復元する過程
%for kk=1:1:P-1%?論値仕様
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

%% ???前後のz平面プ?ット
% figure
% subplot(211)
% zplane([cell2mat(zpkerasys.z)],[cell2mat(zpkerasys.p)])
% legend('before stabilization')
% subplot(212)
% zplane(z_svd,cell2mat(zpkerasys.p))
% legend('after stabilization')

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

%Fb = zeros(1,6);

%オブザーバゲイン決定
Qobs = 1000*eye(6);
%Qobs = sys_era_c.c'*sys_era_c.c;
ObserverNorm = 1;
Robs = 1e-6;
%Robs = 1e-8;
Gt = lqr2(sys_era_c.a',sys_era_c.c',Qobs,Robs);
G = Gt';
%G = zeros(6,1);
%feedback= zeros(1,6);



%% Neural Network
%----------------------
%% データセットの用意
%set the time history
timehist = linspace(0,SimTime,fs*SimTime);
%Numbear of data sets
NumDataSet = 10;
InputDataSet = zeros(NumDataSet,fs*SimTime);
%pink = 0,sine = 1
pinkOrsine = 0;
for i = 1:NumDataSet
    if pinkOrsine == 0
        InputDataSet(i,:) = pinknoise(fs,SimTime);
    else
        f = 60;
        disTime = 0:1:fs*SimTime - 1;
        InputDataSet(i,:) = sin(2*pi*f/fs*disTime);
    end
end
%use the same dataset as the previous one
%InputDataSet = InputTemp;
%clear InputTemp;

%% set properties of NN
%number of neurons in input layer
NumInputNeuro = length(InputDataSet(1,:));
%number of neurons in hidden layer
NumHiddenNeuro = 7000;
%number of neurons in output layer (equals state vector length)
NumOutputNeuro = 6;

h = [NumInputNeuro+1;NumHiddenNeuro+1;NumHiddenNeuro+1;NumOutputNeuro];

%現在のところ推奨の初期化
low = -sqrt(6/(h(1)+h(2)));
high = -low;
W1 = low + (high - low).*rand(h(1),h(2)-1);
W1(1,:) = 0;%initialize by 0 for biases


low = -sqrt(6/(h(2)+h(3)));
high = -low;
W2 = low + (high - low).*rand(h(2),h(3)-1);
W2(1,:) = 0;


low = -sqrt(6/(h(3)+h(4)));
high = -low;
W3 = low + (high - low).*rand(h(3),h(4));
W3(1,:) = 0;

%初期条件を保存したいとき
%--------------------------------
%W1ini = W1;
%W2ini = W2;
%W3ini = W3;

% %特定の初期条件で計算したいとき
% %--------------------------------
% W1 = W1ini;
% W2 = W2ini;
% W3 = W3ini;

%もともとの初期化
%W1 = randn(h(1),h(2)-1);
%W2 = randn(h(2),h(3));

%
yOld = 0;

%% NNをまわして,フィードバックゲインfを得る
iterCount = 1;
eOld = 0;
iterMax = 35;
tune = 1;
tuneOut = 3000;
Fb = zeros(6,1);
jacobi = zeros(P,fs*SimTime+1);
figure
h = animatedline;
movegui(gca,'south')
%警告を非表示
warning('off','all')
clear FbHis
for DataSetNumNow = 1:1:NumDataSet
    iter = 1;
    mse = 1;
    mseLqr = 0;
    eta = 0.000001;
    goNext = 0;
    while iter < iterMax+1
    %while iter < iterMax
       %% forward propagation
        
        %#codegen
        %配列は?cで入れる
        %ニューラルネットの?∮`播で?o力?める?B
        %?o力はフィードバックゲインf
        
        %グ?ーバル変?白闍`
        
        %x,xrefは縦ベクトルでなければならない
        %Outputs of input layer
        Y1 = [1;InputDataSet(DataSetNumNow,:)'];
        
        %Outputs of hidden layer
        %Define activation function at hidden layer here.
        %Y2 = [1;tanh(W1'*Y1)/tune];
        %Y2 = [1;sigmoid(W1'*Y1)/tune];
        %Y2 = [1;softplus(W1'*Y1)/tune];
        Y2 = [1;ReLU(W1'*Y1)/tune];
        
        Y3 = [1;ReLU(W2'*Y2)/tune];
        
        %Outputs of output layer
        %linear
        Y4 = W3'*Y3/tuneOut;
        %tanh
        %Y3 = tanh(W2'*Y2)/100;
        
        %ユニット?o力?i?｡はフィードバックゲイン?jの?ｷ分
        deltaF = Y4 - Fb;
        Fb = Y4;
        
        %得られたFbを適用したときの出力音圧yをシミュレーション。リファレンスyrefも取得
        %フィードバック制御simulinkモデルを用いる
        open_system('SVD_FeedbackForNN','loadonly')
        sim('SVD_FeedbackForNN')
        
        if iter == 1
            yFirst = y;
            Fbtemp = Fb;
            %Fb = FbSVD';
            Fb = 0;
            ytemp = y;
            open_system('SVD_FeedbackForNN','loadonly')
            sim('SVD_FeedbackForNN')
            yLqr = y;
            E = y - yref;
            Ebar = ybar - yref;
            mseLqr = mean(dot(E(:),E(:)))
            mseLqrBar = mean(dot(Ebar(:),Ebar(:)))
            Fb = Fbtemp;
            FbEnerLqr = sum(FbIn.^2);
            eOld = E;
            y = ytemp;
            clear Fbtemp;
            clear ytemp;
        end
        
        %% back propagation
        %learning rate decaying
%         if mod(iter,10) == 0
%             eta = eta/2;
%         end
        
        
        %コスト関数
        
        %ヤコビアンを求める
        deltaY = y-yref - eOld;
        sign = (1./deltaF)*deltaY';
        %かなり雑な近似
%         %signが負のjacobi要素は-1に他は+1
%         jacobi(find(sign>0)) = 1;
%         jacobi(find(sign==0)) = 0;
%         jacobi(find(sign<0)) = -1;
        %jacobi = jacobi*(-1);
        jacobi = sign/max(abs(sign(:)));
        
        %出力層ユニットに対するsigmaK
        %tanh
        %sigK = errorSum*jacobi.*(1-Y3.^2);
        % %linear
        sigK = jacobi*(y - yref)/tuneOut;
       
        %中間層
        sigJ2 = (Y3>0).*(W3*sigK)/tune;
        %tanh
        %sigJ = (1-Y2.^2).*(W2*sigK)/tune;
        %sigmoid
        %sigJ = Y2.*(1-Y2).*(W2*sigK)/tune;
        %softplus
        %sigJ = (1./(1+exp(-Y2))).*(W2*sigK)/tune;
        %ReLU
        sigJ = (Y2>0).*(W2*sigJ2(2:end,:))/tune;
        
        
        
        W1 = W1-eta*(Y1*sigJ(2:end,1)');
        W2 = W2-eta*(Y2*sigJ2(2:end,1)');
        W3dif = Y3*sigK';
        W3dif(:,1) = W3dif(:,1)*eta;
        W3dif(:,2) = W3dif(:,2)*eta;
        W3dif(:,3:6) = W3dif(:,3:6)*eta;
        W3 = W3 - W3dif;
        %W3 = W3-eta*(Y3*sigK');
        
        %eOld = sum(y-yref);
        eOld = y-yref;
        
        E = y - yref;
        Ebar = ybar - yref;
        mse = mean(dot(E(:),E(:)));
        mseBar = mean(dot(Ebar(:),Ebar(:)));
        
%         if mse < 7e+8
%             eta = eta/2;
%         end
        
        %draw error value
        addpoints(h,iterCount,mse);
        drawnow limitrate
        
        %FeedbackGainの時系列保存
        FbHis(iterCount,:) = Fb;
        FbEner(iterCount,:) = sum(FbIn.^2);
        
%         if mod(iter,50) == 0
%             prompt = 'Get bored?';
%             goNext = input(prompt);
%             if goNext == 1
%                 iter = 9999;
%             end
%         end
        iter = iter + 1;
        iterCount = iterCount +1;
        
    end
%     pause
%     disp('enter any key to continue')
end

    
    %%
Press = Press(1:end-1,1)';
PressRef = PressRef(1:end-1,1)';

%%
IrPlant(1,:) = SinSim;
IrSvd(1,:) = PressRef;
IrSvdFb(1,:) = Press;
%i = i + 1;


%%
figure
plot(du)
%% 時間領域でのプ?ット
%SISO(iir)かMIMO(mpap)かで?C表示する?}を?u手動で?v?ﾘり替えてください
figure
set(gcf,'Position',[0 0 600 480])
normalize = 1;
if normalize == 1
    plot(Time,IrPlant(1,:)./max(IrPlant(1,:)),'k','linewidth',1);
    hold on
    %plot(Time,IrPlant(2,:)./max(IrPlant(2,:)),'k--','linewidth',1);
    
    plot(Time,IrSvd(1,:)./max(IrSvd(1,:)),'g','linewidth',2);%iir
    %plot(Time,IrSvd(2,:)./max(IrSvd(2,:)),'g--','linewidth',2);
    
    plot(Time,IrSvdFb(1,:)./max(IrSvdFb(1,:)),'r','linewidth',1);%
    %plot(Time,IrSvdFb(1,:)./max(IrSvdFb(1,:)),'r','linewidth',1);
else
    plot(Time,IrPlant(1,:),'k','linewidth',1);
    hold on
    plot(Time,IrSvd(1,:),'g','linewidth',2);%iir
    plot(Time,IrSvdFb(1,:),'r','linewidth',1);%
end
hold off
xlim([0,0.035])
set(gca,'FontSize',12)
xlabel('Time[sec]','FontSize',12)
ylabel('Amplitude','FontSize',12)
legend('Original','Original with noise','SVD-IIR','SVD-IIR with noise','SVD-IIR-Fb','SVD-IIR-Fb with Noise','Location','NorthEast')%iir
%legend('Original','SVD-MPAP','CAPZ-MPAP','Location','NorthEast');%mpap

%%?}の保存に使う%%
%filename=['IIR_tap50\t_r',num2str(k),];
%saveas(gcf,filename,'epsc')

%% 周波?迫ﾌ域でのプ?ット
%SISO(iir)かMIMO(mpap)かで?C表示する?}を?u手動で?v?ﾘり替えてください

figure
if SignalChange == 0
    %subplot(211)
    ylim([-10 5])
    xlim([20 200])
else
    xlim([20 1000])
end

hold off

if SignalChange == 0
    fftPlot2(IrPlant(1,:),fsSim,'b');
    hold on
    %fftPlot2(IrPlant(2,:),fsSim,'b--');

    set(gcf,'Position',[0 0 600 480])

    fftPlot2(IrSvd(1,:),fsSim,'g');
    %fftPlot2(IrSvd(2,:),fsSim,'g--');

    fftPlot2(IrSvdFb(1,:),fsSim,'r');
    %fftPlot2(IrSvdFb(2,:),fsSim,'r--');

    axis([20 1000 -15 5])
else
    fftPlot(IrPlant(1,:),fsSim,'b');
    hold on
    
    set(gcf,'Position',[0 0 600 480])

    fftPlot(IrSvd(1,:),fsSim,'g');
    fftPlot(IrSvdFb(1,:),fsSim,'r');

    axis([fInput-10 250 -60 0])
end

%xlim([20 200])
set(gca,'FontSize',12)
title('Frequency response')
xlabel('Frequency[Hz]','FontSize',12)
ylabel('Relative SPL[dB]','FontSize',12)
leg=legend('Original','Original with noise','SVD-IIR','SVD-IIR with noise','SVD-IIR-Fb','SVD-IIR-Fb with Noise','Location','Southwest');%iir
set(leg,'FontSize',10)