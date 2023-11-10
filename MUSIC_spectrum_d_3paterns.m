%%%%%%%%% 2つ目のシミュレーションコード %%%%%%%%%%%

%1つ目のシミュレーション同様、疑似アンテナ３本を含むアンテナ数７、ユーザ数５の場合で到来方向推定を行う%
%その際に、アンテナ間隔による推定精度の違いを示したものである%

clear all;
close all hidden;
clc;

%% パラメータ設定
ss=100; %スナップショット
number_ant_z=4; %アンテナ行数
number_ant_x=4; %アンテナ列数
M=number_ant_z*number_ant_x; %総アンテナ数
number_virtual_ant_z=3; %疑似アンテナ行数
number_virtual_ant_x=3; %疑似アンテナ列数
sum_ant_z=number_ant_z+number_virtual_ant_z; %Z軸方向の疑似＋既存アンテナ数
sum_ant_x=number_ant_x+number_virtual_ant_x; %X軸方向の疑似＋既存アンテナ数
virtual_M=number_virtual_ant_z*number_virtual_ant_x; %疑似アンテナ総数
SNR=20; %SNR
sigma=10^(-SNR/10); %雑音分散
L=5; %到来波数
number_user=L; %ユーザ数
theta=[30 60 75 100 150]; %DOA（天頂角）
phi=[30 60 75 100 150]; %DOA（方位角）

derad=pi/180;
c=physconst('LightSpeed'); %光速
f0=100*10^6; %中心周波数
lambda=c/f0; %波長

angle_of_interval=1; %推定角度間隔
s=randn(number_user,ss); %ランダム信号（送信源）

%% アンテナ間隔のパターン（d=lambda/2,lambda,2*lambdaの３パターンで精度比較）
d_paterns=[lambda/2, lambda, 2*lambda];

%% MUSIC spectrum
%theta
plage_theta=(0:angle_of_interval:180);
ltheta=length(plage_theta);
image_theta=zeros(length(d_paterns),ltheta); %MUSIC spectrumの結果を格納

%phi
plage_phi=(0:angle_of_interval:180);
lphi=length(plage_phi);
image_phi=zeros(length(d_paterns),lphi);

%% アンテナ間隔３パターンでループ
for i=1:length(d_paterns)
    
    d=d_paterns(i);
    
    %% real antenna情報

    %%%%%====①天頂角θ推定（y-z平面） %%%%%%
    %既存アンテナの位置座標
    z1=0;
    z2=d;
    z3=2*d;
    z4=3*d;

    %ステアリングベクトル
    A_theta=zeros(number_ant_z,number_user);
    for nu=1:number_user
        psi_theta=2*pi*d*cosd(theta(nu))/lambda;
        a_psi_theta=exp(1i*((1:number_ant_z).'-1)*psi_theta);
        A_theta(:,nu)=a_psi_theta; 
    end

    %noise
    n_theta=sqrt(sigma)*(randn(number_ant_z,ss)+1i*randn(number_ant_z,ss));

    %受信信号
    x_theta=A_theta*s+n_theta;

    %%%%%====②方位角Φ推定（x-y平面）%%%%%%%
    %既存アンテナの位置座標
    x1=0;
    x2=d;
    x3=2*d;
    x4=3*d;

    %既存アンテナによる受信信号
    A_phi=zeros(number_ant_x,number_user);
    for nu=1:number_user
        psi_phi=2*pi*d*cosd(phi(nu))/lambda;
        a_psi_phi=exp(1i*((1:number_ant_x).'-1)*psi_phi);
        A_phi(:,nu)=a_psi_phi; %ステアリングベクトル（方向を特定するベクトル）
    end

    %noise
    n_phi=sqrt(sigma)*(randn(number_ant_x,ss)+1i*randn(number_ant_x,ss));

    %受信信号
    x_phi=A_theta*s+n_theta;

    %% 提案手法（疑似アンテナを用いたMUSIC法）
    for range_virtual=1:2
        if range_virtual==1

            %%%%%====①天頂角θ推定（y-z平面） %%%%%%
            %疑似アンテナの位置座標→既存アンテナの内分から３つの疑似アンテナ位置推定（アンテナの中心に疑似アンテナ配置）
            v_z1=(z1+z2)/(1+1);
            v_z2=(z2+z3)/(1+1);
            v_z3=(z3+z4)/(1+1);

            v_z=[v_z1;v_z2;v_z3]; %アンテナ位置を格納

            %疑似ステアリングベクトル
            v_A_theta=zeros(number_virtual_ant_z,number_user);
            for nu=1:number_user
                v_psi_theta=2*pi*v_z*cosd(theta(nu))/lambda;
                a_v_psi_theta=exp(1i*v_psi_theta);
                v_A_theta(:,nu)=a_v_psi_theta;
            end

            %疑似ノイズ
            variances=var(n_theta,0,2); %ノイズ分散（from real antenna）
            mean_variances=mean(variances); %ノイズ分散の平均

            v_n_theta=sqrt(mean_variances)*(randn(number_virtual_ant_z,ss)+1i*randn(number_virtual_ant_z,ss));

            v_x_theta=v_A_theta*s+v_n_theta; %疑似受信信号
            X_theta=[x_theta;v_x_theta]; %既存＋疑似の受信信号を結合
            Rx_theta=X_theta*X_theta'/ss; %相関行列

            %固有値展開
            [EV,D]=eig(Rx_theta);
            EVA=diag(D)';
            [EVA,I]=sort(EVA);
            EVA=fliplr(EVA);
            EV=fliplr(EV(:,I));

            %雑音部分空間の抽出
            En_theta=EV(:,number_user+1:number_ant_z+number_virtual_ant_z);

            %MUSIC spectrum search（天頂角推定）
            for itheta=1:ltheta
                th=plage_theta(itheta);
                psi_search_th=2*pi*d*cosd(th)/lambda;
                a_psi_search_th=exp(1i*((1:number_ant_z).'-1)*psi_search_th); %既存アンテナのステアリングベクトル
                v_a_psi_search_th=exp(1i*2*pi*v_z*cosd(th)/lambda);
                sum_a_psi_search_th=[a_psi_search_th; v_a_psi_search_th]; %疑似と既存の合わせたステアリングベクトル
                image_theta(i,itheta)=-10*log10(abs(sum_a_psi_search_th'*(En_theta*En_theta')*sum_a_psi_search_th)); %MUSIC spectrum
            end

        elseif range_virtual==2
            %%%%%====②方位角Φ推定（y-z平面） %%%%%%
            %疑似アンテナの位置座標→既存アンテナの内分から３つの疑似アンテナ位置推定
            v_x1=(z1+z2)/(1+1);
            v_x2=(z2+z3)/(1+1);
            v_x3=(z3+z4)/(1+1);

            v_x=[v_x1;v_x2;v_x3];

            %疑似ステアリングベクトル
            v_A_phi=zeros(number_virtual_ant_x,number_user);
            for nu=1:number_user
                v_psi_phi=2*pi*v_x*cosd(phi(nu))/lambda;
                a_v_psi_phi=exp(1i*v_psi_phi);
                v_A_phi(:,nu)=a_v_psi_phi;
            end

            %疑似ノイズ
            variances=var(n_phi,0,2); %ノイズ分散（from real antenna）
            mean_variances=mean(variances); %ノイズ分散の平均

            v_n_phi=sqrt(mean_variances)*(randn(number_virtual_ant_x,ss)+1i*randn(number_virtual_ant_x,ss));

            v_x_phi=v_A_phi*s+v_n_phi; %疑似受信信号
            X_phi=[x_phi;v_x_phi]; %既存＋疑似の受信信号を結合
            Rx_phi=X_phi*X_phi'/ss; %相関行列

            %固有値展開
            [EV,D]=eig(Rx_phi);
            EVA=diag(D)';
            [EVA,I]=sort(EVA);
            EVA=fliplr(EVA);
            EV=fliplr(EV(:,I));

            %雑音部分空間の抽出
            En_phi=EV(:,number_user+1:number_ant_x+number_virtual_ant_x);

            %MUSIC spectrum search（天頂角推定）
            for iphi=1:lphi
                ph=plage_phi(iphi);
                psi_search_ph=2*pi*d*cosd(ph)/lambda;
                a_psi_search_ph=exp(1i*((1:number_ant_x).'-1)*psi_search_ph); %既存アンテナのステアリングベクトル
                v_a_psi_search_ph=exp(1i*2*pi*v_x*cosd(ph)/lambda);
                sum_a_psi_search_ph=[a_psi_search_ph; v_a_psi_search_ph]; %疑似と既存の合わせたステアリングベクトル
                image_phi(i,iphi)=-10*log10(abs(sum_a_psi_search_ph'*(En_phi*En_phi')*sum_a_psi_search_ph)); %MUSIC spectrum
            end
        end
    end
end

%% グラフ
figure(1)

%d=lambda/2,theta
subplot(3,2,1)
plot(plage_theta,image_theta(1,:),'r');
xlim([0 180]);
xlabel('elevation angle \theta [deg.]');
ylabel('MUSIC spectrum [dB]');
title('$d=\frac{\lambda}{2},d^\prime=\frac{\lambda}{4}$', 'Interpreter', 'latex');
grid on;

%d=lambda/2,phi
subplot(3,2,2)
plot(plage_phi,image_phi(1,:),'r');
xlim([0 180]);
xlabel('azimuth angle \phi [deg.]');
ylabel('MUSIC spectrum [dB]');
title('$d=\frac{\lambda}{2},d^\prime=\frac{\lambda}{4}$', 'Interpreter', 'latex');
grid on;

%d=lambda,theta
subplot(3,2,3)
plot(plage_theta,image_theta(2,:),'k');
xlim([0 180]);
xlabel('elevation angle \theta [deg.]');
ylabel('MUSIC spectrum [dB]');
title('$d=\lambda,d^\prime=\frac{\lambda}{2}$', 'Interpreter', 'latex');
grid on;

%d=lambda,phi
subplot(3,2,4)
plot(plage_phi,image_phi(2,:),'k');
xlim([0 180]);
xlabel('azimuth angle \phi [deg.]');
ylabel('MUSIC spectrum [dB]');
title('$d=\lambda,d^\prime=\frac{\lambda}{2}$', 'Interpreter', 'latex');
grid on;

%d=2*lambda,theta
subplot(3,2,5)
plot(plage_theta,image_theta(3,:),'b');
xlim([0 180]);
xlabel('elevation angle \theta [deg.]');
ylabel('MUSIC spectrum [dB]');
title('$d=2\lambda,d^\prime=\lambda$', 'Interpreter', 'latex');
grid on;

%d=2*lambda,phi
subplot(3,2,6)
plot(plage_phi,image_phi(3,:),'b');
xlim([0 180]);
xlabel('azimuth angle \phi [deg.]');
ylabel('MUSIC spectrum [dB]');
title('$d=2\lambda,d^\prime=\lambda$', 'Interpreter', 'latex');
grid on;
