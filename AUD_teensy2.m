addpath('D:\KIT3'); addpath('D:\KIT');
clearvars; %close all;
myKsDir = uigetdir('C:\Users\fedu1\Downloads\teensy_2022\'); % check folder
files2=dir([myKsDir, '\*.wav']); % check wav files in folder
%%
a=1;  Xamp=[]; EODtiend2=[];
Data=[]; Datazs=[];
for i=a:a+19
    try
        [y,Fs] = audioread([myKsDir,'\',files2(i).name]);
        Data=[Data;y];
    catch
        disp(['Archivo ' num2str(i) ' defectuoso' ])
    end
end

[b1, a1] = butter(3, [300/Fs,3000/Fs]*2, 'bandpass');

datr = filter(b1, a1, Data);
Data = datr;

EODti=[];
for j=1:size(Data,2)
    Datazs(:,j)=zscore(Data(:,j));
    [~,sample1_M1]=findpeaks(Datazs(:,j) ,'MINPEAKHEIGHT',1.2,'MINPEAKDISTANCE',70);
    EODti=[EODti;sample1_M1];
end

clear Data datr

EODtisort=sort(EODti); b=2;  AUX2(1,1)=0; AUX2(2:length(EODtisort),1)=diff(EODtisort); figure;hist(AUX2,1000);
EODtiend=[]; EODtiend(1,1)=EODtisort(1);

for i=2:length(EODtisort)
    if AUX2(i)>=20
        EODtiend(b,1)=EODtisort(i);
        b=b+1;
    end
end
Xamp=[Xamp; Datazs(EODtiend,:)];
%% identity 
% epsilon=0.5;
% MinPts=10;
% [IDX, isnoise]=DBSCAN(X,epsilon,MinPts);
 
prompt = 'the best cluster number? ';
result = input(prompt);

[IDX,C,sumD] = kmeans( Xamp, result,... %basic variables data and number of classes
    'distance','sqEuclidean',... %distance measure set to squared euclidean distances
    'maxiter',1000,... % how many iterations before the calculation will be aboarded
    'start','cluster',... % starting positions determined by 10% pre clustering
    'display','final',... % display the number of iterations and the sum of the distance for each reptition
    'replicates', 10,... % The best of 10 repetitions is chosen
    'EmptyAction','singleton'); 

% check that the identity is ok

DD=['b','r','g','y','k','c','m'];
figure; a=0;
 for t=1:8
    plot(Datazs(1:100000,t)+a); hold on;
    a=a+4;
 end
 hold on;
 
 for i=1:size(IDX,1)
     plot(EODtiend(i),1,'.','MarkerSize',30,...
         'MarkerEdgeColor',DD(IDX(i)),...
         'MarkerFaceColor',DD(IDX(i)))
     hold on;
 end

figure;
for h=1:max(IDX)
    AUX=diff(EODtiend(IDX==h)); AUX2=(1./AUX)*Fs; AUX2(end+1)=AUX2(end);
    plot(EODtiend(IDX==h),AUX2,'.','MarkerSize',15,...
         'MarkerEdgeColor',DD(h),...
         'MarkerFaceColor',DD(h))
     hold on;    
end


%%
AUXM=[]; SDS=[]; DD=['b','r','g','y','k','c','m'];
figure(100)
for i=2:1600 %(EODtiend2)
    AUXM(1,1:4)=sqrt(abs(Xamp(i,1:4)))./sum(sqrt(abs(Xamp(i,1:8))));
    AUXM(2,1:4)=sqrt(abs(Xamp(i,[8 7 6 5])))./sum(sqrt(abs(Xamp(i,1:8))));
    w     = 2;   % Size of the sliding window (same number of cols and rows in this case)
    % Extrapolate values for current window
    [Nr,Nc] = size(AUXM);
    Nextra  = 0.5*(w-1);
    Ap      = interp2(1:Nc,1:Nr,AUXM,-Nextra+1:Nc+Nextra,(-Nextra+1:Nr+Nextra).','makima');    % 2D extrapolation must use 'spline' or 'makima' interpolation
    % Smooth data with sliding window
    H  = ones(w)./w^2;                      % The 2D averaging filter
    SDS(:,:,i)  = filter2(H,Ap,'valid');
    MA_control=nan(60,120,1); 
    MA_control(1,1)=SDS(1,1,i); MA_control(1,40)=SDS(1,2,i); MA_control(1,80)=SDS(1,3,i); MA_control(1,120)=SDS(1,4,i);
    MA_control(60,1)=SDS(2,1,i); MA_control(60,40)=SDS(2,2,i); MA_control(60,80)=SDS(2,3,i); MA_control(60,120)=SDS(2,4,i);
    SDS3_1=inpaint_nans(MA_control(:,:)); 
   
    subplot(1,3,1); contour(SDS3_1,100); hold on; contour(SDS3_1,[max(max(SDS3_1))*0.90 max(max(SDS3_1))],'LineColor', 'k'); hold off
    %subplot(1,3,2); plot(Datazs(EODtiend(i)-50:EODtiend(i)+50,1:8));   
    
    [AUXX, AUXY]=find(SDS3_1(:,:)>=max(max(SDS3_1))*0.9); X=round(mean(AUXX)); Y=round(mean(AUXY)); 
    subplot(1,3,2); plot(Y,X,'or'); xlim([0 120]); ylim([0 60]); 
    subplot(1,3,3);
            Color=IDX(i);
    plot(Y,X,'.','MarkerSize',10,...
    'MarkerEdgeColor',DD(Color),...
    'MarkerFaceColor',DD(Color))
    xlim([0 120]); ylim([0 60]); hold on;
    
   %pause(0.0001)   
end