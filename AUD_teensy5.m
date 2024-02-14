addpath('D:\KIT3'); addpath('D:\KIT');
clearvars; %close all;
myKsDir = uigetdir('D:\datosgrid2022\10 a 11_marzo 2023\'); % check folder
files2=dir([myKsDir, '\*.wav']); % check wav files in folder
%%
di=1; 
for i=4880:40:5280%size(files2,1)
    AMP=[];
    try
        [y,Fs] = audioread([myKsDir,'\',files2(i).name]);
        tic
        Data=[]; Data=[Data;y];
        [b1, a1] = butter(3, [300/Fs,2000/Fs]*2, 'bandpass');
        datr = filter(b1, a1, Data);
        AUX=median(datr,1); datr=datr-AUX;
        %datr=datr.^2;
        figure; subplot(8,2,1:2:16); 
        a=0;
        for t=1:8
            plot(datr(:,t)+a); hold on;
            a=a-0.1;
        end
        hold on; title([files2(i).name(7:14),' ', files2(i).name(16:17),':',files2(i).name(18:19),':',files2(i).name(20:21)])
        Datazs=[];  EODti=[];
        a=0;
        for j=1:size(datr,2)
            [~,sample1_M1]=findpeaks(datr(:,j) ,'MINPEAKHEIGHT',0.00020,'MINPEAKDISTANCE',70);
            EODti=[EODti;sample1_M1];
            Xamp=[]; Xamp=datr(sample1_M1,j);
            if length(Xamp)>=60
                AMP(i,j)=nanmedian(Xamp);
            else
                AMP(i,j)=nan;
            end
            
            subplot(8,2,j*2); area(sample1_M1,Xamp); xlim([0 1200000]); ylim([0 1]); %hold on;
            a=a+0.4;
        end
        pause(0.0001)
        if nansum(AMP(i,:))>0
            disp(['Archivo ' num2str(i) ' ' [files2(i).name(7:14),' ', files2(i).name(16:17),':',files2(i).name(18:19),':',files2(i).name(20:21)] ' finalizado en ' num2str(toc)])
        else
            disp(['Archivo ' num2str(i) ' ' [files2(i).name(7:14),' ', files2(i).name(16:17),':',files2(i).name(18:19),':',files2(i).name(20:21)] ' no hay peces' ])
        end
        
   
    catch
        AMP(i,1:8)=nan;
        disp(['Archivo ' num2str(i) ' ' [files2(i).name(7:14),' ', files2(i).name(16:17),':',files2(i).name(18:19),':',files2(i).name(20:21)] ' defectuoso' ])
    end
    
end


