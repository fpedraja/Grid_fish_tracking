addpath('D:\KIT3'); addpath('D:\KIT');
clearvars; %close all;
myKsDir = uigetdir('D:\datosgrid2022\recordings_03_23\'); % check folder
files2=dir([myKsDir, '\*.wav']); % check wav files in folder
%%
%di=1; AMP=[];
v = VideoWriter ('test_first_test3.avi'); open(v);
for i=320:60:size(files2,1) % files to read (1min each)
    
    try
        [y,Fs] = audioread([myKsDir,'\',files2(i).name]);
        
        tic
        Data=[]; Data=[Data;y]; clear y;
        [b1, a1] = butter(3, [300/Fs,2000/Fs]*2, 'bandpass');
        datr = filter(b1, a1, Data); clear Data;
        datr=datr.^2;
        
        figure(50); subplot(8,2,1:2:16)
        a=0;
        for t=1:8
            plot(datr(:,t)); hold on;
            a=a-0.004;
        end
        pause(5)
        close(50)
        %
        EODti=[];
        for j=1:size(datr,2)
            [~,sample1_M1]=findpeaks(datr(:,j) ,'MINPEAKHEIGHT',0.02,'MINPEAKDISTANCE',50);
            EODti=[EODti;sample1_M1];
            %                 Xamp=[]; Xamp=datr(sample1_M1,j);
            %                 if length(Xamp)>=200
            %                     AMP(di,j)=nanmedian(Xamp);
            %                 else
            %                     AMP(di,j)=nan;
            %                 end
            %
            %         subplot(8,2,j*2); area(sample1_M1,Xamp); xlim([0 1200000]); ylim([0 1]); %hold on;
            %         a=a+0.4;
        end
        pause(0.0001)
        %di=di+1;
        if isempty(EODti)==1
            disp(['Archivo ' num2str(i) ' ' [files2(i).name(7:14),' ', files2(i).name(16:17),':',files2(i).name(18:19),':',files2(i).name(20:21)] ' no hay peces' ])
        else
            EODtisort=sort(EODti); b=2;  AUX2(1,1)=0; AUX2(2:length(EODtisort),1)=diff(EODtisort);
            %figure;hist(AUX2,1000);
            EODtiend=[]; EODtiend(1,1)=EODtisort(1);
            
            for h=2:length(EODtisort)
                if AUX2(h)>=20
                    EODtiend(b,1)=EODtisort(h);
                    b=b+1;
                end
            end
            AUXC=datr(EODtiend,:); clear datr;
            
            %% spatial distribution
            %addpath('C:\Users\fedu1\My Drive\training')
            AUXM=nan(3,4); SDS=[]; DD=['b','r','g','y','k','c','m'];
            as=figure(100);
            for p=1:length(EODtiend)
                AUXM(1,1)=sqrt(abs(AUXC(p,1)))./sum(sqrt(abs(AUXC(p,:))));
                AUXM(1,3)=sqrt(abs(AUXC(p,2)))./sum(sqrt(abs(AUXC(p,:))));
                AUXM(2,1:4)=sqrt(abs(AUXC(p,3:6)))./sum(sqrt(abs(AUXC(p,:))));
                AUXM(3,2)=sqrt(abs(AUXC(p,7)))./sum(sqrt(abs(AUXC(p,:))));
                AUXM(3,4)=sqrt(abs(AUXC(p,8)))./sum(sqrt(abs(AUXC(p,:))));
                
                SDS(:,:,p)=AUXM;
                divi=1;
                MA_control=nan(80/divi,120/divi,1);
                
                MA_control(1,1)=SDS(1,1,p);
                MA_control(1,80/divi)=SDS(1,3,p);
                MA_control(40/divi,1)=SDS(2,1,p);
                MA_control(40/divi,40/divi)=SDS(2,2,p);
                MA_control(40/divi,80/divi)=SDS(2,3,p);
                MA_control(40/divi,120/divi)=SDS(2,4,p);
                MA_control(80/divi,40/divi)=SDS(3,2,p);
                MA_control(80/divi,120/divi)=SDS(3,4,p);
                
                MA_control(80/divi,1)=0;
                MA_control(1,40/divi)=0;
                MA_control(1,120/divi)=0;
                MA_control(80/divi,80/divi)=0;
                
                SDS3_1=inpaint_nans(MA_control(:,:));
                
                %subplot(1,3,1); contour(SDS3_1,100); hold on; contour(SDS3_1,[max(max(SDS3_1))*0.90 max(max(SDS3_1))],'LineColor', 'k'); axis equal; hold off
                %subplot(1,3,2); plot(Datazs(EODtiend(i)-50:EODtiend(i)+50,1:8));
                
                [AUXX, AUXY]=find(SDS3_1(:,:)>=max(max(SDS3_1))*0.9); X=round(mean(AUXX)); Y=round(mean(AUXY));
                subplot(1,2,1); plot(Y,X,'or'); axis equal; xlim([0 120/divi]); ylim([0 80/divi]); xticks([0:10:120]); yticks([0:10:80]); grid on;
                title([files2(i).name(7:14),' ', files2(i).name(16:17),':',files2(i).name(18:19),':',files2(i).name(20:21)])
                %                 subplot(1,2,2);
                %                 Color=5;
                %                 plot(Y,X,'.','MarkerSize',10,...
                %                     'MarkerEdgeColor',DD(Color),...
                %                     'MarkerFaceColor',DD(Color))
                %                 axis equal; xlim([0 120/divi]); ylim([0 80/divi]); hold on;
                
                pause(0.000000001)
                
                F=getframe(as);
                writeVideo (v, F);
                
            end
            
            disp(['Archivo ' num2str(i) ' ' [files2(i).name(7:14),' ', files2(i).name(16:17),':',files2(i).name(18:19),':',files2(i).name(20:21)] ' finalizado en ' num2str(toc)])
        end
    catch
        disp(['Archivo ' num2str(i) ' ' [files2(i).name(7:14),' ', files2(i).name(16:17),':',files2(i).name(18:19),':',files2(i).name(20:21)] ' defectuoso' ])
    end
end

close(v);

