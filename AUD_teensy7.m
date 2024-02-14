addpath('D:\KIT3'); addpath('D:\KIT');
clearvars; %close all;
myKsDir = uigetdir('D:\datosgrid2022\10 a 11_marzo 2023\'); % check folder
files2=dir([myKsDir, '\*.wav']); % check wav files in folder
%%
%di=1; AMP=[];
%v = VideoWriter ('test_first_test3.avi'); open(v);
%%
for i=1:1:size(files2,1) % files to read (1min each or 15 seconds?)
    
    y=[];
    try %#ok<TRYNC>
        [y,Fs] = audioread([myKsDir,'\',files2(i).name]);
    end
    
    if isempty(y)==0
        tic
        Data=[]; Data=[Data;y]; clear y;
        [b1, a1] = butter(3, [300/Fs,2000/Fs]*2, 'bandpass');
        datr = filter(b1, a1, Data); clear Data;
        datr=datr.^2;
        %
                figure(50); subplot(8,2,1:2:16)
                a=0;
                for t=1:8
                    plot(datr(:,t)); hold on;
                    a=a-0.004;
                end
                pause(1)
                close(50)
        %
        EODti=[];
        for j=1:size(datr,2)
            [~,sample1_M1]=findpeaks(datr(:,j) ,'MINPEAKHEIGHT',0.028,'MINPEAKDISTANCE',50);
            EODti=[EODti;sample1_M1];
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
            
            if size(EODtiend,1)<30
                disp(['Archivo ' num2str(i) ' ' [files2(i).name(7:14),' ', files2(i).name(16:17),':',files2(i).name(18:19),':',files2(i).name(20:21)] ' no hay peces' ])
            else
                divi=5;  X=[]; Y=[];
                for p=1:length(EODtiend) %nthroot(-27, 3)
                    AUXM=[]; SDS=[]; MA_control=nan(80/divi,120/divi,1);
                    AUXM(1,1)=nthroot(abs(AUXC(p,1)),4)./sum(nthroot(abs(AUXC(p,:)),4));
                    AUXM(1,3)=nthroot(abs(AUXC(p,2)),4)./sum(nthroot(abs(AUXC(p,:)),4));
                    AUXM(2,1:4)=nthroot(abs(AUXC(p,3:6)),4)./sum(nthroot(abs(AUXC(p,:)),4));
                    AUXM(3,2)=nthroot(abs(AUXC(p,7)),4)./sum(nthroot(abs(AUXC(p,:)),4));
                    AUXM(3,4)=nthroot(abs(AUXC(p,8)),4)./sum(nthroot(abs(AUXC(p,:)),4));
                    
                    SDS(:,:,p)=AUXM;
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
                    
                    %subplot(1,3,1); contour(SDS3_1,100); hold on; contour(SDS3_1,[max(max(SDS3_1))*0.95 max(max(SDS3_1))],'LineColor', 'k'); axis equal; hold off
                    %subplot(1,3,2); plot(Datazs(EODtiend(i)-50:EODtiend(i)+50,1:8));
                    
                    [AUXX, AUXY]=find(SDS3_1(:,:)>=max(max(SDS3_1))*0.7); X(p)=(mean(AUXX))* divi; Y(p)=(mean(AUXY))* divi;
                end
                
                
                
                
                
                % Choose the maximum number of clusters to consider
                max_k = 7;
                % Compute the within-cluster sum of squares (WCSS) for different values of k
                wcss = [];
                for k = 1:max_k
                    [~,C,sumd] = kmeans([X; Y]', k, 'EmptyAction', 'singleton');
                    wcss(k)=sum(sumd);
                end
                % Plot the WCSS as a function of the number of clusters
                %plot(1:max_k, wcss); xlabel('Number of clusters'); ylabel('Within-cluster sum of squares');
                % Find the elbow point
                diff_wcss = diff(wcss);
                second_derivative = diff(diff_wcss);
                elbow_point = find(second_derivative > 0, 1, 'first') + 1;
                % Perform k-means clustering using the optimal number of clusters
                num_clusters = elbow_point;
                
                [idx, centroids] = kmeans([X; Y]', num_clusters); % Plot the clustering results scatter(your_data(:,1), your_data(:,2), 10, idx, 'filled'); xlabel('Feature 1'); ylabel('Feature 2'); title(['Clustering Results (', num2str(num_clusters), ' clusters)']);
                
                try
                    for k=1:num_clusters
                        if size(idx(idx==k),1)<=100
                            AUXC(idx==k,:)=[];
                            centroids(k,:)=[];
                            idx(idx==k)=[];
                            
                        end
                    end
                end
                freq=[];
                for j=1:size(centroids,1)
                    freq(j)=nanmedian(1./(diff(EODtiend(idx==j))/Fs));
                end
                
                DD=['b','r','g','y','k','c','m'];
                as=figure(100);
                for p=1:size(centroids,1)
                    subplot(1,3,1); plot(centroids(p,2),centroids(p,1),'.', 'MarkerSize',10,...
                        'MarkerEdgeColor',DD(p),...
                        'MarkerFaceColor',DD(p))
                    
                    axis equal; xlim([0 120]); ylim([0 80]); xticks([0:10:120]); yticks([0:10:80]); grid on;
                    title([files2(i).name(7:14),' ', files2(i).name(16:17),':',files2(i).name(18:19),':',files2(i).name(20:21)])
                    subplot(1,3,2);
                    
                    plot(centroids(p,2),centroids(p,1),'.','MarkerSize',10,...
                        'MarkerEdgeColor',DD(p),...
                        'MarkerFaceColor',DD(p))
                    axis equal; xlim([0 120]); ylim([0 80]); xticks([0:10:120]); yticks([0:10:80]); grid on; hold on;
                    
                    subplot(1,3,3);  plot(i,freq(p),'.','MarkerSize',10,...
                        'MarkerEdgeColor',DD(p),...
                        'MarkerFaceColor',DD(p)); ylim([0 80]); hold on;
                    
                    pause(0.000000001)
                    
                    %                 F=getframe(as);
                    %                 writeVideo (v, F);
                    
                end
                disp(['Archivo ' num2str(i) ' ' [files2(i).name(7:14),' ', files2(i).name(16:17),':',files2(i).name(18:19),':',files2(i).name(20:21)] ' finalizado en ' num2str(toc)])
            end
        end
    else
        disp(['Archivo ' num2str(i) ' ' [files2(i).name(7:14),' ', files2(i).name(16:17),':',files2(i).name(18:19),':',files2(i).name(20:21)] ' defectuoso' ])
    end
end

close(v);