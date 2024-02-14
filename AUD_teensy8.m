addpath('D:\KIT3'); addpath('D:\KIT');
clearvars; %close all;
myKsDir = uigetdir('D:\datosgrid2022\10 a 11_marzo 2023\'); % check folder
files2=dir([myKsDir, '\*.wav']); % check wav files in folder
%%
%di=1; AMP=[];
%v = VideoWriter ('test_first_test3.avi'); open(v);
%%
for i=1280:40:size(files2,1) % files to read (1min each)
    
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
%                 figure; subplot(8,2,1:2:16)
%                 a=0;
%                 for t=1:8
%                     plot(datr(:,t)); hold on;
%                     a=a-0.004;
%                 end
%                 pause(1)
                %close(50)
        %
        EODti=[];
        for j=1:size(datr,2)
            [~,sample1_M1]=findpeaks(datr(:,j) ,'MINPEAKHEIGHT',1e-4,'MINPEAKDISTANCE',50);
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
                X=[]; Y=[]; sigma=0.999; min_x=0; min_y=0; max_x=80; max_y=120; step=1;
                    % Define the grid of candidate emitter locations
                    % where min_x, min_y, max_x, max_y define the boundaries of the search area and step is the grid spacing
                [xgrid, ygrid] = meshgrid(min_x:step:max_x, min_y:step:max_y);
              
                for p=1:length(EODtiend) %nthroot(-27, 3)
                    AUXM=nan(3,4); SDS=[];
                    AUXM(1,1)=nthroot(abs(AUXC(p,1)),4)./sum(nthroot(abs(AUXC(p,:)),4));
                    AUXM(1,3)=nthroot(abs(AUXC(p,2)),4)./sum(nthroot(abs(AUXC(p,:)),4));
                    AUXM(2,1:4)=nthroot(abs(AUXC(p,3:6)),4)./sum(nthroot(abs(AUXC(p,:)),4));
                    AUXM(3,2)=nthroot(abs(AUXC(p,7)),4)./sum(nthroot(abs(AUXC(p,:)),4));
                    AUXM(3,4)=nthroot(abs(AUXC(p,8)),4)./sum(nthroot(abs(AUXC(p,:)),4));
                    
                    AUX3_1=inpaint_nans(AUXM); % 

                    %measurements = [0 0 AUX3_1(1,1); 0 40 AUX3_1(1,2); 0 80 AUX3_1(1,3); 0 120 AUX3_1(1,4); 40 0 AUX3_1(2,1); 40 40 AUX3_1(2,2); 40 80 AUX3_1(2,3); 40 120 AUX3_1(2,4); 80 0 AUX3_1(3,1); 80 40 AUX3_1(3,2); 80 80 AUX3_1(3,3); 80 120 AUX3_1(3,4)]; % where x,y are the coordinates of each measurement location and s is the corresponding signal strength
                     measurements = [0 0 AUX3_1(1,1); 0 80 AUX3_1(1,3); 40 0 AUX3_1(2,1); 40 40 AUX3_1(2,2); 40 80 AUX3_1(2,3); 40 120 AUX3_1(2,4); 80 40 AUX3_1(3,2); 80 120 AUX3_1(3,4)]; % where x,y are the coordinates of each measurement location and s is the corresponding signal strength
                    
        
                    
                    % Compute the weights based on the distance between each measurement location and each candidate emitter location
                    dists = pdist2(measurements(:,1:2), [xgrid(:) ygrid(:)]);
                    weights = exp(-dists.^2 / (2 * sigma^2)); % where sigma is a parameter that controls the shape of the weight function from ~0 to 1
                    %> sigma better will result in smoother weight functions and may lead to more robust estimates but with lower spatial resolution.
                    %< sigma will result in sharper weight functions and may lead to higher spatial resolution but with higher sensitivity to noise and measurement errors.
                    
                    % Compute the weighted spatial average of the signal strengths at each candidate emitter location
                    weighted_sum = bsxfun(@times, measurements(:,3), weights);
                    weighted_avg = sum(weighted_sum, 1) ./ sum(weights, 1);
                    
                    % Find the location with the highest weighted average as the estimated position of the emitter
                    [~, idx] = max(weighted_avg);
                    X(p) = xgrid(idx);
                    Y(p) = ygrid(idx);
                end
                
%          figure; plot(Y,X,'.')
%          axis equal; xlim([0 120]); ylim([0 80]); xticks([0:10:120]); yticks([0:10:80]); grid on;
                
                
                % Choose the maximum number of clusters to consider
                max_k = 10;
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
                
                [jo,je]=sort(freq,'descend'); centroidsNEW=centroids(je,:);
                
                DD=['b','r','g','y','k','c','m'];
                as=figure(100); 
                for p=1:size(centroidsNEW,1)
                    subplot(1,3,1); plot(centroidsNEW(p,2),centroidsNEW(p,1),'.', 'MarkerSize',10,...
                        'MarkerEdgeColor',DD(p),...
                        'MarkerFaceColor',DD(p))
                    
                    axis equal; xlim([0 120]); ylim([0 80]); xticks([0:10:120]); yticks([0:10:80]); grid on;
                    title([files2(i).name(7:14),' ', files2(i).name(16:17),':',files2(i).name(18:19),':',files2(i).name(20:21)])
                    subplot(1,3,2);
                    
                    plot(centroidsNEW(p,2),centroidsNEW(p,1),'.','MarkerSize',10,...
                        'MarkerEdgeColor',DD(p),...
                        'MarkerFaceColor',DD(p))
                    axis equal; xlim([0 120]); ylim([0 80]); xticks([0:10:120]); yticks([0:10:80]); grid on; hold on;
                    
                    subplot(1,3,3);  plot(i,jo(p),'.','MarkerSize',10,...
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

%close(v);