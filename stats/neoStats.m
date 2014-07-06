


%% Compute ROC 

% Create true labels vector
labels = [ones(1,size(sSI,2)),zeros(1,size(sSII,2))];

% Tabulate scores
scores = [sSI,sSII];

% Define positive class
posclass = 1;

% Generate AUC
[X,Y,T,AUC,OPTROCPT] = perfcurve(labels,scores,posclass);
optTPR = OPTROCPT(1,2); %optimal TPR

% more stats
%[X,Y,T,AUC] = perfcurve(labels,scores,posclass,'Xcrit','accu');

% Plot ROC Curve
figure;
plot(X,Y)
xlabel('False positive rate'); ylabel('True positive rate')
title('ROC, using mean landmarks per top 200 matches, basis X-axis')
annotation('textbox',...
    [0.3 0.65 0.2 0.1],...
    'String',{'AUC:',num2str(AUC)},...
    'FontSize',14,...
    'FontName','Arial',...
    'LineStyle','--',...
    'EdgeColor',[0 0 0],...
    'LineWidth',1,...
    'BackgroundColor',[0.9  0.9 0.9],...
    'Color',[0.84 0.16 0]);
axis square; set(gcf,'color','w');
hold on
scatter(OPTROCPT(1,1),OPTROCPT(1,2),80,'*','r');

% Find optimal threshold
[TPR,Y,T,AUC] = perfcurve(labels,scores,posclass,'Xcrit','TPR');
ind = find(TPR==optTPR);
thresh = T(ind);

%% Plot threshold
figure
scatter(1:size(scores,2),scores,'.');
title('searchScore');
set(gcf,'color','w');
hold on
plot(get(gca,'xlim'), [thresh thresh]);
hold off 

