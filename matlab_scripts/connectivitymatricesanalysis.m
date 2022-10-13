%% Import data
connectome1k=importdata("__connmatrix1000.csv");
connectome5k=importdata("__connmatrix5000.csv");
connectome10k=importdata("__connmatrix10000.csv");
connectome50k=importdata("__connmatrix50000.csv");
connectome100k=importdata("__connmatrix100000.csv");
connectome500k=importdata("__connmatrix500000.csv");
connectome1M=importdata("__connmatrix1000000.csv");
connectome5M=importdata("__connmatrix5000000.csv");
connectome10M=importdata("__connmatrix10000000.csv");

connectomecontrol20=importdata("sub-control020_ses-midcycle_mrtrix_matrix.csv");
plotconmat(connectomecontrol20,1e7,"title")

%% Comparison of connectivity matrices mrtrix probablistic
f=figure();
subplot(2,3,1)
plotconmat("connmatrix5k.csv",5000,"5k probablistic")

subplot(2,3,2)
plotconmat("connmatrix10k.csv",10000,"10k probablistic")

subplot(2,3,3)
plotconmat("connmatrix50k.csv",50000,"50k probablistic")

subplot(2,3,4)
plotconmat("connmatrix100k.csv",100000,"100k probablistic")

subplot(2,3,5)
plotconmat("connmatrix.csv",1000000,"1M probablistic")

subplot(2,3,6)
plotconmat("__connmatrix1000000.csv",1000000,"1M deterministic")

sgtitle ("Comparison of matrices of mrtrix probablistic",'interpreter',...
    'latex','FontUnits','points','FontWeight','demi','FontSize',18,'FontName','Times')
f.WindowState = 'maximized';

%% Comparison of connectivity matrices mrtrix deterministic
f=figure();

subplot(2,2,1)
plotconmat("__connmatrix5000.csv",5000,"5k")

% subplot(2,4,2)
% plotconmat("__connmatrix10000.csv",10000,"10k")

subplot(2,2,2)
plotconmat("__connmatrix50000.csv",50000,"50k")

% subplot(2,4,4)
% plotconmat("__connmatrix100000.csv",100000,"100k")

% subplot(2,4,5)
% plotconmat("__connmatrix500000.csv",500000,"500k")

subplot(2,2,3)
plotconmat("__connmatrix1000000.csv",1000000,"1M")

% subplot(2,4,7)
% plotconmat("__connmatrix5000000.csv",5000000,"5M")

subplot(2,2,4)
plotconmat("__connmatrix10000000.csv",10000000,"10M")

sgtitle ("Comparison of matrices of mrtrix deterministic",'interpreter',...
    'latex','FontUnits','points','FontWeight','demi','FontSize',18,'FontName','Times')

f.WindowState = 'maximized';

%% Evolution of certain ROIxROI connectivity

f=figure();

i=1;
for n=1:115
    for m=n+1:116
        %n=randi(116);m=randi(116);if m==n,m=randi(116); end
        a=connectome10M(n,m)./1e7;
        value=[connectome1k(n,m)./1e3 connectome5k(n,m)./5e3 connectome10k(n,m)./1e4 ...
            connectome50k(n,m)./5e4 connectome100k(n,m)./1e5 connectome500k(n,m)./5e5 ...
            connectome1M(n,m)./1e6 connectome5M(n,m)./5e6 connectome10M(n,m)./1e7];
        value=value./a;
        if ~(any(value(:) == 0))
            i=i+1;
            nstreamlines=[1000 5000 10000 50000 100000 500000 1000000 5000000 10000000];
            loglog(nstreamlines,value,'DisplayName',"ROIS "+num2str(min(n,m))+" X "+num2str(max(n,m)));
            xlabel ('Number of streamlines')
            xlim([1000 1e7])
            ylabel ('Ratio')
            hold on
        end
    end
    
end
%legend('show')
%legend('Location','best')

sgtitle('Connectivity Convergence','interpreter','latex','FontUnits','points',...
    'FontWeight','demi','FontSize',18,'FontName','Times')
f.WindowState = 'maximized';
