function [Z4R, Z4Rvect, Z4L, Z4Lvect] = mkDists(k, vMeshGrid, doFig4, indCortMag, BasisPath,scrsz)

start = vMeshGrid(1); step = vMeshGrid(2); stop = vMeshGrid(3);

[X,Y]=meshgrid(start:step:stop);

CM = 9.81.*(sqrt((X).^2+(Y).^2).^-0.83);
ctr = ceil(length(X)/2);
CM(ctr,ctr) = 100; %arbitrarily set origin to 100 instead of infinity

CM = CM./sum(CM(:)); %CM distribution

Z2R = 10*(1./(1+exp(-k*X))); % right hemifield dist
Z2R = Z2R./sum(Z2R(:));

Z3R = (CM.*Z2R);

Z2L = 10*(1-(1./(1+exp(-k*X)))); % left hemifield dist
Z2L = Z2L./sum(Z2L(:));

Z3L = (CM.*Z2L);

if indCortMag == 2 %if no CM, just hemifield distributions, skip step of combining with CM dist
    Z3R = Z2R;
    Z3L = Z2L;
end

Z4R = Z3R./sum(Z3R(:)); %% final-right, normalized
Z4Rvect = Z4R(:);

Z4L = Z3L./sum(Z3L(:)); % final-left, normalized
Z4Lvect = Z4L(:);

if doFig4
    fig = figure('Position',[scrsz(4)/3 1 scrsz(3)/2 scrsz(4)]);
    subplot(3,2,2), s1 = surf(X,Y,CM,'FaceAlpha',1); s1.EdgeColor = 'none';
    axis square, axis tight, colormap parula, colorbar
    subplot(3,2,3), s2 = surf(X,Y,Z2R,'FaceAlpha',1); s2.EdgeColor = 'none';
    axis square, axis tight, colormap parula, zlim([0 0.00002])
    subplot(3,2,4), s3 = surf(X,Y,Z2L,'FaceAlpha',1); s3.EdgeColor = 'none';
    axis square, axis tight, colormap parula, zlim([0 0.00002]),colorbar
    view(37.5,30)
    subplot(3,2,5), s4 = surf(X,Y,Z4R,'FaceAlpha',1); s4.EdgeColor = 'none';
    axis square, axis tight, colormap parula
    subplot(3,2,6), s5 = surf(X,Y,Z4L,'FaceAlpha',1); s5.EdgeColor = 'none';
    axis square, axis tight, colormap parula, colorbar
    view(37.5,30)

    print(fig, [[BasisPath '/Results/'], 'Figure4a'], '-dtiffn', '-r300')
end
