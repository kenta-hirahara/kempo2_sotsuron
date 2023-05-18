function dispersionFigNum = dispersionPlot(kxkyt, cv, wc, newDirAbsolutePath, inputParam, pltColor, EBtex, EBtexCat, startSimulationDatetime, saveName, dispersionFigNum)
    kxkyw = abs(fft(kxkyt, inputParam.ntime, 3)/inputParam.ntime*2);
    clear kxkyt;
    kxkyw = log10(kxkyw);
    kxkyw = fftshift(kxkyw, 2);
    kxkyw = fftshift(kxkyw, 3);

    wmax = 0.5*pi;  
    kmax = pi;
    dkx = 2*kmax/inputParam.nx;
    dky = 2*kmax/inputParam.ny;
    dw = 2*wmax/inputParam.ntime;
    w_axis = [-size(kxkyw, 3)/2+1:size(kxkyw, 3)/2] * dw / abs(wc);
    kx_axis = [0:size(kxkyw, 1)-1] * dkx * cv / abs(wc);
    ky_axis = [-size(kxkyw, 2)/2+1:size(kxkyw, 2)/2] * dky * inputParam.cv / abs(wc);
    
    parameterNames = {'kxkyw_Ex.mat', 'kxkyw_Ey.mat', 'kxkyw_Ez.mat', 'kxkyw_Bx.mat', 'kxkyw_By.mat', 'kxkyw_Bz.mat'};
    save(cell2mat(parameterNames(saveName)), 'startSimulationDatetime','inputParam','pltColor','EBtex','saveName'...
    ,'wmax','kmax','dkx','dky','dw','w_axis','kx_axis','ky_axis','kxkyw', '-v7.3');
    movefile(cell2mat(parameterNames(saveName)), newDirAbsolutePath);

    fig = figure(dispersionFigNum);
    fig.Name = 'Dispersion Relation';
    fig.Position = [0, 100, 1000, 350];

    ax(1) = subplot(1,2,1);
    im = imagesc(kx_axis, w_axis, squeeze(kxkyw(:,size(kxkyw, 2)/2,:))');
    colormap(pltColor.map); c(1) = colorbar; shading flat;
    caxis([-9, -2]);
    c(1).Title.FontSize = inputParam.Fontsize;
    c(1).Label.Interpreter = 'latex';
    c(1).Label.String = ['$\log_{10}|' cell2mat(EBtexCat(saveName)) '|$'];
    ax(1).YDir='normal';
    ax(1).Title.Interpreter = 'latex';
    ax(1).Title.FontSize = inputParam.Fontsize;
    ax(1).Title.String = cell2mat(EBtex(saveName));
    ax(1).XLabel.Interpreter = 'latex';
    ax(1).XLabel.FontSize = inputParam.Fontsize;
    ax(1).XLabel.String = '$k_{x}c\Omega_{e}^{-1}$';
    ax(1).YLabel.Interpreter = 'latex';
    ax(1).YLabel.FontSize = inputParam.Fontsize;
    ax(1).YLabel.String = '$\omega\Omega_{e}^{-1}$';
    ax(1).XLim = [0,10];
    ax(1).YLim = [0,10];

    ax(2) = subplot(1,2,2); im = imagesc(ky_axis, w_axis, squeeze(kxkyw(1,:,:))');
    colormap(pltColor.map); c(2) = colorbar; shading flat;
    % caxis([-9, -2]);
    c(1).Title.FontSize = inputParam.Fontsize;
    c(2).Label.Interpreter = 'latex';
    c(2).Label.String = ['$\log_{10}|' cell2mat(EBtexCat(saveName)) '|$'];
    ax(2).YDir='normal';
    ax(2).Title.Interpreter = 'latex';
    ax(2).Title.String = cell2mat(EBtex(saveName));
    ax(2).Title.FontSize = inputParam.Fontsize;
    ax(2).XLabel.Interpreter = 'latex';
    ax(2).XLabel.String = '$k_{y}c\Omega_{e}^{-1}$';
    ax(2).XLabel.FontSize = inputParam.Fontsize;
    ax(2).YLabel.Interpreter = 'latex';
    ax(2).YLabel.String = '$\omega\Omega_{e}^{-1}$';
    ax(2).YLabel.FontSize = inputParam.Fontsize;
    ax(2).XLim = [0,10];
    ax(2).YLim = [0,10];

    dispersionParallelFigName = strcat('dispersionParallel_', cell2mat(parameterNames(saveName)), '.fig');
    savefig(fig, dispersionParallelFigName);
    movefile(dispersionParallelFigName, newDirAbsolutePath);
    dispersionFigNum = dispersionFigNum + 1;
end