function  figNum = oblique(nkx, nky, kxkyw, kx_axis, w_axis, inputParam, pltColor, EBstring, EB, figNum, app.UIAxes)
    kx_div_ky = nkx/abs(nky); % int/int
    ky_div_kx = nky/nkx; % int/int
    % krw = zeros(size(kxkyw, 1), size(kxkyw, 3)/2+1); %wの正の部分のみ
    numOfPoints = min(ceil(size(kxkyw,1)/nkx), ceil(size(kxkyw,1)/abs(nky)));
    krw = zeros(numOfPoints, size(kxkyw, 3)); %wの正の部分のみ
    for w=1:size(kxkyw, 3)
      for i=1:numOfPoints
        krw(i, w) = kxkyw(1+(i-1)*nkx, size(kxkyw,1)+(i-1)*nky-1, w);
      end
    end
    
    k_norm = sqrt(nkx^2+nky^2)/nkx;
    fig  = figure(figNum);
    fig.Name = 'Oblique Dispersion Relation';
    fig.Position = [0, 100, 800, 600];
    ax = app.UIAxes;
    im_krw = imagesc(kx_axis(1:nkx:numOfPoints*nkx)*k_norm, w_axis, krw');
    colormap(pltColor.map); c = colorbar; shading flat;
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = inputParam.Fontsize;
    c.Label.String = ['$\log_{10}|' cell2mat(EBstring(EB.number)) '|$'];
    ax.XLabel.Interpreter = 'latex';
    ax.XLabel.String = '$k$';
    ax.XLabel.FontSize = inputParam.Fontsize;
    ax.YLabel.Interpreter = 'latex';
    ax.YLabel.String = '$\omega$';
    ax.YLabel.FontSize = inputParam.Fontsize;
    ax.YDir = 'normal';
    ax.Title.String = sprintf('%s \n %3.3f degree oblique mode', cell2mat(EBstring(EB.number)), rad2deg(atan(ky_div_kx))-inputParam.phi);
    ax.XLim = [0,25];
    ax.YLim = [0,25];
    figNum = figNum + 1;
  end