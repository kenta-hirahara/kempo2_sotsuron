videoDuration = 40;
minFramerate = 20;
if check.xyPlot
  xyExyBvideoFileName = strcat(fileWithoutDotM, '_xy_', num2str(endTime), '_', datetimePath, '.mp4');
  v_xyEorB = VideoWriter(xyExyBvideoFileName, 'MPEG-4');
  v_xyEorB.FrameRate = max(minFramerate, floor(ntime/ndskip/videoDuration));
  % v_xyEorB.FrameRate = 30;
  open(v_xyEorB);
end
if check.kxkyPlot
  kxkyEorBvideoFileName = strcat(fileWithoutDotM, '_kxky_', num2str(endTime), '_', datetimePath, '.mp4');
  v_kxkyEorB = VideoWriter(kxkyEorBvideoFileName, 'MPEG-4');
  v_kxkyEorB.FrameRate = max(minFramerate, floor(ntime/ndskip/videoDuration));
  % v_kxkyEorB.FrameRate = 30;
  open(v_kxkyEorB);
end
if check.veloDistPlot
  velocitydistVideoFileName = strcat(fileWithoutDotM, '_velocitydist_', num2str(endTime), '_', datetimePath, '.mp4');
  v_velocitydist = VideoWriter(velocitydistVideoFileName, 'MPEG-4');
  v_velocitydist.FrameRate = max(minFramerate, floor(ntime/ndskip/videoDuration));
  % v_velocitydist.FrameRate = 30;
  open(v_velocitydist);
end
if check.EJ
  EJvideoFileName = strcat(fileWithoutDotM, '_EJ_', num2str(endTime), '_', datetimePath, '.mp4');
  v_EJ = VideoWriter(EJvideoFileName, 'MPEG-4');
  v_EJ.FrameRate = max(minFramerate, floor(ntime/ndskip/videoDuration));
  % v_EJ.FrameRate = 30;
  open(v_EJ);
end
if check.BJ
  BJvideoFileName = strcat(fileWithoutDotM, '_BJ_', num2str(endTime), '_', datetimePath, '.mp4');
  v_BJ = VideoWriter(BJvideoFileName, 'MPEG-4');
  v_BJ.FrameRate = max(minFramerate, floor(ntime/ndskip/videoDuration));
  % v_BJ.FrameRate = 30;
  open(v_BJ);
end