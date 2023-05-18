%Courant Condition Alert
if app.dxEditField.Value / app.dtEditField.Value <= sqrt(2) * app.cvEditField.Value                
    uialert(app.KEMPO2UIFigure, app.courantCondMsg, 'Warning');
    return
end

%ntime confirmation
if log2(app.ntimeEditField.Value) ~= ceil(log2(app.ntimeEditField.Value)) && app.wkxkyCheckBox.Value
    selection = uiconfirm(app.KEMPO2UIFigure,'We recommend ntime is a power of 2 since we use FFT for dispersion relations',...
        'Confirm ntime',...
            'Icon','warning');
            return
end