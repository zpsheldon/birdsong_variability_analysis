function audioapps_closereq(src,evnt)
% User-defined close request function
% to display a question dialog box
selection = questdlg('Close application?',...
    'Close Request Function',...
    'Yes','No','Yes');
switch selection,
    case 'Yes',
        delete(gcf)
    case 'No'
        return
end
