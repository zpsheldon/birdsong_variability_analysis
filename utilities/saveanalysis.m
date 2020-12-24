function [] = saveanalysis(strct)

[filename,pathname] = uiputfile('*.mat','Choose file name');

if ~filename
    return;
end

fl = [pathname filesep filename];
appendopt = 0;
if exist(fl)
    an = questdlg('Append data?');
    
    switch an 
        case 'Yes'
        appendopt = 1;
        case 'Cancel'
            return;
    end
        
end

varnm = inputdlg('Variable name');
varnm = varnm{1};
eval([varnm '= strct;'])

if appendopt
    save(fl,varnm,'-append')
else
    save(fl,varnm)
end

