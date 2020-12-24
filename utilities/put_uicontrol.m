function varargout = put_uicontrol(h_parent,objstrings,objtags,objstyles,orientation)

if nargin < 4
    objstyles = repmat({'pushbutton'},1,length(objstrings));
end

if nargin < 5
    orientation = 'vertical';
end


for i = 1:length(objstrings)

    cb = [];
    switch orientation

        case 'vertical'

            switch objstyles{i}
                case 'text'
                    pos = [.05 1-i/(length(objstrings)+.5) .9 .7/(length(objstrings)+.5)];
                case 'edit'
                    pos = [.2 1-i/(length(objstrings)+.5)+.2/(length(objstrings)+.5) .6 .7/(length(objstrings)+.5)];
                otherwise
                    pos = [.05 1-i/(length(objstrings)+.5) .9 .9/(length(objstrings)+.5)];
                    cb = '1';
            end

        case 'horizontal'

            switch objstyles{i}
                case 'text'
                    pos = [(i-1)/(length(objstrings)+.2)+.01 .1 .7/(length(objstrings)+.2) .6];
                case 'edit'
                    pos = [(i-1)/(length(objstrings)+.2)+.01 .1 .7/(length(objstrings)+.2) .6];
                otherwise
                    pos = [(i-1)/(length(objstrings)+.2)+.01 .1 .7/(length(objstrings)+.2) .8];
                    cb = '1';
            end

    end
    
    htmp = uicontrol(...
        'Parent',h_parent,...
        'Units','normalized',...
        'Callback',cb,...
        'Position',pos,...
        'String',objstrings{i},...
        'Style',objstyles{i},...
        'Tag',objtags{i});
end