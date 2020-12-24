function h_mainfig = wavbrowser_mainfig()
% policy - create a new figure or use a singleton. 'new' or 'reuse'.

h_mainfig = figure(...
'Units','normalized',...
'Color',[0.941176470588235 0.941176470588235 0.941176470588235],...
'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
'IntegerHandle','off',...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'MenuBar','none',...
'Name','wavbrowser',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'Position',[.01 .05 .98 .9],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','on',...
'HandleVisibility','callback',...
'Tag','mainfig',...
'UserData',[],...
'Behavior',get(0,'defaultfigureBehavior'),...
'CloseRequestFcn',@audioapps_closereq,...
'Visible','on' );

h_title = uicontrol(...
'Parent',h_mainfig,...
'Units','normalized',...
'Position',[.01 .895 .97 .025],...
'String','',...
'Tag','title_bar',...
'Style','text',...
'FontSize',11,...
'FontWeight','normal',...
'HorizontalAlignment','center');

h_daytimeaxis = axes(...
'Parent',h_mainfig,...
'Units','normalized',...
'Position',[.03 .725 .94 .07],...
'FontSize',8,...
'Tag','daytime_axis',...
'XTick',[],...
'YTick',[]);

h_dayaxis = axes(...
'Parent',h_mainfig,...
'Units','normalized',...
'Position',[.03 .825 .94 .07],...
'Tag','day_axis',...
'FontSize',8,...
'XTick',[],...
'YTick',[]);

h_filemenu = uimenu('parent',h_mainfig, 'label','File');
uimenu(h_filemenu,'Label','Load directory','Tag','load_wav_dir','Callback','1');
uimenu(h_filemenu,'Label','Load templates','Tag','load_templates','Callback','1');
uimenu(h_filemenu,'Label','Load matches','Tag','load_matches','Callback','1');

h_controlpanel = uipanel(...
'Parent',h_mainfig,...
'Units','normalized',...
'Position',[.01 .94 .6 .05],...
'Tag','control_panel');


objstrings = {'Prev','Next','Go to','1','Play','Slow fact.','1','Extract','Plot','Zoom mode'};
objtags = {'prev','next','go_to_lab','go_to','play','slow_fact_lab','slow_factor','extract','plot_spec_labs','zoom_mode'};
objstyles = {'pushbutton','pushbutton','text','edit','pushbutton','text','edit','pushbutton','pushbutton','radiobutton'};

put_uicontrol(h_controlpanel,objstrings,objtags,objstyles,'horizontal');


h_feature = axes(...
'Parent',h_mainfig,...
'Units','normalized',...
'Position',[.03 .025 .94 .2],...
'FontSize',8,...
'Tag','feature',...
'XTick',[],...
'YTick',[]);

h_spectrogram = axes(...
'Parent',h_mainfig,...
'Units','normalized',...
'FontSize',8,...
'Position',[.03 .275 .94 .35],...
'Tag','spectrogram',...
'XTick',[],...
'YTick',[]);

h_matchaxis = axes(...
'Parent',h_mainfig,...
'Units','normalized',...
'Position',[.03 .6325 .94 .06],...
'FontSize',8,...
'Tag','match_axis',...
'XTick',[],...
'YTick',[]);

