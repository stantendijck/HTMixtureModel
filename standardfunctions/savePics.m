function savePics(Nme, saveOn, purpose, width, height)

if nargin ~= 1 && ~saveOn
    return
end

if nargin < 3
    purpose = 'paper';
end

if nargin < 4
    switch purpose
        case 'paper'
            width = 10/2;     % Width in inches
            height = 10*2/3/2;    % Height in inches
        case 'beamer'
            width = 10;
            height = 10*2/3;
        case 'poster'
            width = 10;
            height = 10*2/3;
    end
end
papersize = [width,height];

if ~strcmp(get(gcf,'WindowStyle'),'docked')
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
end

set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
set(gcf, 'PaperSize',papersize);
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
fontname='Helvetica';
set(findall(gcf,'type','text'),'fontname',fontname)
set(findall(gcf,'Type','axes'),'fontname',fontname)
% print(Nme,'-djpeg','-r150');
print(Nme,'-dpdf','-r300');
% print(Nme,'-dpdf','-r200');
%print(Nme,'-depsc2','-r600');