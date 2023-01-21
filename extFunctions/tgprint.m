function tgprint(chat_id, hFig, options)
% TGPRINT send an image to a Telegram bot
%
% Use tgprintf() in the same way as sprintf()
% Example: figure(1); plot(x,y);
%          tgprint();
% 
% There are two sending modes:
% (1) tgprint('photo'): send an image w/ compression using sendPhoto
% (2) tgprint('document'): send an image w/o compression using sendDocument
% 
% Define token and chat_id before use, 
% which are the authorization token of the target Telegram bot 
% and the identifier or username of the target chat
%
% Please refer the following post 
% "Creating a Telegram bot for personal notifications"
% https://www.forsomedefinition.com/automation/creating-telegram-bot-notifications/
% 
% This also uses urlreadpost by Dan Eills
% https://www.mathworks.com/matlabcentral/fileexchange/27189-urlreadpost-url-post-method-with-binary-file-uploading
% 
% Seongsik Park
% seongsikpark@postech.ac.kr

% if nargin > 0
%     options = varargin{1};
% end

filename = 'temp.png';
% if exist(filename,'file')
%     warning('temporay file ''%s'' wil be overwritten!',filename);
% end

% if vargin
% if hFig.
% print(hFig, filename, '-dpng');
% saveas(hFig, filename, 'png');
exportgraphics(hFig, filename)

f = fopen(filename,'rb');
d = fread(f,Inf,'*uint8')';
fclose(f);

% default token and chat_id
token = '1727085776:AAGm7Vfig68zsUR0rxXeP2-dTF4jnQeXY3Q';


if strcmpi(options,'photo')
    sendstr = ['https://api.telegram.org/bot',token,'/sendPhoto'];
    urlreadpost(sendstr,{'chat_id',chat_id,'photo',d}); 
else
    sendstr = ['https://api.telegram.org/bot',token,'/sendDocument'];
    urlreadpost(sendstr,{'chat_id',chat_id,'document',d}); 
end

end