function ret = tgprintf(str)
% TGPRINTF send a message to a Telegram bot
%
% Use tgprintf() in the same way as sprintf()
% Example: tgprintf('Hello, World!');
%          tgprintf('%d + %d = %d',1,2,1+2);
% 
% Define token and chat_id before use, 
% which are the authorization token of the target Telegram bot 
% and the identifier or username of the target chat
%
% default token and chat_id
%
% Bot name: 'LBIS_MATLAB'
%
% Change the chat_id to the ID of your own group with the bot.
%
token = '1727085776:AAGm7Vfig68zsUR0rxXeP2-dTF4jnQeXY3Q';
chat_id = '-512325870';
% str = sprintf(varargin{:});
% print to MATLAB command window
% fprintf(str);
% convert MATLAB string to url query string
sendstr = urlencode(str);
sendstr = ['https://api.telegram.org/bot',token,...
           '/sendMessage?chat_id=',chat_id,...
           '&text=',sendstr];
% send a message   
ret = webread(sendstr); 
% assert(ret.ok);
% append human readable datetime to results [Set TimeZone value to desired time zone]
ret.result.datetime=datetime(ret.result.date,'ConvertFrom','posixtime','TimeZone','Asia/Seoul');
end
