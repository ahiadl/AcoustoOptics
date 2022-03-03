function in = dialogBox(ref)
%     d = uifigure('Position',[100 100 366 270]);
%     d = dialog('Position',[300 300 250 150],'Name','My Dialog');
% 
%     txt = uicontrol('Parent',d,...
%                'Style','text',...
%                'Position',[20 80 210 40],...
%                'String', sprintf('Change reference to %.2f', ref) );
% 
%    
%     btn = uicontrol('Parent',d,...
%                'Position',[85 20 70 25],...
%                'String','Close',...
%                'Callback','delete(gcf)');
%            
%     in = uieditfield(d,'numeric',...
%       'Position',[100 175 100 22],...
%       'ValueChangedFcn',@(txt,event) textChanged(txt,lbl));
 
%     btn = uibutton(fig,'push',...
%                'Position',[420, 218, 100, 22],...
%                'ButtonPushedFcn', @(btn,event) plotButtonPushed(btn,ax));
    beep()
    prompt = {sprintf('Change the Gain of the PMT to %d.\n Enter its DC level(mV):', ref)};
    dlgtitle = 'Gain Change';
    definput = {'0'};
    dims = [1 40];
%     opts.Interpreter = 'tex';
    opts = struct();
    answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
    in = str2double(answer);  
%     uiwait(d);
end
% 
% btn = uibutton(fig,'push',...
%                'Position',[420, 218, 100, 22],...
%                'ButtonPushedFcn', @(btn,event) plotButtonPushed(btn,ax));
% end

% % Create the function for the ButtonPushedFcn callback
% function plotButtonPushed(btn,ax)
%         
% end