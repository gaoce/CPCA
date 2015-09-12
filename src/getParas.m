function [delta, corr] = getParas()
% Get parameters from user

prompt = {'Proportion of Variance Explained(%)'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'50'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

delta = str2double(answer)/100;
if delta > 1 || delta < 0
    errordlg(['Invalid data: ', num2str(delta)]);
    
    % Set error value and return
    delta = -1;
    corr = -1;
    return
end
    
corr = 1;

end