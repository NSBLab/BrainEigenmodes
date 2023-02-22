% This script help produce a new 'jet'-like colormap based on other RGB reference colors

% ------- I WAS ASKED ---------------
% "is there a chance that you could add a diverging map going from blue to green to red as in jet, 
% but using the red and blue from your RdBu map and the third darkest green from your RdYlGn map?""
%
% ANSWER:
%   You should construct the new colormap based on the existing RGB values of 'jet'
%   but projecting these RGB values on your new RGB basis. 
% -----------------------------------

% load colormaps
jet=colormap('jet');
RdBu=cbrewer('div', 'RdBu', 11);
RdYlGn=cbrewer('div', 'RdYlGn', 11);

% Define the new R, G, B references (p stands for prime)
Rp=RdBu(1,:);
Bp=RdBu(end, :);
Gp=RdYlGn(end-2, :);
RGBp=[Rp;Gp;Bp];

% construct the new colormap based on the existing RGB values of jet
% Project the RGB values on your new basis
newjet = jet*RGBp;

% store data in a strcuture, easier to handle
cmap.jet=jet;
cmap.newjet=newjet;
cnames={'jet', 'newjet'};

% plot the RGB values
fh=figure();
colors={'r', 'g', 'b'};
for iname=1:length(cnames)
   subplot(length(cnames),1,iname)
   dat=cmap.(cnames{end-iname+1});
   for icol=1:size(dat,2)
      plot(dat(:,icol), 'color', colors{icol}, 'linewidth', 2);hold on;
   end % icol
   title([' "' cnames{end-iname+1} '" in RGB plot'])
end

% plot the colormaps
fh=figure();
for iname=1:length(cnames)
   F=cmap.(cnames{iname});
   ncol=length(F);
   fg=1./ncol; % geometrical factor
   X=fg.*[0 0 1 1];
   Y=0.1.*[1 0 0 1]+(2*iname-1)*0.1;

   for icol=1:ncol
      X2=X+fg.*(icol-1);
      fill(X2,Y,F(icol, :), 'linestyle', 'none')
      hold all
   end % icol
   text(-0.1, mean(Y), cnames{iname}, 'HorizontalAlignment', 'right', 'FontWeight', 'bold', 'FontSize', 10, 'FontName' , 'AvantGarde')
   xlim([-0.4, 1])
   axis off
   set(gcf, 'color', [1 1 1])
   ylim([0.1 1.05.*max(Y)]);
 end % iname

