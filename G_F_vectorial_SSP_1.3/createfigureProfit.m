function createfigure(zdata1)
%CREATEFIGURE(zdata1)
%  ZDATA1:  surface zdata

%  Auto-generated by MATLAB on 20-Jun-2020 10:55:06

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create surf
surf(zdata1,'Parent',axes1);

% Create zlabel
zlabel({'Profit'});

% Create ylabel
ylabel({'Infrastructure Type'});

% Create xlabel
xlabel({'Service Type'});

view(axes1,[-37.5 30]);
grid(axes1,'on');
saveas(gca,'C:\Users\balio\Dropbox\Apps\ShareLaTeX\computing at the edge --Journal\figures\profits.eps','epsc');

