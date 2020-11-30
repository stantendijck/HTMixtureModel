function updated_boxplot(samples, labels, p)

d = length(samples);
if nargin < 2
    labels = cell(d,1);
    for i = 1:d
        labels{i} = sprintf('%d',i);
    end
end
if nargin < 3
    p = [0.025,0.975];
end

n = zeros(length(samples)+1,1);
for i = 1:d
    n(i+1) = n(i) + length(samples{i});
end

% samples = [-log(1-rand(n(1),1));-log(-log(rand(n(2),1)))];
% labels = [vertcat(repmat({'exp'},n(1),1));vertcat(repmat({'gumb'},n(2),1))];

old_boxplot_format_samples = cell2mat(samples);
old_boxplot_format_labels = [];
for i = 1:d
    vec = repmat({labels{i}},n(i+1)-n(i),1);
    old_boxplot_format_labels = [old_boxplot_format_labels;vertcat(vec(:))];
end


% create axes object and handle
hAx = axes;

% draw a box plot
boxplot(old_boxplot_format_samples, old_boxplot_format_labels)

% get handles to the lines in the HGGroup object
lines = hAx.Children;

% get handles
uw = findobj(lines, 'tag', 'Upper Whisker');
uav = findobj(lines, 'tag', 'Upper Adjacent Value');
lw = findobj(lines, 'tag', 'Lower Whisker');
lav = findobj(lines, 'tag', 'Lower Adjacent Value');
out = findobj(lines, 'tag', 'Outliers');

for i = 1:d
    j = uw(i).XData(1);
    
    samp = samples{j};
    q_lower = quantile(samp,p(1));
    q_upper = quantile(samp,p(2));
    
    uw(i).YData(1,2) = q_upper;
    uav(i).YData(:) = q_upper;
    lw(i).YData(1,1) = q_lower;
    lav(i).YData(:) = q_lower;
    I = samp < q_lower | samp > q_upper;
    out(i).YData = samp(I); out(i).XData = uw(i).XData(1)*ones(size(out(i).YData));
end
        