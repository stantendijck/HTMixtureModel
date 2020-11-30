function PltHT(HT, noc)
%plotHT model
if nargin < 2
    noc = 2;
end

for i = noc
    Alp = zeros(length(HT{i}.p),i);
    Bet = zeros(length(HT{i}.p),i);
    Mu = zeros(length(HT{i}.p),i);
    Sig = zeros(length(HT{i}.p),i);
    for j = 1:length(HT{i}.p)
        Alp(j,:) = HT{i}.p{j}(1,:);
        Bet(j,:) = HT{i}.p{j}(2,:);
        Mu(j,:) = HT{i}.p{j}(3,:);
        Sig(j,:) = HT{i}.p{j}(4,:);
    end
    
    figure(i); clf
    subplot(2,2,1);
    plot(Alp);
    ylabel('\alpha')
    subplot(2,2,2);
    plot(Bet);
    ylabel('\beta');
    subplot(2,2,3);
    plot(Mu);
    ylabel('\mu');
    subplot(2,2,4);
    plot(Sig);
    ylabel('\sigma');
end
end