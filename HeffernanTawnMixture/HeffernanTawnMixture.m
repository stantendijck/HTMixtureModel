classdef HeffernanTawnMixture
    %HeffernanTawn is a class that allows to model data following the
    %Heffernan-Tawn framework
    
    properties
        x               % (n x (nDmn + 1)) data matrix
        n               % size of the data
        nDmn            % dimension of the model
        MrgTp           % margin type: Gumbel or Laplace (not important)
        TwoParamFlag    % Two parameter flag in fitting model
        p               % exceedance probabilities
        Prm             % MLE of the parameters
        Rsd             % Residuals of the model fit
        cp              % mixture probability/change point
        NLOGL           % Negative log-likelihood of the fit
    end % properties
    
    methods
        function HT = HeffernanTawnMixture(x,p,B,TwoParamFlag)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 2
                p = 0.9;
            end
            if nargin < 3
                B = 10;
            end
            if nargin < 4
                TwoParamFlag = false;
            end
            
            bp = zeros(B,1);
            bp(1) = p;
            bp(2:B) = rand(B-1,1)*(1-p) + (3/2*p-1/2);
            
            Prm = cell(B,1); cp = cell(B,1); %#ok<*PROPLC>
            new_x = cell(B,1); new_n = cell(B,1);
            for iB = 1:B
                new_x{iB} = x(x(:,1)>quantile(x(:,1),bp(iB)),:);
                new_n{iB} = size(new_x{iB},1);
                
                HT.x = new_x{iB};
                HT.n = new_n{iB};
                HT.nDmn = size(HT.x,2)-1;
                HT.TwoParamFlag = TwoParamFlag;
                [prm,NLOGL] = FitHT_M(HT);
                Prm{iB} = prm{1}; cp{iB} = prm{2};
            end
            
            HT.x = new_x{1};
            HT.n = new_n{1};
            HT.p = bp;
            HT.Prm = Prm;
            HT.cp = cp;
            HT.NLOGL = NLOGL;
        end % constructor
        
        function [Prm,Rsd,NLOGL] = FitHT(HT,nBoot,keef)
            % Fit a 2/4 parameter Heffernan and Tawn model conditional on X(:,1) given (no)
            % constraints and with starting values p0
            % INPUTS:
            % X (n x (nDmn+1)),     datamatrix
            % nBoot,                  number of bootstraps
            % TwoParamFlag,           if true, HT model is fitted with two free
            %                         parameters, mu and sigma are estimated via
            %                         theoretical MLE estimates.
            %                         if false, HT model is fitted with 4 free
            %                         parameters
            %
            % OUTPUTS:
            % Prm,                    fitted H&T parameters for iBt'th bootstrap
            % Rsd,                    residuals
            % NLOGL,                  model likelihood
            
            % Defaults
            if nargin < 2
                nBoot = 1;
            end
            if nargin < 3
                keef = true;
            end

            
            opts=optimset('display','off');
            
            % Pre allocation
            Prm = cell(nBoot,1);
            NLOGL = nan(nBoot,1);
            I = (1:HT.n)';
            
            for iBoot = 1:nBoot % loop over bootstraps
                if iBoot == 1
                    if HT.TwoParamFlag
                        p0 = repmat([0;0],1,HT.nDmn); %starting value (alpha, beta)
                    else
                        p0 = repmat([0;0;0;1],1,HT.nDmn); %starting value (alpha, beta, mu, sigma)
                    end
                end
                
                %% Bootstrap the data
                if iBoot ~= 1
                    I = randi(HT.n,HT.n,1);
                end
                bX = cell(HT.n,1);
                for i = 1:HT.n
                    bX{i} = HT.x(I(i),:);
                end
                bX = cell2mat(bX);
                
                tX=bX(:,1);             %exceedances for iBoot'th bootstrap
                tY=bX(:,2:HT.nDmn+1);    %conditioning value given exceedance in conditioned value
                
                %% fit model: find optimal parameter values for HT model via MLE
                [tHTPrm,NLOGL(iBoot)]=fminsearch(@(p)HeffernanTawn.HTlike(tX,tY,p,'Laplace',HT.TwoParamFlag,keef),p0,opts);
                
                if HT.TwoParamFlag
                    if HT.nDmn ~= size(p0,2)
                        alp = tHTPrm(1,:);
                        bet = tHTPrm(2,:);
                        for iHTnDmn = size(p0,2):HT.nDmn
                            alp(iHTnDmn) = alp(1)^iHTnDmn; % works only for HT(1) as null hypothesis
                            bet(iHTnDmn) = bet(1);
                        end
                    else
                        alp = tHTPrm(1,:);
                        bet = tHTPrm(2,1);
                        for iHTnDmn = 2:HT.nDmn
                            bet(iHTnDmn) = bet(1);
                        end
                    end
                    tHTPrm = [alp;bet];
                    % Compute residuals and MLE estimates for mean and std
                    Res = (tY - tX.*tHTPrm(1,:))./(tX.^tHTPrm(2,:));
                    M = mean(Res);
                    S = sqrt(var(Res));
                else
                    M = tHTPrm(3,:); S = tHTPrm(4,:);
                end
                
                Prm{iBoot}=[tHTPrm(1:2,:);M;S];
                
                if iBoot == 1
                    %% keep residuals
                    tAlp = Prm{1}(1,:);
                    tBet = Prm{1}(2,:);
                    tMu = Prm{1}(3,:);
                    tSig = Prm{1}(4,:);
                    
                    tXpB=tX.^tBet;
                    mn=tY-tAlp.*tX-tMu.*tXpB;
                    std=tSig.*tXpB;
                    
                    mn = reshape(mn,[],HT.nDmn);
                    std = reshape(std,[],HT.nDmn);
                    Rsd = nan(HT.n,HT.nDmn);
                    for iDmn = 1:HT.nDmn
                        Rsd(:,iDmn)=mn(:,iDmn)./std(:,iDmn);
                    end
                end
            end %fit
            % fprintf('\n')
        end % FitHT
        
        function [Prm,NLOGL] = FitHT_M(HT)
            
            % First fit HT model
            [Prm,~,~] = FitHT(HT);
            
            % Initialise
            eps = [1e-1,1e-4];
            alpha_0 = max(0,min(1,Prm{1}(1) + [-1,1]*eps(1)));
            alpha_0 = [0.15,0.85];
            beta_0 = [0,0];
            mu_0 = [0,0];
            sigma_0 = [1,1];
            
            curr_alpha = alpha_0;
            curr_beta = beta_0;
            curr_mu = mu_0;
            curr_sigma = sigma_0;
            curr_cp = 0.5;
            
            cnt = 1; curr_diff = 10;
            while (cnt == 1 || curr_diff > eps(2)) && cnt < 100
                %% Allocate the observations
                fprintf('%.2f, %d\n',curr_diff,cnt)
                tx = HT.x(:,1); ty = HT.x(:,2);
                tx_beta = tx.^(curr_beta(cnt,:));
                Rsd1 = (ty - curr_alpha(cnt,1)*tx - curr_mu(cnt,1)*tx_beta(:,1))./(curr_sigma(cnt,1)*tx_beta(:,1));
                Rsd2 = (ty - curr_alpha(cnt,2)*tx - curr_mu(cnt,2)*tx_beta(:,2))./(curr_sigma(cnt,2)*tx_beta(:,2));
                
                L1 = normpdf(Rsd1);
                L2 = normpdf(Rsd2);
                
                fun = @(alp,bet,m,s)(sum(-log(L1./(L1+L2) .* exp(-HeffernanTawnMixture.HTlike(tx,ty,[alp(1);bet(1);m(1);s(1)],'Laplace',HT.TwoParamFlag,true)) + L2./(L1+L2) .* exp(-HeffernanTawnMixture.HTlike(tx,ty,[alp(2);bet(2);m(2);s(2)],'Laplace',HT.TwoParamFlag,true)))));
                [curr_prm,curr_fval] = fminsearch(@(p)(fun(p(1,:),p(2,:),p(3,:),p(4,:))),[curr_alpha(cnt,:);curr_beta(cnt,:);curr_mu(cnt,:);curr_sigma(cnt,:)],optimset('MaxFunEvals',100000));
                
                curr_alpha(cnt+1,:) = curr_prm(1,:);
                curr_beta(cnt+1,:) = curr_prm(2,:);
                curr_mu(cnt+1,:) = curr_prm(3,:);
                curr_sigma(cnt+1,:) = curr_prm(4,:);
                curr_cp(cnt+1) = mean(L1./(L1+L2));
                
                curr_diff = abs(fun(curr_alpha(cnt+1,:),curr_beta(cnt+1,:),curr_mu(cnt+1,:),curr_sigma(cnt+1,:)) - fun(curr_alpha(cnt,:),curr_beta(cnt,:),curr_mu(cnt,:),curr_sigma(cnt,:)));
                
                cnt = cnt + 1;
            end
%             figure(1); clf; plot(curr_alpha)
%             figure(2); clf; plot(curr_beta)
%             figure(3); clf; plot(curr_mu)
%             figure(4); clf; plot(curr_sigma)
%             figure(5); clf; plot(curr_cp)
%             figure(6); clf; plot(HT.x(:,1),HT.x(:,2),'k.'); hold on;
%             XX = min(HT.x(:,1)):0.001:max(HT.x(:,1));
%             plot(XX,XX*curr_prm(1,1) + curr_prm(3,1)*XX.^curr_prm(2,1))
%             plot(XX,XX*curr_prm(1,2) + curr_prm(3,2)*XX.^curr_prm(2,2))
            
            Prm = {curr_prm,curr_cp(end)};
            NLOGL = curr_fval;
        end
    end % methods
    
    methods(Static)
        function NLOGL = HTlike(X,Y,p,MrgTp,TwoParamFlag,keef)
            %Compute Heffernan and Tawn likelihood for data on MrgTp scale.
            %INPUT
            % X,            n x 1 conditioned value
            % Y,            n x nDmn conditioning value
            % p,            1 x (nDmn-1) parameter values
            % MrgTp,        Gumbel or Laplace
            % TwoPrmFlag,   TwoParameters HT Model if true, else 4 parameter
            % keef,         Check keef constraints
            if nargin < 5
                TwoParamFlag = false;
            end
            if nargin < 6
                keef = false;
            end
            
            % Get parameters
            Alp=p(1,:); %nBin x nD-1
            Bet=p(2,:); %1 x nD-1
            if ~TwoParamFlag
                Mu = p(3,:);
                Sig = p(4,:);
                if any(Sig<0)
                    NLOGL=Inf;
                    return
                end
            end
            nDmn = size(Y,2);
            
            % Parameter validation checking
            if any(Alp(:)>1) || any(Alp(:)<-1) || any(Bet>1)
                NLOGL=Inf;
                return
            end
            
            if isequal(MrgTp,'Gumbel')
                if any(Alp(:)<0)
                    NLOGL=Inf;
                    return
                end
            end
            
            kcond = ones(2,1);
            if keef
                kcond(1) = HeffernanTawn.KeefConstraints(X,Y,[Alp;Bet],1);
                kcond(2) = HeffernanTawn.KeefConstraints(X,Y,[Alp;Bet],0);
            end
            
            if ~all(kcond)
                NLOGL = Inf;
                return
            end
            
            
            % Compute the likelihood:
            Xb=bsxfun(@power,X,Bet);  %X^b
            
            Y = reshape(Y,[],nDmn);
            if TwoParamFlag
                Rsd = (Y-Alp.*X)./Xb;
                Mu = mean(Rsd);        %MLE estimate for mean
                Sig = std(Rsd);        %MLE estimate for std
            end
            
            Std=bsxfun(@times,Sig,Xb);  %standard deviation
            NLOGL =0.5*(log(2*pi)+(Y-Alp.*X-bsxfun(@times,Mu,Xb)).^2./(Std).^2)+log(Std);
%             NLOGL=sum(tNLOGL);
        end % HTlike
        
        function kcond = KeefConstraints(X,Y,p,q)
            v = max(X)+1;
            nDmn = size(Y,2);
            
            alpha = p(1,:); beta = p(2,:);
            
            z_plus = Y - X;
            z_min = Y + X;
            step1 = Y - alpha.*X;
            step2 = X.^beta;
            z = step1./step2;
            
            if q == 1
                z_plus_q = max(z_plus);
                z_min_q = max(z_min);
                z_q = max(z);
            elseif q == 0
                z_plus_q = min(z_plus);
                z_min_q = min(z_min);
                z_q = min(z);
            else
                z_plus_q = quantile(z_plus,q);
                z_min_q = quantile(z_min,q);
                z_q = quantile(z,q);
            end
            
            min_vec1 = ones(length(q),3);
            min_vec2 = ones(length(q),3);
            
            for id = 1:nDmn
                min_vec1(:,2) = 1 - beta(id)*v^(beta(id)-1)*z_q(:,id);
                min_vec1(:,3) = 1 - v^(beta(id)-1)*z_q(:,id) + z_plus_q(:,id)/v;
                
                pos_quant1 = (1 - 1/beta(id))*(1 - alpha(id))^(-beta(id)/(1-beta(id)))*(beta(id)*z_q(:,id)).^(1/(1-beta(id))) + z_plus_q(:,id);
                
                cond1_flag1 = alpha(id) <= min(min_vec1,[],2);
                cond1_flag2 = min_vec1(:,2) < alpha(id) & pos_quant1 > 0;
                
                cond1 = cond1_flag1 | cond1_flag2;
                
                min_vec2(:,2) = 1 + beta(id)*v^(beta(id)-1)*z_q(:,id);
                min_vec2(:,3) = 1 + v^(beta(id)-1)*z_q(:,id) - z_min_q(:,id)/v;
                
                pos_quant2 = (1 - 1/beta(id))*(1 + alpha(id))^(-beta(id)/(1-beta(id)))*(-beta(id)*z_q(:,id)).^(1/(1-beta(id))) - z_min_q(:,id);
                
                cond2_flag1 = -alpha(id) <= min(min_vec2,[],2);
                cond2_flag2 = min_vec2(:,2) < -alpha(id) & pos_quant2 > 0;
                
                cond2 = cond2_flag1 | cond2_flag2;
            end
            
            kcond = (cond1 & cond2);
        end % KeefConstraints
        
        function NLOGL = HT_M_like(X,Y,p)
            noc = size(p,2);
            L = zeros(length(X),noc);
            for ic = 1:noc
                L(:,ic) = exp(-HeffernanTawnMixture.HTlike(X,Y,p(:,ic),'Laplace'));
            end
            U = rand(length(X),1);
            A = sum([zeros(size(L,1),1),cumsum(L./sum(L,2),2)] < U,2);
%             A = 0*(A - 1) + 1;
            
            I = cell(noc,1);
            Rsd = zeros(length(X),noc);
            for ic = 1:noc
                I{ic} = A == ic;
                
                ty = Y;
                tx = X;
                tx_b = X.^p(2,ic);
                
                Rsd(:,ic) = (ty - p(1,ic)*tx - p(3,ic)*tx_b)./(p(4,ic)*tx_b);
            end
            
            Rsd_density_estimates = zeros(length(X),noc);
            
            for ic = 1:noc
                for jc = 1:noc
                    if ~isempty(Rsd(I{jc},jc))
                        Rsd_density_estimates(:,ic) = ksdensity(Rsd(I{jc},jc),Rsd(:,ic));%,'BandWidth',0.1607);
                    end
                end
            end
            
            NLOGL = sum(-log(sum((L./sum(L,2)) .* Rsd_density_estimates,2)));
        end % HT_M_like
        
        function kcond = M_KeefConstraints(X, Y, p)
            
            noc = size(p,2);
            v = max(X)+1000;
            
            alpha = p(1,:);
            beta = p(2,:);
            mu = p(3,:);
            sigma = p(4,:);
            
            Rsd(:,1) = (Y - alpha(1)*X - mu(1)*X.^beta(1))./(sigma(1)*X.^beta(1));
            Rsd(:,2) = (Y - alpha(2)*X - mu(2)*X.^beta(2))./(sigma(2)*X.^beta(2));
            
            L = zeros(length(X),2);
            L(:,1) = 1/sqrt(2*pi) * 1./(sigma(1)*X.^beta(1)) .* exp( - Rsd(:,1).^2/2 );
            L(:,2) = 1/sqrt(2*pi) * 1./(sigma(2)*X.^beta(2)) .* exp( - Rsd(:,2).^2/2 );
            
            mixture_probs = mean(L./sum(L,2));
            
            [~,I] = sort(alpha);
            
            alpha = alpha(I);
            beta = beta(I);
            mu = mu(I);
            sigma = sigma(I);
            mixture_probs = mixture_probs(I);
            cs_mixture_probs = [0,cumsum(mixture_probs)]; cs_mixture_probs(end) = 1;
            
            Rsd = Rsd(:,I);
            L = L(:,I);
            
            B = 1;
            A = cell(noc,B);
            tz_q_0 = zeros(noc,B);
            tz_q_1 = zeros(noc,B);
            for b = 1:B
                I = sum([zeros(size(L,1),1),cumsum(L./sum(L,2),2)] < rand(size(L,1),1),2);
                for ic = 1:noc
                    A{ic,b} = I > -1; %I == ic;
                    tz_q_0(ic,b) = min(Rsd(A{ic,b},ic));
                    tz_q_1(ic,b) = max(Rsd(A{ic,b},ic));
                end
            end
            
            z_q_0 = min(tz_q_0,[],2);
            z_q_1 = max(tz_q_1,[],2);
            
            z_q = reshape([z_q_0,z_q_1]',1,[]);
            
            z_min = Y + X;
            t_q = quantile(z_min,cs_mixture_probs(2:end-1));
            z_min_q = [min(z_min),reshape(repmat(t_q,2,1),[],1)',max(z_min)];
            
            z_plus = Y - X;
            t_q = quantile(z_plus,cs_mixture_probs(2:end-1));
            z_plus_q = [min(z_plus),reshape(repmat(t_q,2,1),[],1)',max(z_plus)];
            
            kcond = 1;
            for iq = 0:2*noc-1
                prm_id = floor(iq/2) + 1;
                
                min_vec1 = ones(3,1);
                min_vec2 = ones(3,1);
            
                min_vec1(2) = 1 - beta(prm_id)*v^(beta(prm_id)-1)*z_q(iq+1);
                min_vec1(3) = 1 - v^(beta(prm_id)-1)*z_q(iq+1) + z_plus_q(iq+1)/v;
                
                pos_quant1 = (1 - 1/beta(prm_id))*(1 - alpha(prm_id))^(-beta(prm_id)/(1-beta(prm_id)))*(beta(prm_id)*z_q(iq+1)).^(1/(1-beta(prm_id))) + z_plus_q(iq+1);
                
                cond1_flag1 = alpha(prm_id) <= min(min_vec1);
                cond1_flag2 = min_vec1(2) < alpha(prm_id) & pos_quant1 > 0;
                
                cond1 = cond1_flag1 | cond1_flag2;
                
                min_vec2(2) = 1 + beta(prm_id)*v^(beta(prm_id)-1)*z_q(iq+1);
                min_vec2(3) = 1 + v^(beta(prm_id)-1)*z_q(iq+1) - z_min_q(iq+1)/v;
                
                pos_quant2 = (1 - 1/beta(prm_id))*(1 + alpha(prm_id))^(-beta(prm_id)/(1-beta(prm_id)))*(-beta(prm_id)*z_q(iq+1)).^(1/(1-beta(prm_id))) - z_min_q(iq+1);
                
                cond2_flag1 = -alpha(prm_id) <= min(min_vec2);
                cond2_flag2 = min_vec2(2) < -alpha(prm_id) & pos_quant2 > 0;
                
                cond2 = cond2_flag1 | cond2_flag2;
                kcond = kcond & (cond1 & cond2);
            end
        end % M_KeefConstraints
    end % methods(Static)
end % HeffernanTawn

