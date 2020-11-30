classdef HeffernanTawn
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
        NLOGL           % Negative log-likelihood of the fit
    end % properties
    
    methods
        function HT = HeffernanTawn(x,p,B,TwoParamFlag)
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
            
            Prm = cell(B,1); Rsd = cell(B,1); NLOGL = zeros(B,1);
            new_x = cell(B,1); new_n = cell(B,1);
            for iB = 1:B
                new_x{iB} = x(x(:,1)>quantile(x(:,1),bp(iB)),:);
                new_n{iB} = size(new_x{iB},1);
                
                HT.x = new_x{iB};
                HT.n = new_n{iB};
                HT.nDmn = size(HT.x,2)-1;
                HT.TwoParamFlag = TwoParamFlag;
                [prm,r,nl] = FitHT(HT);
                Prm{iB} = prm{1}; Rsd{iB} = r; NLOGL(iB) = nl;
            end
            
            HT.x = new_x{1};
            HT.n = new_n{1};
            HT.p = bp;
            HT.Prm = Prm;
            HT.Rsd = Rsd;
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
            % keef                    if true, keef constraints are taken
            %                         into account, if false then not.
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
        
        function ModelDiagnostics(HT, saveOn, saveNme)
            if nargin < 2
                saveOn = false;
            end
            if nargin < 3
                saveNme = '';
            end
            
            %% Quantile regression
            tau = 0.1:0.05:0.9;
            slopes_regression = zeros(1,length(tau));
            slopes_regression_beta = zeros(1,length(tau));
            
            figure(1); clf;
            colmap = jet(length(tau));
            plot(HT.x(:,1),HT.x(:,2),'k.'); hold on;
            s_x = sort(HT.x(:,1));
            for itau = 1:length(tau)
                pp = HeffernanTawn.quantreg(HT.x(:,1),HT.x(:,2),tau(itau),1,1);
                plot(s_x,polyval(pp,s_x),'-','Color',colmap(itau,:))
                slopes_regression(itau) = pp(1);
            end
            title('Quantile regression');
            xlabel('X');
            ylabel('Y');
            legend('data','1st order quantile fit','location','best');
            if saveOn
                savePics(strcat('figures/quantile_regression_',saveNme,'.pdf'));
            end
            
            figure(15); clf;
            colmap = jet(length(tau));
            plot(HT.x(:,1),HT.x(:,2),'k.'); hold on;
            s_x = sort(HT.x(:,1));
            for itau = 1:length(tau)
                pp = HeffernanTawn.quantreg_beta2(HT.x(:,1),HT.x(:,2),tau(itau),HT.Prm{1}(2,1),1);
                plot(s_x,pp(1)*s_x + pp(2)*s_x.^pp(3),'-','Color',colmap(itau,:))
                slopes_regression_beta(itau) = pp(1);
            end
%             title('Quantile regression');
            xlabel('X');
            ylabel('Y');
            legend('data','q_{\tau}(X)','location','Northwest');
            if saveOn
                savePics(strcat('figures/conditional_quantile_regression_',saveNme,'.pdf'));
            end
            
            tic()
            B = 500;
            tau2 = 0.1:0.1:0.9;
            slopes_test = zeros(B,length(tau2));
            for b = 1:B
                I = randsample(HT.n,HT.n,true);
                new_tx = [HT.x(:,1),HT.x(I,2) - HT.x(I,1) + HT.x(:,1)];
                for itau = 1:length(tau2)
                    pp = HeffernanTawn.quantreg(new_tx(:,1),new_tx(:,2),tau2(itau),1,1);
                    slopes_test(b,itau) = pp(1);
                end
            end
            toc()
            
            
            if 0
                B = 500;
                tau2 = 0.1:0.1:0.9;
                slopes_beta_test = zeros(B,length(tau2));
                for b = 1:B
                    I = randsample(HT.n,HT.n,true);
                    new_tx = [HT.x(:,1),HT.x(I,2) - HT.x(I,1) + HT.x(:,1)];
                    for itau = 1:length(tau2)
                        pp = HeffernanTawn.quantreg(new_tx(:,1),new_tx(:,2),tau2(itau),1,1);
                        slopes_beta_test(b,itau) = pp(1);
                    end
                end
            end
            
            figure(2); clf;
            plot(tau,slopes_regression); hold on;
            plot(tau2,quantile(slopes_test,[0.025,0.975]),'k--')
            title('Slopes of quantile regression plots');
            xlabel('Percentile');
            ylabel('Slope');
            ylim([-0.2 max(max(max(1.2,slopes_regression)),max(max(quantile(slopes_test,[0.025,0.975]))))])
            
            if 0
                figure(2); clf;
                plot(tau,slopes_regression_beta); hold on;
                plot(tau2,quantile(slopes_beta_test,[0.025,0.975]),'k--')
                title('Slopes of quantile regression plots');
                xlabel('Percentile');
                ylabel('Slope');
                %             ylim([-0.2 max(max(max(1.2,slopes_regression_beta)),max(max(quantile(slopes_beta_test,[0.025,0.975]))))])
            end
            if saveOn
                savePics(strcat('figures/quantile_regression_slopes_',saveNme,'.pdf'));
            end
            
            %% Conditional mean and confidence interval according to HT
            q = quantile(HT.Rsd{1}(:,1),[0.025,0.975]);
            
            figure(3); clf;
            plot(HT.x(:,1),HT.x(:,2),'k.');
            hold on;
            s_x = sort(HT.x(:,1));
            plot(s_x,HT.Prm{1}(1,1)*s_x + s_x.^HT.Prm{1}(2,1)*HT.Prm{1}(3,1),'b-')
            plot(s_x,HT.Prm{1}(1,1)*s_x + s_x.^HT.Prm{1}(2,1)*HT.Prm{1}(3,1) + q.*s_x.^HT.Prm{1}(2,1)*HT.Prm{1}(4,1),'r--')
            xlabel('X');
            ylabel('Y');
            title('Heffernan-Tawn model');
            legend('data','Model mean','95% confidence interval');
            if saveOn
                savePics(strcat('figures/conditional_mean_',saveNme,'.pdf'));
            end
            
            %% Conditional residual distribution
            figure(4); clf;
            plot(log(HT.x(:,1)),HT.Rsd{1}(:,1),'k.');
            title('Conditional residual distribution');
            xlabel('X');
            ylabel('Residuals');
            if saveOn
                savePics(strcat('figures/conditional_residuals_',saveNme,'.pdf'));
            end
            
            %% Marginal residual distribution
            figure(5); clf;
            histogram(HT.Rsd{1}(:,1));
            title('Marginal residual distribution');
            xlabel('Residuals');
            ylabel('Counts');
            if saveOn
                savePics(strcat('figures/marginal_residuals_',saveNme,'.pdf'));
            end
            
            %% Estimates of alpha and beta as function of the exceedance-q
            alpha = zeros(length(HT.Prm),HT.nDmn);
            beta = zeros(length(HT.Prm),HT.nDmn);
            for iB = 1:length(HT.Prm)
                alpha(iB,:) = HT.Prm{iB}(1,:);
                beta(iB,:) = HT.Prm{iB}(2,:);
            end
            figure(6); clf;
            subplot(1,2,1);
            plot(HT.p,alpha,'k.'); ylim([0 1])
            ylabel('\alpha');
            xlabel('quantile');
            subplot(1,2,2);
            plot(HT.p,beta,'k.'); ylim([-1 1])
            ylabel('\beta');
            xlabel('quantile');
            if saveOn
                savePics(strcat('figures/parameter_estimates_',saveNme,'.pdf'))
            end
            
            %% Keef constraints
            alpha_vec = -1:0.02:1;
            beta_vec = -1:0.02:1;
            
            f = ones(length(alpha_vec),length(beta_vec));
            
            for q = 0:1
                if q
                    for i_beta = 1:length(beta_vec)
                        if beta_vec(i_beta) > 0
                            for i_alpha = 1:length(alpha_vec)
                                if i_alpha > 1 && f(i_alpha - 1, i_beta) == 0
                                    f(i_alpha, i_beta) = 0;
                                elseif f(i_alpha,i_beta)
                                    f(i_alpha, i_beta) = f(i_alpha, i_beta) * HeffernanTawn.KeefConstraints(HT.x(:,1), HT.x(:,2), [alpha_vec(i_alpha); beta_vec(i_beta)], q);
                                end
                            end
                        else
                            i_alpha = 1; flag = true;
                            while i_alpha <= length(alpha_vec) && flag
                                if i_alpha > 1 && f(i_alpha - 1, i_beta) == 1
                                    flag = false;
                                elseif f(i_alpha,i_beta)
                                    f(i_alpha, i_beta) = f(i_alpha, i_beta) * HeffernanTawn.KeefConstraints(HT.x(:,1), HT.x(:,2), [alpha_vec(i_alpha); beta_vec(i_beta)], q);
                                end
                                i_alpha = i_alpha + 1;
                            end
                        end
                    end
                else
                    for i_alpha = 1:length(alpha_vec)
                        i_beta = 1; flag = true;
                        while i_beta <= sum(beta_vec <= 0) && flag
                            if i_beta > 1 && f(i_alpha, i_beta - 1)
                                flag = false;
                            end
                            f(i_alpha, i_beta) = f(i_alpha, i_beta) * HeffernanTawn.KeefConstraints(HT.x(:,1), HT.x(:,2), [alpha_vec(i_alpha); beta_vec(i_beta)], q);
                            i_beta = i_beta + 1;
                        end
                    end
                end
            end
            
            green_coords = cell(sum(f(:)),1); green_counter = 1;
            red_coords = cell(numel(f) - sum(f(:)),1); red_counter = 1;
            for i_alpha = 1:length(alpha_vec)
                for i_beta = 1:length(beta_vec)
                    if f(i_alpha,i_beta)
                        green_coords{green_counter} = [alpha_vec(i_alpha),beta_vec(i_beta)];
                        green_counter = green_counter+1;
                    else
                        red_coords{red_counter} = [alpha_vec(i_alpha),beta_vec(i_beta)];
                        red_counter = red_counter+1;
                    end
                end
            end
            
            green_coords_mat = cell2mat(green_coords);
            red_coords_mat = cell2mat(red_coords);
            
            figure(16); clf;
            plot(green_coords_mat(:,1),green_coords_mat(:,2),'g.'); hold on;
            plot(red_coords_mat(:,1),red_coords_mat(:,2),'r.');
            plot(alpha,beta,'k*');
            [non_keef_prm,~,~] = FitHT(HT,1,false);
            plot(non_keef_prm{1}(1,1),non_keef_prm{1}(2,1),'ko');
            xlabel('\alpha');
            ylabel('\beta');
            if saveOn
                savePics(strcat('figures/keef_constraints_',saveNme,'.pdf'))
            end
            
        end % ModelDiagnostics
        
        function Plot(HT)
            figure(1); clf;
            for i = 1:HT.nDmn
                subplot(1,HT.nDmn,i);
                plot(HT.x(:,1),HT.Rsd{1}(:,i),'k.');
            end
            
            figure(2); clf;
            for i = 1:HT.nDmn
                subplot(1,HT.nDmn,i);
                plot(HT.x(:,1),HT.x(:,i+1),'k.');
                hold on;
                xx = xlim;
                xl = linspace(min(HT.x(:,1)),xx(2),1000);
                plot(xl,HT.Prm{1}(1,i)*xl + xl.^HT.Prm{1}(2,i).*HT.Prm{1}(3,i),'r-');
            end
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
            tNLOGL =sum(0.5*(log(2*pi)+(Y-Alp.*X-bsxfun(@times,Mu,Xb)).^2./(Std).^2)+log(Std));
            NLOGL=sum(tNLOGL);
        end % HTlike
        
        function kcond = KeefConstraints(X,Y,p,q,X_beta)
            if size(q,1) < size(q,2)
                q = q';
            end
            v = max(X)+1;
            nDmn = size(Y,2);
            
            alpha = p(1,:); beta = p(2,:);
            
            if any(abs(alpha) > 1) || any(beta > 1)
                kcond = zeros(size(q));
            end
            
            z_plus = Y - X;
            z_min = Y + X;
            if nargin < 5
                X_beta = X.^beta;
            end
            z = (Y - alpha.*X)./X_beta;
            
            if q == 1
                z_plus_q = max(z_plus);
                z_min_q = max(z_min);
                z_q = max(z);
            elseif q == 0
                z_plus_q = min(z_plus);
                z_min_q = min(z_min);
                z_q = min(z);
            elseif all(q == [0;1])
                z_plus_q = [min(z_plus);max(z_plus)];
                z_min_q = [min(z_min);max(z_min)];
                z_q = [min(z);max(z)];
            else
                z_plus_q = quantile(z_plus,q);
                z_min_q = quantile(z_min,q);
                z_q = quantile(z,q);
            end
            
            min_vec1 = ones(length(q),3);
            min_vec2 = ones(length(q),3);
            
            kcond = 1;
            for id = 1:nDmn
                min_vec1(:,2) = 1 - beta(id)*v^(beta(id)-1)*z_q(:,id);
                min_vec1(:,3) = 1 - v^(beta(id)-1)*z_q(:,id) + z_plus_q(:,id)/v;
                
                a_beta = (1 - alpha(id))^(-beta(id)/(1-beta(id)));
                bz_beta = exp(log(beta(id)*z_q(:,id))/(1-beta(id)));
                pos_quant1 = (1 - 1/beta(id)) * a_beta * bz_beta + z_plus_q(:,id);
                
                cond1_flag1 = alpha(id) <= min(min_vec1,[],2);
                cond1_flag2 = min_vec1(:,2) < alpha(id) & pos_quant1 > 0;
                
                cond1 = cond1_flag1 | cond1_flag2;
                
                min_vec2(:,2) = 1 + beta(id)*v^(beta(id)-1)*z_q(:,id);
                min_vec2(:,3) = 1 + v^(beta(id)-1)*z_q(:,id) - z_min_q(:,id)/v;
                
                a_beta = (1 + alpha(id))^(-beta(id)/(1-beta(id)));
                bz_beta = exp(log(-beta(id)*z_q(:,id))/(1-beta(id)));
                pos_quant2 = (1 - 1/beta(id)) * a_beta * bz_beta - z_min_q(:,id);
                
                cond2_flag1 = -alpha(id) <= min(min_vec2,[],2);
                cond2_flag2 = min_vec2(:,2) < -alpha(id) & pos_quant2 > 0;
                
                cond2 = cond2_flag1 | cond2_flag2;
                
                kcond = kcond & (cond1 & cond2);
            end
            
            
        end % KeefConstraints
    end % methods(Static)
end % HeffernanTawn

