classdef QuantileRegression
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x       % (n x 1) data
        y       % (n x 1) data
        p       % (2 x max_noc) parameters
        f       % (2 x max_noc) minima values of optimised function
        q       % (2 x max_noc) changepoints
        B       % number of bootstraps
        default % struct containing defaults
        model   % struct containing model specifications
        BS      % struct containing bootstrapped results
    end % properties
    
    methods
        function QR = QuantileRegression(model, default, smlOn, fitOn)
            %Create QuantileRegression object
            % INPUT
            % model     struct with the following components
            % model.name, model.n, model.tau, model.constraint, model.tau_cp
            % default   struct with the following components
            % default.PltOn, default.saveOn, default.noPrint
            
            %% Check inputs
            if nargin < 1
                return;
            end
            
            if nargin < 3
                smlOn = true;
            end
            
            if nargin < 4
                fitOn = true;
            end
            
            %% Defaults
            QR.model = model;
            QR.default = default;
            
            if ~smlOn
                return
            end
            
            %% Simulate data
            [x, y] = DataSml(model.name, model.n, model.tau_cp, default.PltOn);
            figure(1); clf; plot(x,y,'k.');
            
            QR.x = x;
            QR.y = y;
            
            if ~fitOn
                return
            end
            
            QR.B = model.B;
            
            %% Fit the model
            start_time = tic();
            [p, f, ~, q] = Fit_Model(x, y, model.tau, model.max_n_cp, model.constraint, default.noPrint);
            toc(start_time);
            
            QR.p = p;
            QR.f = f;
            QR.q = q;
            
            %% Bootstrap
            copy_QR = QR;
            
            QR.BS = cell(QR.B,1);
            for i = 1:QR.B
                QR.BS{i} = copy_QR;
            end
            
            bX = cell(QR.B,1); bY = cell(QR.B,1);
            bP = cell(QR.B,1); bF = cell(QR.B,1); bQ = cell(QR.B,1);
            parfor iB = 2:QR.B
                fprintf('Resample %d/%d commenced\n',iB,QR.B);
                [tx, ty] = DataSml(model.name, model.n, model.tau_cp, default.PltOn);
                bX{iB} = tx;
                bY{iB} = ty;
                [bp, bf, ~, bq] = Fit_Model(bX{iB}, bY{iB}, model.tau, model.max_n_cp, model.constraint, default.noPrint);
                bP{iB} = bp; bF{iB} = bf; bQ{iB} = bq;
            end
            
            for iB = 2:QR.B
                QR.BS{iB}.x = bX{iB};
                QR.BS{iB}.y = bY{iB};
                QR.BS{iB}.p = bP{iB};
                QR.BS{iB}.f = bF{iB};
                QR.BS{iB}.q = bQ{iB};
            end            
        end % constructor
        
        function Plot(QR, noc, legendOn)
            %% Plot
            if nargin < 2
                noc = 2;
            end
            
            switch QR.default.purpose
                case {'poster','beamer'}
                    fontsize = 18;
                    legend_flag = false;
                case 'paper'
                    fontsize = 11;
                    legend_flag = true;
            end
            if nargin < 4
                legendOn = legend_flag;
            end
            
            % Legend for plot
            legend_str = cell(length(QR.model.tau)+1,1);
            legend_str{1} = 'data';
            
            FigCnt = 2;
            for j_noc = noc
                figure(FigCnt); FigCnt = FigCnt + 1;
                clf;
                hold on;
                plot(QR.x, QR.y, 'k.')
                
                tx = linspace(min(QR.x),max(QR.x),1000);
                
                colors = jet(length(QR.model.tau));
                for k_quantile = 1:length(QR.model.tau)
                    if isempty(QR.q{1,j_noc})
                        curr_comp = 1; curr_quantile_index = k_quantile;
                    else
                        curr_comp = sum(QR.q{1,j_noc} <= QR.model.tau(k_quantile)) + 1;
                        if curr_comp == 1
                            curr_quantile_index = k_quantile;
                        else
                            curr_quantile_index = sum(QR.model.tau(1:k_quantile) > QR.q{1,j_noc}(curr_comp - 1));
                        end
                    end
                    
                    curr_alpha = QR.p{1, j_noc}{curr_comp}(1);
                    curr_c = QR.p{1, j_noc}{curr_comp}(2);
                    curr_beta = QR.p{1, j_noc}{curr_comp}(3);
                    curr_z = QR.p{1, j_noc}{curr_comp}(3 + curr_quantile_index);
                    
                    plot(tx, curr_c + curr_alpha * tx + curr_z * tx .^ curr_beta, 'LineWidth', 3, 'Color', colors(k_quantile,:))
                    legend_str{k_quantile+1} = sprintf('\\tau = %d%%',round(100*QR.model.tau(k_quantile)));  %sprintf('%.2f + %.2f*x + %.2f*x^{%.2f}', curr_c, curr_alpha, curr_z, curr_beta);
                    
                end
                if legendOn
                    legend(legend_str,'Location','best','fontsize',fontsize/2);
                end
                xlabel('X','fontsize',fontsize);
                ylabel('Y','fontsize',fontsize);
                title(sprintf('Best fit with %d mixture(s)',j_noc),'fontsize',fontsize);
                set(gca,'fontsize',fontsize);
            end
            plot(QR.x, QR.y, 'k.')
        end % Plot
    end % methods
end

