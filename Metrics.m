function [Metric] = Metrics(obs,est)

%% compute metrics
    % find intersections
        % Correlation
        CORR  = corr(obs,est,'rows', 'complete');
        % beta (Bias)
        beta  = nanmean(est)/nanmean(obs);
        % gamma (Variability)
        % gamma = std(est) / std(obs);
        gamma = (nanstd(est)/nanmean(est))/(nanstd(obs)/nanmean(obs));
        % KGE
        P1 = (CORR - 1)^2;
        P2 = (beta - 1)^2;
        P3 = (gamma -1)^2;
    wr = 1;
    wb = 1;
    wg = 1;
    
    KGE = 1 - sqrt(wr*(CORR-1)^2 + wb*(beta-1)^2 + wg*(gamma-1)^2);


        % KGE = 1 - sqrt(P1+P2+P3); clear P1 P2 P3
        % Mean Bias (MBias)
        MBias = nansum(est)/nansum(obs);
        % Relative Bias (RBias)
        RBias = 100*nansum(est - obs)/nansum(obs);
        % MBE (Mean Bias Error)
        MBE = nansum(est - obs)/length(obs);
        % MAE (Mean Abs Error)
        MAE = nansum(abs(est - obs))/length(obs);
        % RCV (Relative Coefficient of Variations)
        RCV = 100.*(nanstd(est)/nanmean(est))/(nanstd(obs)/nanmean(obs));
        % RMSE (Root Mean Square Error)
        RMSE = sqrt(nanmean((obs - est).^2));
        % NRMSE (Normalized Root Mean Square Error)
        NRMSE = sqrt(nanmean((obs - est).^2))/nanstd(obs);
        % Systematic Error
        MSE = nansum((est - obs).^2)/length(est);
        [p] = polyfit(obs,est,1);
        E_star = p(1)*obs + p(2);
        MSE_s = 100*(nansum((E_star - obs).^2)/length(E_star))/MSE;
        MSE_r = 100*(nansum((est - E_star).^2)/length(E_star))/MSE;
        % HSS
        % % % obs_ = obs;
        % % % est_ = est;
        % % % obs_(obs_<1) = 0;
        % % % est_(est_<1) = 0;
        % % % a = length(find(obs_>0 & est_>0));
        % % % b = length(find(obs_==0 & est_>0));
        % % % c = length(find(obs_>0 & est_==0));
        % % % d = length(find(obs_==0 & est_==0));
        % % % HSS = 2*(a*d - b*c)/[(a+c)*(c+d) + (a+b)*(b+d)];
        % NSE
        NSE = 1 - nansum((obs-est).^2)/nansum((obs-nanmean(obs)).^2);
        % CSNSE

        obs_res=obs-nanmean(obs);
        CSNSE = 1 - nansum((obs-est).^2)/nansum((obs_res).^2);

        Error = est - obs;
        IQR = iqr(Error);

            Q1 = prctile(Error,25);
            Q3 = prctile(Error,75);
            IQR = Q3 - Q1;
            
            idx_iqr = Error >= Q1 & Error <= Q3;
            IQR_RMSE = sqrt(nanmean(Error(idx_iqr).^2));



Metric.Corr  = CORR;
Metric.KGE   = KGE;
Metric.beta  = beta;
Metric.gamma = gamma;
Metric.NSE   = NSE;
Metric.CSNSE = CSNSE;
Metric.RMSE  = RMSE;
Metric.NRMSE = NRMSE;
Metric.MBias = MBias;
Metric.RBias = RBias;
Metric.MBE   = MBE;
Metric.MAE   = MAE;
Metric.RCV   = RCV;
Metric.MSE_s = MSE_s;
Metric.MSE_r = MSE_r;
Metric.IQR  = IQR;
Metric.IQR_RMSE = IQR_RMSE;
% Metric.HSS   = HSS;


end

