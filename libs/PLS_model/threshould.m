function plsda_thres=threshould(yprev_cal,ycal)

rsize = 100;
for g=1:size(yprev_cal,2)
    class_in = ones(size(ycal,1),1);
    class_in(find(ycal ~= g)) = 2;
    count = 0;
    y_in = yprev_cal(:,g);
    miny = min(y_in);
    thr = max(y_in);
    step = (thr - miny)/rsize;
    spsn = [];
    while thr > miny
        count = count + 1;
        class_calc_in = ones(size(ycal,1),1);
        thr = thr - step;
        sample_out_g = find(y_in < thr);
        class_calc_in(sample_out_g) = 2;
        cp = calc_class_param(class_calc_in,class_in);
        sp(count,g) = cp.specificity(1);
        sn(count,g) = cp.sensitivity(1);
        thr_val(count,g) = thr;       
    end
end

% find best thr based on bayesian discrimination threshold
for g=1:max(ycal)
    P_g = yprev_cal(find(ycal==g),g);
    P_notg = yprev_cal(find(ycal~=g),g);
    m_g = mean(P_g); s_g = std(P_g);
    m_notg = mean(P_notg); s_notg = std(P_notg);
    stp = abs(m_g - m_notg)/1000;
    where = [m_notg:stp:m_g];
    % fit normal distribution
    % npdf_g = normpdf(where,m_g,s_g);
    x_g = (where - m_g) ./ s_g;
    npdf_g = exp(-0.5 * x_g .^2) ./ (sqrt(2*pi) .* s_g);
    %npdf_notg = normpdf(where,m_notg,s_notg);
    x_notg = (where - m_notg) ./ s_notg;
    npdf_notg = exp(-0.5 * x_notg .^2) ./ (sqrt(2*pi) .* s_notg);
    minval = NaN;
    for k=1:length(where)
        diff = abs(npdf_g(k)-npdf_notg(k));
        if isnan(minval)|diff < minval
             minval = diff;
             class_thr(g) = where(k);
        end
    end
    if isnan(minval)
        class_thr(g) = mean([m_g m_notg]);
    end  
end

res.class_thr = class_thr;
res.sp = sp;
res.sn = sn;
res.thr_val = thr_val;

plsda_thres=res;


