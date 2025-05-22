
short_bond_cf = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Calculations', 'Range', 'F9:F10');
short_bond_date = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Calculations', 'Range', 'E9:E10');
long_bond_date = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Calculations', 'Range', 'J9:J18');
long_bond_cf = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Calculations', 'Range', 'K9:K18');
portfolio_dates = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Answer', 'Range', 'A2:A6');
swap_rates = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Timeseries', 'Range', 'AO12:AX1055');

swap_rates = swap_rates/100;
[dates, maturities] = size(swap_rates);
r_eff = zeros(dates,maturities);
r_cont = zeros(dates,maturities);
forward_rates = zeros(dates,maturities);
%Bootstrap

for t = 1:dates
    for T = 1:maturities
        if T ==1
        r_eff(t,T) = (1 + swap_rates(t,T)^(1/T)-1);

        else
        sum_term = 0;
            for i = 1:(T-1)
                sum_term = sum_term + swap_rates(t,T) / (1+r_eff(t,i))^i;
            end
            r_eff  (t,T) = (((1+ swap_rates(t,T))/(1-sum_term))^(1/T)-1);
        end
    end 
end

for i=1:dates
    for j= 1:maturities
        r_cont(i,j) = log(1+r_eff(i,j));
    end
end

for t = 1:dates
    forward_rates(t,1) = r_cont(t,1);
    for T1 = 1:(maturities-1)
        T2 = T1 + 1;
        forward_rates(t,T2) = (r_cont(t,T2) * T2 - r_cont(t,T1) * T1) / (T2 - T1);
    end
end


diff_r = diff(r_cont);
diff_f = diff(forward_rates);

cov_zero_coupon = cov(diff_r);
cov_forward = cov(diff_f);

[eigenvectors_zero_coupon, eigenvalues_zero_coupon] = eig(cov_zero_coupon);
eigenvalues_zero_coupon=diag(eigenvalues_zero_coupon);
[eigenvalues_zero_coupon_sorted, i] = sort(eigenvalues_zero_coupon, 'descend');
eigenvectors_zero_coupon_sorted = eigenvectors_zero_coupon(:, i);
% 


[eigenvectors_forward, eigenvalues_forward] = eig(cov_forward);
eigenvalues_forward=diag(eigenvalues_forward);
[eigenvalues_forward_sorted, j] = sort(eigenvalues_forward, 'descend');
eigenvectors_forward_sorted = eigenvectors_forward(:, j);


Number_of_PC = 3;
PC_r = zeros(maturities,Number_of_PC);
PC_f = zeros(maturities,Number_of_PC);

for i = 1:Number_of_PC
    PC_r(:,i) = eigenvectors_zero_coupon_sorted(:,i);
    PC_f(:,i) = eigenvectors_forward_sorted(:,i);
end

writematrix(PC_r,'portfolioManagerV4Projekt5.xls', 'Sheet', 'Answer', 'Range', 'C9:E18')
writematrix(PC_f,'portfolioManagerV4Projekt5.xls', 'Sheet', 'Answer', 'Range', 'F9:H18')
portfolio_time = length(portfolio_dates);

time_to_maturity_short = zeros(portfolio_time, length(short_bond_date));
time_to_maturity_long = zeros(portfolio_time, length(long_bond_date));

interpolated_eigenvectors_short = zeros(length(short_bond_date),Number_of_PC);
interpolated_eigenvectors_long = zeros(length(long_bond_date),Number_of_PC);

interpolated_zero_coupon_short = zeros(portfolio_time, length(short_bond_date));
interpolated_zero_coupon_long = zeros(portfolio_time, length(long_bond_date));

%Calculate time to each cash flow
for i = 1:portfolio_time
    for j = 1:length(long_bond_date)

    time_to_maturity_long(i, j) = yearfrac(portfolio_dates(i), long_bond_date(j));
    maturity = time_to_maturity_long(1, j); % Use the first portfolio date's maturities
    year1 = floor(maturity);
    year2 = ceil(maturity);
    year1 = max(1, min(year1, size(PC_r, 1)));
    year2 = max(1, min(year2, size(PC_r, 1)));

    if year1 == year2
        interpolated_eigenvectors_long(j, :) = PC_r(year1, :);
        interpolated_zero_coupon_long(i, j) = r_eff(i, j);
    else
        fraction = maturity - year1;
        interpolated_eigenvectors_long(j, :) = (1 - fraction) * PC_r(year1, :) + fraction * PC_f(year2, :);
        interpolated_zero_coupon_long(i, j) = (year1 * r_eff(i, year1) + fraction* forward_rates(i,year2))/(year2);
    end

        if j <= length(short_bond_date)
            time_to_maturity_short(i, j) = yearfrac(portfolio_dates(i), short_bond_date(j));
            
           maturity = time_to_maturity_short(1, j); % Use the first portfolio date's maturities
            year1 = floor(maturity);
            year2 = ceil(maturity);
            year1 = max(1, min(year1, size(PC_r, 1)));
            year2 = max(1, min(year2, size(PC_r, 1)));
        
            if year1 == year2
                interpolated_eigenvectors_short(j, :) = PC_r(year1, :);
                interpolated_zero_coupon_short(i, j) = r_eff(i, j);
            else
                fraction = maturity - year1;
                interpolated_eigenvectors_short(j, :) = (1 - fraction) * PC_r(year1, :) + fraction * PC_f(year2, :);
                interpolated_zero_coupon_short(i, j) = (year1 * r_eff(i, year1) + fraction* forward_rates(i,year2))/(year2);
            end
                
        end

    end
end

exposure_cf_long = zeros(length(portfolio_dates)-1,Number_of_PC); % Portfolio exposure to each PC
exposure_cf_short = zeros(length(portfolio_dates)-1,Number_of_PC);
total_exposure = zeros(length(portfolio_dates)-1,Number_of_PC);

for i = 1:length(portfolio_dates)-1
    for j = 1:length(long_bond_date)
       %Equation 1: Derivative of discount factors 
        r_t = interpolated_zero_coupon_long(i+1, :);
        T_long = time_to_maturity_long(i+1, :);
        d_t_long = exp(-diag(r_t) * T_long'); % Discount factors for long bond
        disc_factor_long = -diag(T_long) * diag(d_t_long) * interpolated_eigenvectors_long;

        % Equation 2: Derivative of discount factors 
        r_t = interpolated_zero_coupon_long(i, :);
        d_t_long = exp(-diag(r_t) * T_long'); % Discount factors for long bond
        gradient_disc_factor_long = -diag(T_long) * diag(d_t_long) * interpolated_eigenvectors_long;

        %Equation 3: Exposure of each cash flow 
        exposure_cf_long(i,:) = gradient_disc_factor_long' * long_bond_cf; % Element-wise multiplication

        if j <= length(short_bond_date)
        %Equation 1: Derivative of discount factors 
        r_t_short = interpolated_zero_coupon_short(i+1, :);
        T_short = time_to_maturity_short(i+1, :);
        d_t_short = exp(-diag(r_t_short) * T_short');
        disc_factor_short = -diag(T_short) * diag(d_t_short) * interpolated_eigenvectors_short;

        % Equation 2: Derivative of discount factors 
        r_t_short = interpolated_zero_coupon_short(i, :);
        d_t_short = exp(-diag(r_t_short) * T_short');
        gradient_disc_factor_short = -diag(T_short) * diag(d_t_short) * interpolated_eigenvectors_short;

        % Equation 3: Exposure of each cash flow 
        exposure_cf_short(i,:) = gradient_disc_factor_short' * short_bond_cf; % Element-wise multiplication

        end
        % Equation 4: Portfolio exposure
        total_exposure(i,:) = exposure_cf_short(i,:) + exposure_cf_long(i,:);
    end
   
end

   writematrix(total_exposure,'portfolioManagerV4Projekt5.xls', 'Sheet', 'Answer', 'Range', 'C2:E5')

    %Hedging
  OIS_cf = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Calculations', 'Range', 'U8:U9');
  OIS_r = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Calculations', 'Range', 'T8:T9');
  OIS_maturity = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Calculations', 'Range', 'R8:R9');
  OIS_years = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Calculations', 'Range', 'S8:S9');
  %no_bonds = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Calculations', 'Range', 'C8:C8');
  number_of_OIS = length(OIS_maturity);
  
  interpolated_eigenvectors_OIS = zeros(number_of_OIS, Number_of_PC);
    for j = 1:number_of_OIS
        maturity = OIS_years(j);
        year1 = floor(maturity);
        year2 = ceil(maturity);
        year1 = max(1, min(year1, size(PC_r, 1)));
        year2 = max(1, min(year2, size(PC_r, 1)));
    
        if year1 == year2
            interpolated_eigenvectors_OIS(j, :) = PC_r(year1, :);
        else
            interpolated_eigenvectors_OIS(j, :) = (PC_r(year1, :) + PC_r(year2, :)) / 2;
        end
    end

  OIS_exposure = zeros(length(portfolio_dates)-2, 3, number_of_OIS); % Exposure to PC1 and PC3

for i = 2:length(portfolio_dates)-1
    for j = 1:number_of_OIS
       %Equation 1: Derivative of discount factors 
        r_t =  OIS_r(j);
        T_long =  daysact(portfolio_dates(i),OIS_maturity(j))/365;
        d_t_long = exp(-diag(r_t) * T_long'); % Discount factors for long bond
        disc_factor_long = -diag(T_long) * diag(d_t_long) * interpolated_eigenvectors_OIS(j,:);

        % Equation 2: Derivative of discount factors 
        r_t =  OIS_r(j);
        d_t_long = exp(-diag(r_t) * T_long'); % Discount factors for long bond
        gradient_disc_factor_long = -diag(T_long) * diag(d_t_long) *  interpolated_eigenvectors_OIS(j,:);

        %Equation 3: Exposure of each cash flow 
        exposure_cf_long(i,:) = gradient_disc_factor_long' * OIS_cf(j); % Element-wise multiplication
        % Equation 4: Portfolio exposure
        OIS_exposure(i,:,j) = exposure_cf_long(i,:);
    end
    if i == 2
    %A_hedge_t = OIS_exposure(i, [1, 3],1) + OIS_exposure(i, [1, 3],2); % OIS exposures to PC1 and PC3 on day t
    A_hedge_t = squeeze(OIS_exposure(i, [1, 3],:));
    
    v_hedge_t = total_exposure(i, [1, 3])'; % Portfolio exposure on day t 
    hedge_positions = A_hedge_t \ -v_hedge_t;
    end
end

total_exposure_hedge = OIS_exposure(:,:,1) * hedge_positions(1,1) + OIS_exposure(:,:,2) * hedge_positions(2,1);
writematrix(total_exposure_hedge,'portfolioManagerV4Projekt5.xls', 'Sheet', 'Answer', 'Range', 'F2:H5')

function yearFraction = yearfrac(startDate, endDate)
    yearFraction = daysact(startDate, endDate) / 365.25; % Using actual day count
end
