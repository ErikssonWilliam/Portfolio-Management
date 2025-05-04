
cov = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Covariance Matrix', 'Range', 'BJ3:BS12');
short_bond_cf = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Calculations', 'Range', 'F9:F10');
short_bond_date = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Calculations', 'Range', 'E9:E10');
long_bond_date = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Calculations', 'Range', 'J9:J18');
long_bond_cf = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Calculations', 'Range', 'K9:K18');
portfolio_dates = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Answer', 'Range', 'A2:A5');
zero_coupon_rates = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Rates', 'Range', 'C7:L10');
forward_rates = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Rates', 'Range', 'O7:X10');

[eigenvectors, eigenvalues] = eig(cov);
eigenvalues=diag(eigenvalues);

time_to_maturity_short = zeros(length(portfolio_dates), length(short_bond_date));
time_to_maturity_long = zeros(length(portfolio_dates), length(long_bond_date));

interpolated_eigenvectors_short = zeros(length(short_bond_date),length(eigenvalues));
interpolated_eigenvectors_long = zeros(length(long_bond_date),length(eigenvalues));

interpolated_zero_coupon_short = zeros(length(portfolio_dates), length(short_bond_date));
interpolated_zero_coupon_long = zeros(length(portfolio_dates), length(long_bond_date));

%Calculate time to each cash flow
for i = 1:length(portfolio_dates)
    for j = 1:length(long_bond_date)

    time_to_maturity_long(i, j) = yearfrac(portfolio_dates(i), long_bond_date(j));
    maturity = time_to_maturity_long(1, j); % Use the first portfolio date's maturities
    year1 = floor(maturity);
    year2 = ceil(maturity);
    year1 = max(1, min(year1, size(eigenvectors, 1)));
    year2 = max(1, min(year2, size(eigenvectors, 1)));

    if year1 == year2
        interpolated_eigenvectors_long(j, :) = eigenvectors(year1, :);
        interpolated_zero_coupon_long(i, j) = zero_coupon_rates(i, j);
    else
        fraction = maturity - year1;
        interpolated_eigenvectors_long(j, :) = (1 - fraction) * eigenvectors(year1, :) + fraction * eigenvectors(year2, :);
        interpolated_zero_coupon_long(i, j) = (year1 * zero_coupon_rates(i, year1) + fraction* forward_rates(i,year2))/(year2);
    end

        if j <= length(short_bond_date)
            time_to_maturity_short(i, j) = yearfrac(portfolio_dates(i), short_bond_date(j));
            
           maturity = time_to_maturity_short(1, j); % Use the first portfolio date's maturities
            year1 = floor(maturity);
            year2 = ceil(maturity);
            year1 = max(1, min(year1, size(eigenvectors, 1)));
            year2 = max(1, min(year2, size(eigenvectors, 1)));
        
            if year1 == year2
                interpolated_eigenvectors_short(j, :) = eigenvectors(year1, :);
                interpolated_zero_coupon_short(i, j) = zero_coupon_rates(i, j);
            else
                fraction = maturity - year1;
                interpolated_eigenvectors_short(j, :) = (1 - fraction) * eigenvectors(year1, :) + fraction * eigenvectors(year2, :);
                interpolated_zero_coupon_short(i, j) = (year1 * zero_coupon_rates(i, year1) + fraction* forward_rates(i,year2))/(year2);
            end
                
        end

    end
end

exposure_long = zeros(length(long_bond_date),length(eigenvalues)); % Portfolio exposure to each PC
exposure_short = zeros(length(short_bond_date),length(eigenvalues));

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
        exposure_cf_long = gradient_disc_factor_long' * long_bond_cf; % Element-wise multiplication

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
        exposure_cf_short = gradient_disc_factor_short' * short_bond_cf; % Element-wise multiplication

        end
        % Equation 4: Portfolio exposure
        total_exposure = exposure_cf_short + exposure_cf_long;
    end
   
end

%Here we shall interpolate the eigenvectors%

function yearFraction = yearfrac(startDate, endDate)
    yearFraction = daysact(startDate, endDate) / 365.25; % Using actual day count
end
