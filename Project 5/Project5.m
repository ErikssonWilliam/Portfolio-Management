
cov = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Covariance Matrix', 'Range', 'BJ3:BS12');
short_bond_cf = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Calculations', 'Range', 'F9:F10');
short_bond_date = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Calculations', 'Range', 'E9:E10');
long_bond_date = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Calculations', 'Range', 'J9:J19');
long_bond_cf = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Calculations', 'Range', 'K9:K10');
portfolio_dates = readmatrix('portfolioManagerV4Projekt5.xls', 'Sheet', 'Answer', 'Range', 'A2:A5');

[eigenvectors, eigenvalues] = eig(cov);
eigenvalues=diag(eigenvalues);

time_to_maturity_short = zeros(length(portfolio_dates), length(short_bond_date));
time_to_maturity_long = zeros(length(portfolio_dates), length(long_bond_date));

%Calculate time to each cash flow
for i = 1:length(portfolio_dates)
    for j = 1:length(long_bond_date)
        time_to_maturity_long(i, j) = yearfrac(portfolio_dates(i), long_bond_date(j));

        if j <= length(short_bond_date)
            time_to_maturity_short(i, j) = yearfrac(portfolio_dates(i), short_bond_date(j));
        end

    end
end

%Here we shall interpolate the eigenvectors%

function yearFraction = yearfrac(startDate, endDate)
    yearFraction = daysact(startDate, endDate) / 365.25; % Using actual day count
end
