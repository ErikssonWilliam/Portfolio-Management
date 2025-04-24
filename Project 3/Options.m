
t = readmatrix('portfolioManagerV4Projekt3.xlsm', 'Sheet', 'Refinitiv', 'Range', 'T4:T9');
S = readmatrix('portfolioManagerV4Projekt3.xlsm', 'Sheet', 'Refinitiv', 'Range', 'U4:U9');
q = readmatrix('portfolioManagerV4Projekt3.xlsm', 'Sheet', 'Refinitiv', 'Range', 'V4:V9');
r = readmatrix('portfolioManagerV4Projekt3.xlsm', 'Sheet', 'Refinitiv', 'Range', 'W4:W9');
K = readmatrix('portfolioManagerV4Projekt3.xlsm', 'Sheet', 'Refinitiv', 'Range', 'X4:X9');
Maturity = readmatrix('portfolioManagerV4Projekt3.xlsm', 'Sheet', 'Refinitiv', 'Range', 'W1:W1');
BS_Call = readmatrix('portfolioManagerV4Projekt3.xlsm', 'Sheet', 'Refinitiv', 'Range', 'Y4:Y9');
BS_Put = readmatrix('portfolioManagerV4Projekt3.xlsm', 'Sheet', 'Refinitiv', 'Range', 'Z4:Z9');

red_days = [
    '2025-04-17'; % Maundy Thursday (Skærtorsdag)
    '2025-04-18'; % Good Friday (Langfredag)
    '2025-04-21'; % Easter Monday (2. Påskedag)
    '2025-05-01'; % Labour Day (Første Maj)
];

t_dates = datetime(t, 'ConvertFrom', 'excel');
Maturity_date = datetime(Maturity, 'ConvertFrom', 'excel');

% Initialize T vector
T = zeros(size(t));
IV = zeros(size(t));
vega = zeros(size(t));
gamma = zeros(size(t));
delta_call = zeros(size(t));
delta_put = zeros(size(t));

% Calculate business days for each date in t
for i = 1:length(t)
    T(i) = days252bus(t_dates(i), Maturity_date, red_days)/252;
    IV(i)= Implicit_sigma(BS_Call(i),S(i),K(i),r(i),q(i),T(i),0.00001, 0.2);
    [delta, vega_, gamma_] = CalculateGreeks(S(i), K(i), r(i), q(i), T(i), IV(i), true);
    delta_call(i) = delta;
    vega(i) = vega_;
    gamma(i) = gamma_;

    [delta, vega_, gamma_] = CalculateGreeks(S(i), K(i), r(i), q(i), T(i), IV(i), false);
    delta_put(i) = delta;
end

    %writematrix(IV,'portfolioManagerV4Projekt3.xlsm', 'Sheet', 'Answer', 'Range', 'L8:L13')
    %writematrix(delta_call,'portfolioManagerV4Projekt3.xlsm', 'Sheet', 'Answer', 'Range', 'M8:M13')
    %writematrix(delta_put,'portfolioManagerV4Projekt3.xlsm', 'Sheet', 'Answer', 'Range', 'N8:N13')
    %writematrix(gamma,'portfolioManagerV4Projekt3.xlsm', 'Sheet', 'Answer', 'Range', 'O8:O13')
    %writematrix(vega,'portfolioManagerV4Projekt3.xlsm', 'Sheet', 'Answer', 'Range', 'P8:P13')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 5 in step by step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noOfSimulations = 500;      % Number of scenarios to generate
confidence_level = 0.95;    % For 95% VaR
risk_limit = 1e6;           % 1 million DKK VaR limit

portfolio_values = zeros(noOfSimulations, 1);
option_values_call = zeros(noOfSimulations, 1);
option_values_put = zeros(noOfSimulations, 1);


index_prices = readmatrix('portfolioManagerV4Projekt3.xlsm', 'Sheet', 'Refinitiv', 'Range', 'AF6:AF2086');

log_returns = diff(log(index_prices));
VaR_95 = zeros(length(t),1);
ratio_puts_to_calls = zeros(length(t),1);
num_calls = zeros(length(t),1);
num_puts = zeros(length(t),1);
value_options = zeros(length(t),1);
hedge_delta = zeros(length(t),1);

portfolio_history = struct();
nCall_fixed = 0;
nPut_fixed = 0;

for day = 1:length(t)

    current_price = S(day);
    strike = K(day);
    rate = r(day);
    div_yield = q(day);
    time_to_maturity = T(day);

    % Each day, you simulate 500 scenarios of the index (based on past returns and current index.
    random_indices = randi(numel(log_returns), noOfSimulations, 1);
    simulated_log_returns = log_returns(random_indices);
 
    % Convert to prices
    simulated_prices = current_price * exp(simulated_log_returns);

    call_values = zeros(noOfSimulations, 1);
    put_values = zeros(noOfSimulations, 1);
    call_deltas = zeros(noOfSimulations, 1);
    put_deltas = zeros(noOfSimulations, 1);

    %Calculate the values of the portfolio of call/put based on the simulated index returns.
    
    for i = 1:noOfSimulations    
        call_values(i) = BlackScholes(simulated_prices(i), strike, rate, div_yield, time_to_maturity, IV(day), true);
        put_values(i) = BlackScholes(simulated_prices(i), strike, rate, div_yield, time_to_maturity, IV(day), false);
    end
    
    %Calculate delta
    [current_call_delta, ~, ~] = CalculateGreeks(current_price, strike, rate, div_yield, time_to_maturity, IV(day), true);
    [current_put_delta, ~, ~] = CalculateGreeks(current_price, strike, rate, div_yield, time_to_maturity, IV(day), false);

    % Start with equal quantities (will be scaled)
    
    ratio_puts_to_calls(day) = -current_call_delta / current_put_delta;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate VaR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Current portfolio value
    current_call_price = BlackScholes(current_price, strike, rate, div_yield, time_to_maturity, IV(day), true);
    current_put_price = BlackScholes(current_price, strike, rate, div_yield, time_to_maturity, IV(day), false);

    if day == 1
        % Determine optimal nCall and nPut on day 1 only
        max_nCall = 0;
        best_VaR = 0;

        for test_nCall = 0:1:100000
            test_nPut = round(test_nCall * ratio_puts_to_calls(day));

            current_port_value = test_nCall * current_call_price + test_nPut * current_put_price;
            scenario_port_values = -test_nCall * call_values - test_nPut * put_values;
            pnl = scenario_port_values - current_port_value;

            sorted_losses = sort(pnl);
            test_VaR = -sorted_losses(round((1 - confidence_level) * noOfSimulations));

            if test_VaR <= risk_limit && test_nCall > max_nCall
                max_nCall = test_nCall;
                best_nPut = test_nPut;
                best_VaR = test_VaR;
            end
        end

        nCall_fixed = max_nCall;
        nPut_fixed = best_nPut;

        nCall = max_nCall;
        nPut = best_nPut;
    else
        % Use the same options position from day 1
        nCall = nCall_fixed;
        nPut = nPut_fixed;
    end

    % Always compute VaR
    current_port_value = nCall * current_call_price + nPut * current_put_price;
    scenario_port_values = -nCall * call_values - nPut * put_values;
    pnl = scenario_port_values - current_port_value;
    sorted_losses = sort(pnl);
    VaR_95(day) = -sorted_losses(round((1 - confidence_level) * noOfSimulations));

    portfolio_history(day).date = t_dates(day);
    num_calls(day) = nCall;
    num_puts(day) = nPut;
    portfolio_history(day).VaR = VaR_95(day);
    hedge_delta(day) = nCall * current_call_delta + nPut * current_put_delta;
    value_options(day) = nPut * current_put_price + nCall * current_call_price;
    
end

   % 
   % writematrix(VaR_95,'portfolioManagerV4Projekt3.xlsm', 'Sheet', 'Answer', 'Range', 'R8:R13')
   % writematrix(num_calls,'portfolioManagerV4Projekt3.xlsm', 'Sheet', 'Answer', 'Range', 'F8:F13')
   % writematrix(num_puts,'portfolioManagerV4Projekt3.xlsm', 'Sheet', 'Answer', 'Range', 'G8:G13')
   % writematrix(value_options,'portfolioManagerV4Projekt3.xlsm', 'Sheet', 'Answer', 'Range', 'E8:E13')
   % writematrix(hedge_delta,'portfolioManagerV4Projekt3.xlsm', 'Sheet', 'Answer', 'Range', 'Q8:Q13')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function implicit_volatility = Implicit_sigma(BS, S0, K, r, q, T, tol, guess_sigma)


CallOptionPrice = BlackScholes(S0, K, r, q, T, guess_sigma, true);

upper_sigma = 1;
lower_sigma = 0.01;
current_sigma = guess_sigma;

while (abs(CallOptionPrice - BS) > tol)
    if CallOptionPrice < BS
        lower_sigma = current_sigma;
    else 
        upper_sigma = current_sigma;
    end
    current_sigma = (upper_sigma + lower_sigma) / 2;
    CallOptionPrice = BlackScholes(S0, K, r, q, T, current_sigma,true);
end

implicit_volatility = current_sigma;

end

function BS = BlackScholes(S0, K, r, q, T, current_sigma, isCall)
    
d1 = (log(S0/K) + (r -q + (current_sigma^2)/2)*T) / (current_sigma*sqrt(T));
d2 = d1 - current_sigma*sqrt(T);

if isCall
    BS = S0*exp(-q*T)*normcdf(d1) - K*exp(-r*T)*normcdf(d2);
else
    BS = K*exp(-r*T)*normcdf(d2) - S0*exp(-q*T)*normcdf(d1);
end

end

function [delta, vega, gamma] = CalculateGreeks(S, K, r, q, T, sigma, isCall)
    % Calculate d1 and d2 for Black-Scholes formulas
        d1 = (log(S/K) + (r -q + (sigma^2)/2)*T) / (sigma*sqrt(T));
    
    % Delta calculation
    if isCall
        delta = exp(-q*T) * normcdf(d1);
    else
        delta = -exp(-q*T) * normcdf(-d1);
    end
    
    % Vega calculation (same for calls and puts)
    vega = S * exp(-q*T) * normpdf(d1) * sqrt(T);
    
    % Gamma calculation (same for calls and puts)
    gamma = exp(-q*T) * normpdf(d1) / (S * sigma * sqrt(T));
end
