% The function task(id) creates problem instance id with randomized values 
% (stored in d), the structure of the solution is stored in s
% The function task(id, d, s) checks if the solution s is correct

result = '\nThis is the result that should be appended to the report:\n';

[statusString, d, s] = task(1);  %  Randomize data
m = d.m; % m contains the market data
p = d.p; % p contains the portfolio data
datestr(m.dates); % contains the dates where the last 5 days are the corresponding trading dates. 
datestr(p.dates); % contains the trading dates that buy p.buySellEquity. 

fprintf('Market data:\n'); fprintf(m.readMe);
fprintf('\nPortfolio data:\n'); fprintf(p.readMe);
fprintf('\nTask 1:\n'); fprintf(statusString); fprintf('\n');

% Check the structure of the solution in s!

proportionalCommissionEquity  = d.proportionalCommissionEquity; % Commission
% Get structure of the required calculation
% cashEquityTradesForeign = 5 trading days x 5 assets 
% commissionEquityTradesForeign = 5 trading days x 5 commissions 
cashEquityTradesForeign       = s.cashEquityTradesForeign; 
commissionEquityTradesForeign = s.commissionEquityTradesForeign; 

% Enter code for task 1 here 

for day = 1:length(p.dates)
    for asset = 1:size(p.buySellEquity, 2)
        antal = p.buySellEquity(day, asset);
        marketIdx = find(m.dates == p.dates(day));
        % Skapa en map mellan datum och radindex

        if antal > 0
            % Buy at ask price
            price = m.askPrices(marketIdx, asset);
            cashFlow = -antal * price;
        elseif antal < 0
            % Sell at bid price
            price = m.bidPrices(marketIdx, asset);
            cashFlow = -antal * price;  % antal < 0, so cashFlow is positive
        else
            cashFlow = 0;
            price = 0;
        end

        commission = abs(cashFlow) * d.proportionalCommissionEquity;

        cashEquityTradesForeign(day, asset) = cashFlow;
        commissionEquityTradesForeign(day, asset) = commission;
    end
end

% After that assign the results to s

s.cashEquityTradesForeign       = cashEquityTradesForeign;
s.commissionEquityTradesForeign = commissionEquityTradesForeign;

% Verify the results
[statusString] = task(1, d, s);
result = strcat(result, statusString);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[statusString, d, s] = task(2); %  Randomize data

m = d.m; % m contains the market data
p = d.p; % p contains the portfolio data

fprintf('\nTask 2:\n'); fprintf(statusString); fprintf('\n');

% Get structure of the required calculation
% which business day (fxForwardMaturityDate) is 1 month from the trading date (p.date)
% what is the Actual/365 to maturity
% F = S exp((i_d - i_f)t) where F,S are the forward and spot FX calculated in m. 
fxForwardMaturityDate   = s.fxForwardMaturityDate; % Work with the simplification that spot rate is tlength(p.dates);= length(p.dates);oday
tAct365                 = s.tAct365; 
impliedInterestRateDiff = s.impliedInterestRateDiff; 


% Enter code for task 2 here



for i = 1:length(p.dates)
    tradeDate =  datetime(p.dates(i), 'ConvertFrom', 'datenum');

    % 1. Bestäm förfallodatum (lägg till 1 månad)
    maturityDate = tradeDate + calmonths(1);

    % Om helg, flytta till måndag
    while isweekend(maturityDate)
        maturityDate = maturityDate + 1;
    end
    fxForwardMaturityDate(i) = datenum(maturityDate);

    % 2. Beräkna t (faktisk tid i år)
    tAct365(i) = days(maturityDate - tradeDate) / 365;

    % 3. Hitta index i m.dates som motsvarar tradeDate
    marketIdx = find(m.dates == p.dates(i));

    % 4. Spotkurs = snitt mellan bid/ask
    spotFX = (m.bidExchangeRate(marketIdx) + m.askExchangeRate(marketIdx)) / 2;

    % 5. Spread = snitt mellan bid/ask (i pips → konvertera till faktisk kurs)
    spreadPips = (m.bidFxForward1M(marketIdx) + m.askFxForward1M(marketIdx)) / 2;
    forwardFX = spotFX + spreadPips * 1e-4;

    % 6. Räkna ut implicita ränteskillnaden
     impliedInterestRateDiff(i) = log(forwardFX / spotFX) / tAct365(i);
   
end


% After that assign the results to s
s.fxForwardMaturityDate   = fxForwardMaturityDate;
s.tAct365                 = tAct365;
s.impliedInterestRateDiff = impliedInterestRateDiff;


[statusString] = task(2, d, s);
result = strcat(result, statusString);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[statusString, d, s] = task(3); % Randomize data
m = d.m;
p = d.p;
fprintf('\nTask 3:\n'); fprintf(statusString); fprintf('\n');
% --- Läs in data ---

impliedInterestRateDiff = d.impliedInterestRateDiff; % 5x1 (en diff per dag i)
fxForwardMaturityDate = d.fxForwardMaturityDate; % 4x1 (en maturity per kontrakt j)
fxForwardK = d.fxForwardK; % 4x1 (strike för kontrakt j)
nDays = length(p.dates); % oftast 5
nContracts = 4; % forward-kontrakt fre-ons

% --- Initiera output ---
interestRateSek1Mcont = zeros(nDays,1);
fxForwardPrices = zeros(nDays,nContracts);
fxForwardValue = zeros(nDays,nContracts);
% --- Loopa över varje värderingsdag (i) ---

for i = 1:nDays
        % 1) Hitta index i marknadsdata för dag i
        marketIdx_i = find(m.dates == p.dates(i));
        % 2) Spot på dag i
        spot_i = (m.bidExchangeRate(marketIdx_i) + m.askExchangeRate(marketIdx_i)) / 2;
        % 3) Räkna om simple -> kontinuerlig SEK-ränta (360 -> 365) för dag i
        tradeDate_i = datetime(p.dates(i), 'ConvertFrom', 'datenum');
        maturityDate_i = tradeDate_i + calmonths(1);
    while isweekend(maturityDate_i)
     maturityDate_i = maturityDate_i + 1;
    end
        dagar = days(maturityDate_i - tradeDate_i);
        t360 = dagar / 360;
        t365 = dagar / 365;
        rSimple_i = m.interestRateSek1M(marketIdx_i);
        rCont_i = log(1 + rSimple_i * t360) / t365;
        interestRateSek1Mcont(i) = rCont_i;
        % 4) Ränteskillnad för dag i
        rateDiff_i = impliedInterestRateDiff(i);
        % --- Loopa över alla forwardkontrakt j (fre-ons) ---
    for j = 1:nContracts
    % Om j > i => kontraktet är inte skrivet ännu
        if j > i
            continue;
        end
    % a) Räkna ut T = tid från dag i -> kontraktets förfall
    maturityDate_j = datetime(fxForwardMaturityDate(j), 'ConvertFrom', 'datenum');
    T = days(maturityDate_j - tradeDate_i) / 365;
        if T > 0
            % b) Forwardpris = spot_i * exp( rateDiff_i * T )
            F = spot_i * exp(rateDiff_i * T);
            fxForwardPrices(i,j) = F;
            % c) Värde = (F - K_j) * e^(-rCont_i * T)
            K_j = fxForwardK(j);
            fxForwardValue(i,j) = (F - K_j) * exp(-rCont_i * T);
        else
            % Kontraktet har förfallit (T <= 0)
            % => värdet = spot_i - K_j (vid förfall)
            K_j = fxForwardK(j);
            fxForwardPrices(i,j) = spot_i;
            fxForwardValue(i,j) = spot_i - K_j;
        end
    end
end
% --- Tilldela struktur ---
s.interestRateSek1Mcont = interestRateSek1Mcont;
s.fxForwardPrices = fxForwardPrices;
s.fxForwardValue = fxForwardValue;
% --- Verifiera ---
[statusString] = task(3, d, s);
result = strcat(result, statusString);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[statusString, d, s] = task(4);%  Randomize data

m = d.m; % m contains the market data
p = d.p; % p contains the portfolio data (Friday - Thursday)

fprintf('\nTask 4:\n'); fprintf(statusString); fprintf('\n');

lotSizeFutures     = d.lotSizeFutures;

% Get structure of the required calculation
% equityBeta = beta of 5 days x 5 assets
equityBeta         = s.equityBeta;
equityValueForeign = s.equityValueForeign;
equityWeights      = s.equityWeights;
equityHoldingBeta  = s.equityHoldingBeta;
% On Tuesday (3rd day), using future to hedge 
nEquityFutures     = s.nEquityFutures;
% On Wednesday (4th day) settlement of future contract
settlementFutures  = s.settlementFutures;


% Enter code for task 4 here

% After that assign the results to s
s.equityBeta         = equityBeta;
s.equityValueForeign = equityValueForeign;
s.equityWeights      = equityWeights;
s.equityHoldingBeta  = equityHoldingBeta;
s.nEquityFutures     = nEquityFutures;
s.settlementFutures  = settlementFutures;


[statusString] = task(4, d, s);
result = strcat(result, statusString);

fprintf(result);





% Enter code for task 1 here



