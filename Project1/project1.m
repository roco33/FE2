% initial stock price
S0 = 90.74;
% risk-free rate
r = 0.01;
% dividend yield
div = 0.02;
% volatility
sigma = 0.21;
% down side threshold level 
K = 72.592;
% redemption price 
red = [95.28, 99.81, 104.35];
% number of steps
N = 900;
% stock tree
VStock = zeros(1000);
% option tree 
VOption = zeros(1000);
% determination date
barr = N/12;


% time step 
delta = 3/N;

u = exp((r-div)*delta + sigma * sqrt(delta));
d = exp((r-div)*delta - sigma * sqrt(delta));

qu = (exp((r-div)*delta)-d)/(u-d);
qd = 1- qu;

for i = 1:N+1
	VStock(N+1,i) = S0 * u^(i-1) * d^(N+1-i);
	if VStock(N+1,i) >= K
		VOption(N+1,i) = 10.2125;
	else
		VOption(N+1,i) = 10 * VStock(N+1,i) / S0;
	end
end
for j = N:-1:1
	for i = j:-1:1
		VStock(j,i) = S0 * u^(i-1) * d^(j+1-i);
		VOption(j,i) = qu * VOption(j+1,i+1) + qd * VOption(j+1,i);
		for k = 1:3
			if ismember(j-1,(k-1)*N/3+[barr*1,barr*2,barr*3,barr*4])
				if VStock(j,i) > K && VStock(j,i) < red(k)
					VOption(j,i) = VOption(j,i) + 0.2125;
				elseif VStock(j,i) > red(k)
					VOption(j,i) = 10.2125;
				end
			end
		end
	end
end

VOption(1,1)
	