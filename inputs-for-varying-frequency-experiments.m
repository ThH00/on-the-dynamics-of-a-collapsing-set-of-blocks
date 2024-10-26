% low frequency
omega = pi*(0:0.25:5)';
n_osc = 300/(2*pi)*omega;

omega = pi*(0.05:0.05:0.2)';
n_osc = 300/(2*pi)*omega;

omega = pi*(6:1:10)';
n_osc = 300/(2*pi)*omega;

% high frequency
omega2 = pi*(0:5:100)';
n_osc2 = 120/(2*pi)*omega2;