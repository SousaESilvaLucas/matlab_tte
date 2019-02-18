function [Z] = noiseTTE( a,imp_2,num,Pn_dB )

% This function generates a TTE noise vector Z = (X + Y).*exp(1i*Theta), where X is the 
% rayleigh random variable, Y is the Weibull random variable and Theta is a
% uniform random variable defined in [0,2*pi]

% a = weibull shape parameter (0 < a < 2)
% imp_2 = square of the impulsivity parameter (imp_2 > 0)
% imp_2 = E(X.^2)/E(Y.^2)
% Pn_dB = Noise Power in dBW
% num = Length of the noise vector

% --- calculates the noise power in W:
	
Pn = 10.^(Pn_dB/10);

% --- calculates R, the scale parameter of the Weibull random variable

KTE = gamma(1 + 2/a)/imp_2 + sqrt(pi)*gamma(1 + 1/a)*sqrt(gamma(1 + 2/a))/sqrt(imp_2) + gamma(1 + 2/a);
R = sqrt(Pn/KTE);

% --- calculates sigma_0, the scale parameter of the Rayleigh random variable

R_0 = R*sqrt(gamma(1 + 2/a))/sqrt(imp_2);
sigma_0 = R_0/sqrt(2);

% --- generates the TTE noise vector

X = raylrnd(sigma_0,num,1);
Y = wblrnd(R,a,num,1);
Theta = 2*pi*rand(num,1);
Z = (X + Y).*exp(1i*Theta);

% potX = mean(X.^2);
% 
% potY = mean (Y.^2);
% 
% P_teoricoX = R_0.^2;
% P_teoricoY = gamma(1 + 2/a)*(R.^2);
% potZ = mean (Z.^2);

end


