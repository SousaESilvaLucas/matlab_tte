function [ h ] = complexfirminphase(b)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Make sure b is a row
b = b(:).';

% Find the roots of the filter polynomial
r = roots(b);

% Find the unit circle zeros
    tol = 1e-15;
    [ru,tol,Nsz] = findUnitCircleZeros(b,r,tol);

% Find the zeros strictly inside the unit circle
ri = removeUnitCircleZeros(r,Nsz);

% Form vector of single roots on unit circle and all roots inside
rm = [ru;ri];

% Form polynomial
h = poly(leja(rm));

% Apply gain correction
desired = sqrt(max(abs(freqz(b))));
actual = max(abs(freqz(h)));
h = h * (desired/actual);

function ri = removeUnitCircleZeros(r,Nsz)
% Get zeros strictly inside unit circle

% Get total number of zeros
N = length(r);

% Sort roots by magnitude
[dummy,k] = sort(abs(r));

% Get N/2 - Nsz smallest roots
rs = r(k);
ri = rs(1:N/2-Nsz);
end

function [ru,tol,Nsz] = findUnitCircleZeros(b,r,tol,Nsz)
% Find single multiplicity zeros on the unit circle


% Compute the zeros of the derivative
% This will contain the same zeros on the unit circle as b, but not double
rd = roots(polyder(b));

% Initialize flag indicating whether to use the zeros from the derivative
% polynomial or not. Don't use derivative zeros by default.
derivflag = 0;

% Initialize flag
flag = 1;

while tol < 1e-3 && flag,
	% Find the zeros on the unit circle
	ru = r(find(abs(abs(r)-1) < tol));
		
	% Find the zeros on the unit circle corresponding to the derivative
	rdu = rd(find(abs(abs(rd)-1) < tol));
	
    if nargin > 3,
        % Check if number of single zeros found corresponds to number
        % given number of zeros
        if length(rdu) >= Nsz,
            rdu = rdu(1:Nsz);
            flag = 0;
            derivflag = 1;
        elseif ~rem(length(ru),2) && length(ru)/2 == Nsz 
            flag = 0;
        else
            % Decrease the tolerance
            tol = 10*tol;
        end
    elseif ~rem(length(ru),2) && length(rdu) == length(ru)/2 && length(ru) ~= 0,
        % Check that the length of the derivative zeros half of the length of
        % the original zeros on the upper unit circle. If not, decrease
        % tolerance and try again
		flag = 0;
        Nsz = length(rdu);
        if max(abs(polyval(b,ru))) > max(abs(polyval(b,rdu)))
            derivflag = 1;
        end
	else
		% Decrease the tolerance
		tol = 10*tol;
	end
end

% If max tolerance has been reached and number of unit circle zeros has not
% been specified, use the zeros on the unit circle found from the original
% polynomial (not the derivative)
if tol >= 1e-3 && nargin < 4,
    Nsz = length(ru)/2;
end  

if ~isempty(ru)
    if derivflag,
        % Use unit circle zeros from derivative
        ru = rdu;
    else
        ru = getsinglezeros(ru,Nsz);
    end
    
    if tol < 1e-3,
        % Force length of zeros to be exactly one, keep the angle the same
        ru = exp(i*angle(ru));
    end
else
    Nsz = 0;
end

end

end

