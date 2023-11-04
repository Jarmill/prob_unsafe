function moments = NormalMom(dmin, dmax)    
%moments 0..d of univariate unit (0,1) normal distribution
    if nargin < 2
        dmax = dmin;
        dmin = 0;
    end
    mom_hits = dmin:dmax;
    moments = zeros(length(mom_hits),1);
    for i = 1:numel(moments)
        j = mom_hits(i);
        if mod(j, 2)
            moments(i)=0;
        else
            k = j/2;
            %https://math.stackexchange.com/questions/1796806/calculate-ex2n-where-x-is-normal-0-1?noredirect=1&lq=1
            moments(i) = (nchoosek(2*k, k)*factorial(k)/2^(k));
        end
    end
end