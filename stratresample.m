function I = stratresample(p , N)

% Stratified resampling for sequential Monte Carlo. 
% Sample N times with repetitions from a distribution on [1:length(p)] with probabilities p.

% See http://www.cornebise.com/julien/publis/isba2012-slides.pdf

p = p/sum(p);  % normalize, just in case...

cdf = cumsum(p) ;
cdf(end) = 1 ;
I = zeros(N,1) ;
U = rand(N,1) ;
U = U/N + (0 : (N - 1))'/N;
index = 1 ;
for k = 1 :N
    while (U(k) > cdf(index))
      index = index + 1 ;
    end
    I(k) = index ;
end

end

