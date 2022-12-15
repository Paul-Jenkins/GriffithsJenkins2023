function sigmax = diffusionsigma(xin) % x = [x11 x12; x21; x22]
% Function to compute the diffusion sigma on p19 of the paper. Input: a 2x2
% matrix x.

x = reshape(xin,[size(xin,1)*size(xin,2),1]); % x = [x11; x21; x12; x22];
assert(all(size(x) == [4,1]));
sigmax = diag(sqrt(x)) - diag(x)*ones(length(x))*diag(sqrt(x)); % Pal (2011, eq (5))
assert(all(size(sigmax) == [4,4]));
end