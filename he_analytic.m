syms p q lambda1 lambda2 mu nu beta
private_cost(p, lambda1, mu) = 1 / (mu - lambda1 + lambda1 * p);
public_cost(p, q, lambda1, lambda2, mu, nu) = (1 / (nu - lambda1 * p - lambda2* q));
cost(p, q, lambda1, lambda2, mu, nu) = p * public_cost(p, q, lambda1, lambda2, mu, nu) + (1-p) * private_cost(p, lambda1 , mu);
he_term(p, q, lambda1, lambda2, mu, nu, beta) = beta/2 * (cost(q, p, lambda2, lambda1 , mu, nu) - cost(p, q, lambda1, lambda2, mu, nu));
he_cost(p, q, lambda1, lambda2, mu, nu, beta) = cost(p, q, lambda1, lambda2, mu, nu) + he_term(p, q, lambda1, lambda2, mu, nu, beta);

diff_he_cost(p, q, lambda1, lambda2, mu, nu, beta) = diff(he_cost(p, q, lambda1, lambda2, mu, nu, beta), p);
sol = solve(diff_he_cost(p, q, lambda1, lambda2, mu, nu, beta) == 0, p);
pstar(q, lambda1, lambda2, mu, nu, beta) = sol(2);
qstar(p, lambda1, lambda2, mu, nu, beta) = pstar(p, lambda2, lambda1, mu, nu, beta);

solve(pstar(qstar(p, lambda1, lambda2, mu, nu, beta), lambda1, lambda2, mu, nu, beta) - p == 0, p)