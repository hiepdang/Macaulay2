needsPackage "Schubert2"
-- computes the degree of Fano schemes of k-dimensional linear subspaces
-- on a general hypersurface of degree d in P^n.

degreeFanoScheme = (k,d,n) -> (
     p := (k+1)*(n-k) - binomial(d+k,d);
     G := flagBundle {n-k,k+1};
     s := chern symmetricPower_d last bundles G;
     q := chern_1 last bundles G;
     integral(s*q^p)
     )
