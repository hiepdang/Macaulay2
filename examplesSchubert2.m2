needsPackage "Schubert2"

degreeFanoScheme = (k,d,n) -> (
-- computes the degree of Fano schemes of k-dimensional linear subspaces
-- on a general hypersurface of degree d in P^n.
    p := (k+1)*(n-k) - binomial(d+k,d);
    G := flagBundle {n-k,k+1};
    s := chern symmetricPower_d last bundles G;
    q := chern_1 last bundles G;
    integral(s*q^p)
    )

GWLines = (k,a,b) -> (
-- computes the Gromov-Witten invariants for lines as in the paper of Katz
-- http://arxiv.org/abs/alg-geom/9301006
    c := k - a - b;
    G := flagBundle {2,k};
    s := chern symmetricPower_(k+2) dual first bundles G;
    qa := chern_(a-1) last bundles G;
    qb := chern_(b-1) last bundles G;
    qc := chern_(c-1) last bundles G;
    integral(s*qa*qb*qc)
    )
