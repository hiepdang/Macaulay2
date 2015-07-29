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

--Calculate the number of lines in P^3 meeting each of 4 lines
restart
loadPackage "Schubert2"
G = flagBundle({2,2})            -- Grassmannian of lines in P^3
sigma_1 = schubertCycle(G,{1,0}) -- Create schubert cycle sigma_1
integral sigma_1^4               -- Fourfold self intersection

--Calculate the number of lines on a cubic surface
restart
loadPackage "Schubert2"
G = flagBundle({2,2})            -- Grassmannian of lines in P^3
(S,Q) = G.Bundles                -- Tautological sub- and quotient-bundles
B = symmetricPower(3,Q)          -- 3rd symmetric power
integral chern B                 -- top chern class

--Calculate the number of lines on a quintic threefold
restart
loadPackage "Schubert2"
G = flagBundle({3,2})            -- Grassmannian of lines in P^4
(S,Q) = G.Bundles                -- Tautological bundles					    
B = symmetricPower(5,Q)          -- 5th symmetric power
integral chern B                 -- top chern class

--In a pencil of quartic surfaces in P^3, how many contain a line?
restart
loadPackage "Schubert2"
G = flagBundle({2,2})
(S,Q) = G.Bundles
B = symmetricPower(4,Q)
integral chern B					    

--Calculate the number of conic curves in P^3 meeting each of 8 lines
restart
loadPackage "Schubert2"
G = flagBundle({1,3})                                --Projective space of planes in P^3
h = schubertCycle(G,{1,0,0})                         --hyperplane class
(S,Q) = G.Bundles                                    --Tautological bundles  
B = symmetricPower(2,Q)                              --2nd symmetric power
H = projectiveBundle(dual B,VariableNames =>{,{z}})  --Space H of conics, with z = c(O(1))
integral (2*h+z)^8                                   --8-fold self int.

--Calculate the number of conic curves in P^4 meeting each of 11 2-planes
restart
loadPackage "Schubert2"
G = flagBundle({2,3})                                -- Grassmannian of 2-planes in P^4
sigma_1 = schubertCycle(G,{1,0,0})                   -- Schubert cycle
(S,Q) = G.Bundles                                    -- Tautological bundles
B = symmetricPower(2,Q)                              -- 2nd symmetric power
H = projectiveBundle(dual B,VariableNames =>{,{z}})  -- Projective bundle of B/G; tautological class z
integral (2*sigma_1+z)^11                            -- 11-fold self-int. of relevant class

--Calculate the number of singular members in a pencil of degree 3 curves
restart
loadPackage "Schubert2"
G = flagBundle({2,1})               -- Projective space P^2
(S,Q) = G.Bundles                   -- Taut. bundles
F = (dual(tangentBundle G))*(Q^**3) -- Bundles filtering E
G = Q^**3                           
E = F+G                             -- Can think of E as F+G
integral chern E                    -- Top chern class of E

--Calculate the number of singular members in a pencil of degree 4 curves
restart
loadPackage "Schubert2"
G = flagBundle({2,1})               -- Projective space P^2
(S,Q) = G.Bundles                   -- Taut. bundles
F = (dual(tangentBundle G))*(Q^**4) -- Bundles filtering E
G = Q^**4
E = F+G                             -- Can think of E as F+G
integral chern E                    -- Top chern class of E

--Calculate the number of lines meeting a given general quintic surface in P^3 in only one point.
restart
loadPackage "Schubert2"
G = flagBundle({2,2})				     -- Grassmannian of lines in P^3
(S,Q) = G.Bundles
Phi = projectiveBundle(dual Q,VariableNames =>{s,q}) -- Universal line over grassmannian
(S1,Q1) = Phi.Bundles				     -- Taut. bundles on Phi
vTBd = dual Hom(S1,Q1)				     -- Dual of the vertical tangent bundle
O5 = Q1^**5					     -- O_Phi(5)
L0 = O5						     -- Bundles filtering E
L1 = vTBd * O5
L2 = vTBd^**2 * O5
L3 = vTBd^**3 * O5
L4 = vTBd^**4 * O5
E = L0 + L1 + L2 + L3 + L4			     -- E can be thought of as the sum of the line bundles in the filtration
integral chern E				     -- Top chern class of E

--Calculate the number of lines meeting a given general septic hypersurface in P^4 in only one point
restart
loadPackage "Schubert2"
G = flagBundle({3,2})
(S,Q) = G.Bundles
Phi = projectiveBundle(dual Q,VariableNames =>{s,q})
(S1,Q1) = Phi.Bundles
vTBd = dual Hom(S1,Q1)
O7 = Q1^**7
L0 = O7
L1 = vTBd * O7
L2 = vTBd^**2 * O7
L3 = vTBd^**3 * O7
L4 = vTBd^**4 * O7
L5 = vTBd^**5 * O7
L6 = vTBd^**6 * O7
E = L0 + L1 + L2 + L3 + L4 + L5 + L6
integral chern E
