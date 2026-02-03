-- PFD Algorithm Applied to Wavefunction Coefficient
-- define polynomial ring
R = QQ[x_1,x_2,x_3,x_4,y_14,y_24,y_34];

-- define denominator hyperplanes and the number of them, n
n = 11;
l_1 = x_1 + y_14;
l_2 = x_2 + y_24;
l_3 = x_3 + y_34;
l_4 = x_4 + y_14 + y_24 + y_34;
l_5 = x_1 + x_2 + x_3 + x_4;
l_6 = x_1 + x_4 + y_34 + y_24;
l_7 = x_2 + x_4 + y_14 + y_34;
l_8 = x_3 + x_4 + y_14 + y_24;
l_9 = x_1 + x_2 + x_4 + y_34;
l_10 = x_3 + x_4 + x_1 + y_24;
l_11 = x_3 + x_4 + x_2 + y_14;

-- define numerator
f = 8*y_14*y_24*y_34*(l_7*l_8*l_10*l_11+l_6*l_8*l_10*l_11+l_7*l_8*l_9*l_11+l_6*l_7*l_9*l_11+l_6*l_8*l_9*l_10+l_6*l_7*l_9*l_10);

-- check if numerator has any initial spurious poles and divide out by any that exist
spuriousls = {};
for q from 1 to n do(if isMember(f,ideal(l_q))==true then (spuriousl = append(spuriousls,l_q), f = lift(f/l_q,R)))
#spuriousls

-- find maximum d for which f is in the ideal I_{d,L}
for d from 1 to n do(
    Igens = {}, for lineindices in subsets(n,d) do(
        lprod = 1, for lineindex in lineindices do(lprod = lprod*l_(lineindex+1)),
         Igens = append(Igens,lprod)), I = ideal(Igens), if isMember(f,I) == false then break else (print d, finalideal=I, finald = d))

-- generating list of non-zero numerators and matching them with their corresponding l_{i_1}*...*l_{i_finald} factors
numeratorlist = flatten(entries(f // (gens finalideal)));
entrynumber = -1;
labellednumeratorlist = {};
for lineindices in subsets(n,finald) do(
    entrynumber = entrynumber + 1, if numeratorlist#(entrynumber) != 0 then (labellednumeratorlist = append(labellednumeratorlist,{lineindices,numeratorlist#(entrynumber)})))

labellednumeratorlist