restart
needsPackage "Matroids"

--input: field kk, number of variables r
--output: polynomial ring
defineRing = (kk, r) -> (
    kk[x_1..x_r]
)


--input: (r x n) matrix of coefficients for the n lines
--output: list of the n lines
defHypArr = M -> (
    y = matrix{gens R};
    line = y * M;
    linelist = flatten entries line
)

--input: (r x n) matrix of coefficients for the n lines
--output: matroid associated to hyperplane arrangement
defMatroid = M -> (
    matroid M
)

--input: (r x n) matrix of coefficients for the n lines, d
--output: ideal I_L,d
defILd = (M,d) -> (
    indL = subsets(n,d);
    L = defHypArr(M);
    n = #L;
    Igens = {};
    for S in indL do(
	f = 1;
	for j in S do(
	    f = f * (L_j));
	Igens = append(Igens, f););
    ideal(Igens)
)

--input: (r x n) matrix of coefficients for the n lines, d
--output: checks primary decomposition in Corollary 2.9
checkPriDecomp = (M,d) -> (
    L = defHypArr(M);
    n = #L;
    r = dim R;
    ILd = defILd(M,d);
    matM = defMatroid(M);
    goodF = {};
    for F in flats(matM) do(
	if #F > (n-d) then goodF = append(goodF, F)
	);
    flatI = {};
    for F in goodF do(
	Jgens = {};
	for i in toList(F) do(
	    Jgens = append(Jgens, L_i););
	J = (ideal(Jgens))^(d-n+#F);
	flatI = append(flatI, J););
    intersect(flatI) == ILd
)    

--Algorithm 1
--input: (r x n) matrix of coefficients for the n lines, numerator of rational function f
--output: new list of lines, numerator
reduceInput = (M, f) -> (
    L = defHypArr(M);
    newL = {};
    for i in L do(
	I = ideal(i);
	if isMember(f,I) then f = f // i
	else newL = append(newL, i)
	);
    newL, f
)

--Algorithm 2
--input: reduced rational function in the form of (r x n) matrix of coefficients for the n lines,
--numerator of rational function f
--output: maximum d for which PFD of degree d exists, associated PFD via hashtable where keys are indices of linear factors
--appearing in the denominator, associated value is the corresponding numerator
PFD = (M, f) -> (
    L = defHypArr(M);
    n = #L;
    d = 0;
    for D from 1 to n do(
	ILd = defILd(M,D);
	if isMember(f,ILd) then(
	    d = D)
	else break
	);
    if d>0 then(
	ILd = defILd(M,d);
	coeffs = flatten entries (f // (gens ILd));
	denoms = apply(subsets(n,d), S->toList(0..n-1) - set S);
	coeffTable = {};
	for i from 0 to (binomial(n,d)-1) do(
	    coeffTable = append(coeffTable, denoms_i => coeffs_i)
	    );
	coeffHash = new HashTable from coeffTable;
	return (d, coeffHash)
	)
    else return "no PFD exists"
)
    


---example usage

R = defineRing(QQ,3);
M = transpose matrix{{1,0,0},{0,1,0},{0,0,1},{1,-1,0},{1,0,-1},{0,1,-1}}

checkPriDecomp(M,4) --confirms primary decomposition of I_L,4 given in Corollary 2.9

f = random(4,ideal(gens R)) --random degree 4 homogeneous polynomial
PFD(M,f) --returns PFD of degree d, with d maximal
