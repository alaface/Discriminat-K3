// DiscriminantGroup
// Input: A matrix with integer coefficients
// Output: Its discriminant group and the quadratic form in Q/2Z

function DiscriminantGroup(M)
    // Compute the Smith Normal Form of the input matrix
    S, A, B := SmithForm(M);
    // Extract nontrivial diagonal elements of S and their positions
    l := [[S[i,i], i] : i in [1..NumberOfColumns(S)] | S[i,i] notin {0,1}];
    // Submatrix of B corresponding to nontrivial elements
    sA := Matrix(Rationals(), ColumnSubmatrixRange(B, l[1][2], l[#l][2]));
    // Scale columns by the reciprocal of the diagonal entries
    for i in [1..#l] do
        MultiplyColumn(~sA, 1 / l[i][1], i);
    end for;
    // Compute the quadratic form in the scaled basis
    Q := Transpose(sA) * Matrix(Rationals(), M) * sA;
    // Reduce entries modulo 2 for the quadratic form
    for i, j in [1..NumberOfColumns(Q)] do
        if i ne j then
            Q[i, j] := Q[i, j] - Floor(Q[i, j]);
        else
            Q[i, j] := Q[i, j] - Floor(Q[i, j]) + (Floor(Q[i, j]) mod 2);
        end if;
    end for;
    // Return the discriminant group and the quadratic form
    return [l[i][1] : i in [1..#l]], Q;
end function;

// MatrixMod2
// Input: A matrix with rational coefficients
// Output: The matrix with its entries reduced modulo 2

MatrixMod2 := function(Q);
    for i, j in [1..Nrows(Q)] do
        if i ne j then 
            Q[i, j] := Q[i, j] - Floor(Q[i, j]);
        else 
            Q[i, j] := Q[i, j] - 2 * Floor(Q[i, j] / 2);
        end if;
    end for;
    return Q;
end function;

// IsotropicElements
// Input: A matrix with integer coefficients
// Output: The null elements of its discriminant group

IsotropicElements := function(M);
    // Compute the discriminant group and quadratic form
    v, U := DiscriminantGroup(M);
    Q := Rationals();
    // Define an abelian group with orders from the discriminant group
    A := AbelianGroup(v);
    // Return elements that are isotropic with respect to the quadratic form
    return [Eltseq(a) : a in A | 
            MatrixMod2(Matrix(Q, 1, #v, Eltseq(a)) * U * Matrix(Q, #v, 1, Eltseq(a)))[1, 1] eq 0];
end function;

// CompareDiscriminants
// Input: Two matrices M and Q with integer coefficients
// Output: True if the discriminant group and quadratic form of M and Q are equivalent, false otherwise

CompareDiscriminants := function(M, Q)
    // Compute the discriminant group and quadratic form for M
    v, U := DiscriminantGroup(M);
    // Compute the discriminant group and quadratic form for Q
    w, D := DiscriminantGroup(Q);

    // If the discriminant groups have different orders, return false
    if v ne w then 
        return false;
    end if;

    // Create an abelian group from the discriminant group
    A := AbelianGroup(v);
    // Compute the automorphism group of A
    Aut := AutomorphismGroup(A);
    // Obtain the permutation representation of the automorphism group
    f, G := PermutationRepresentation(Aut);
    h := Inverse(f); // Inverse map for automorphisms

    // Compute the matrices representing automorphisms
    ll := [Matrix(Rationals(), [Eltseq((h(g))(A.i)) : i in [1..Ngens(A)]]) : g in G];
    // Apply each automorphism to the quadratic form and reduce modulo 2
    dd := [MatrixMod2(a * U * Transpose(a)) : a in ll];

    // Check if the quadratic form of Q is equivalent to any transformed form of M
    return D in dd;
end function;
