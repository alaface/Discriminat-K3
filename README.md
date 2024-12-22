# Discriminant Group Computation

This repository contains functions for computing the discriminant group of a matrix, along with associated operations, implemented in the Magma language.

## Functions

### `DiscriminantGroup(M)`
- **Input:** A matrix with integer coefficients.
- **Output:** 
  1. The discriminant group as a list of integers.
  2. The quadratic form in $ \mathbb{Q}/2\mathbb{Z} $.

This function calculates the discriminant group of a given integer matrix using the Smith Normal Form, and computes the quadratic form in a reduced basis.

### `MatrixMod2(Q)`
- **Input:** A matrix with rational coefficients.
- **Output:** The same matrix with all entries reduced modulo 2.

This utility function reduces the entries of a matrix modulo 2, ensuring compatibility with quadratic form computations.

### `IsotropicElements(M)`
- **Input:** A matrix with integer coefficients.
- **Output:** The list of null elements of its discriminant group.

This function identifies the isotropic elements of the discriminant group by testing quadratic form values for nullity modulo 2.

### `CompareDiscriminants(M, Q)`
- **Input:** 
  - `M`: A matrix with integer coefficients.
  - `Q`: Another matrix with integer coefficients.
- **Output:** `true` if the discriminant group and quadratic form of `M` and `Q` are equivalent, `false` otherwise.

This function compares the discriminant groups and quadratic forms of two matrices. It determines whether there exists an automorphism of the discriminant group that transforms the quadratic form of `M` into that of `Q`.

## Usage

These functions are designed to work together to compute properties of discriminant groups derived from integer matrices. The workflow typically involves:
1. Using `DiscriminantGroup` to obtain the discriminant group and quadratic form.
2. Utilizing `IsotropicElements` to find isotropic vectors in the group.
3. Using `CompareDiscriminants` to check equivalence between two matrices.

## Example

```magma
// Example matrices
M := Matrix(Integers(), [[0, 1, 0], [1, 0, 0], [0, 0, -12]]);
Q := Matrix(Integers(), [[-2, 1, 0], [1, -2, 0], [0, 0, 4]]);

// Compare discriminants and quadratic forms
isEquivalent := CompareDiscriminants(M, Q);
if isEquivalent then
    print "The discriminant groups and quadratic forms are equivalent.";
else
    print "The discriminant groups and quadratic forms are not equivalent.";
end if;
