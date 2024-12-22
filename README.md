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

## Usage

These functions are designed to work together to compute properties of discriminant groups derived from integer matrices. The workflow typically involves:
1. Using `DiscriminantGroup` to obtain the discriminant group and quadratic form.
2. Utilizing `IsotropicElements` to find isotropic vectors in the group.

## Example

```magma
// Example matrix
M := Matrix(Integers(), [[2, 1], [1, 2]]);

// Compute discriminant group and quadratic form
DG, QF := DiscriminantGroup(M);
print "Discriminant Group:", DG;
print "Quadratic Form:", QF;

// Find isotropic elements
iso := IsotropicElements(M);
print "Isotropic Elements:", iso;
