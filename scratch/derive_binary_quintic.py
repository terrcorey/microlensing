"""Symbolic derivation of the binary-lens quintic polynomial coefficients
used by `lc_models._quintic_coefficients`.

Not part of the regular pipeline and not a project dependency -- lives in
scratch/ with the project's other one-time/dev scripts, run manually (in a
throwaway venv with `pip install sympy`) only if you need to regenerate or
audit those coefficients.

Lens equation (Witt 1990 form), with z1/z2 real (center-of-mass frame,
lenses on the real axis) so conj(z1) = z1, conj(z2) = z2:

    zeta = z - m1/conj(z - z1) - m2/conj(z - z2)     ... (1)

Conjugating (1):

    conj(zeta) = conj(z) - m1/(z - z1) - m2/(z - z2)  ... (2)

Solving (2) for conj(z) and substituting into (1), then clearing
denominators, gives a degree-5 polynomial in z. Its roots include every
true image position plus spurious extras introduced by clearing the
conjugate -- `binary_images` filters those out by checking each root
against the original equation (1).
"""
import sympy as sp

z, zeta, zetab, z1, z2, m1 = sp.symbols("z zeta zetab z1 z2 m1")
m2 = 1 - m1

zb_expr = zetab + m1 / (z - z1) + m2 / (z - z2)  # conj(z) from (2)
eq = z - m1 / (zb_expr - z1) - m2 / (zb_expr - z2) - zeta  # (1), rearranged to = 0

num, den = sp.fraction(sp.together(eq))
poly = sp.Poly(sp.expand(num), z)
assert poly.degree() == 5
coeffs = poly.all_coeffs()  # highest degree first

replacements, reduced = sp.cse(coeffs, optimizations="basic")

print("def _quintic_coefficients(zeta, zetab, m1, z1, z2):")
for sym, expr in replacements:
    print(f"    {sym} = {sp.printing.pycode(expr)}")
print("    return [")
for expr in reduced:
    print(f"        {sp.printing.pycode(expr)},")
print("    ]")
