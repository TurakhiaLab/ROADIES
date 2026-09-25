# Supported Model Scope

This is a current implementation constraint, not a list of models that can be
enabled without additional likelihood and validation work. The public CLI in
the root `README.md` is authoritative.

## Maintained model

- alphabet: 4-state DNA
- substitution model: GTR
- rate heterogeneity: 1--8 discrete-Gamma categories
- category weights: equal, `1 / ncat`
- invariant-site proportion: `pinv = 0`

## Rate heterogeneity

- Production support is intentionally limited to discrete-Gamma categories
  with equal weights.
- FreeRate models, mixture models, and custom rate-category weights are
  rejected by the public parsing paths.
- Keep the internal category-weight arrays: likelihood kernels need the
  uniform `1 / ncat` factors when summing discrete-Gamma categories.

## Unsupported model shapes

The public validation paths reject amino-acid models, non-GTR substitution
models, invariant-site likelihoods, FreeRate, mixture models, and custom
rate-category weights. Supporting one of these requires coordinated changes to
input validation, model installation, PMAT construction, likelihood kernels,
optimization, and regression coverage; widening a vector or removing a parser
check is not sufficient.
