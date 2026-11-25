# Dummy model function for the micro-dumux surrogate.
# We do not wrap the original DuMuX model because we will directly provide
# snapshots (computed by the Micro Manager) to BayesValidRox.
def model(samples):
    return None
