# creates a color that reproduces transparency
def fakealpha(rgb, alpha):
    return [alpha * c1 + (1 - alpha) * c2
            for (c1, c2) in zip(rgb, (1, 1, 1))]
