from universalcovering.hyperbolic_deck_transform import *
def  compute_decks_from_bp_h(bp, mc):

    ns = len(bp)
    decks = [list()]*ns
    for i in range(ns):
        i2 = i + 1
        if i2 >= ns:
            i2 = i2 - ns

        j = mc[i]
        j2 = j + 1
        if j2 >= ns:
            j2 = j2 - ns

        deck = hyperbolic_deck_transform(bp[i], bp[i2], bp[j2], bp[j])
        decks[i] = deck
    return decks
