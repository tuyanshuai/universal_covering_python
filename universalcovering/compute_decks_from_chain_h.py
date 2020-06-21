from universalcovering.hyperbolic_deck_transform import *


def compute_decks_from_chain_h(z,chain,mc):
    ns = len(chain)
    decks = [list()]*ns
    for i in range(ns):
        ch1 = chain[i]
        ch2 = chain[mc[i]]
        deck = hyperbolic_deck_transform(z[ch1[0]],z[ch1[-1]],z[ch2[-1]],z[ch2[0]])
        decks[i] = deck

    return decks
