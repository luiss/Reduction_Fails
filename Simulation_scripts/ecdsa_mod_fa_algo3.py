import math
from random import randrange, seed
from ecpy.curves import Curve
from sympy.ntheory.modular import crt
from collections import namedtuple
import argparse

from ilp_solver import build_factor_lists, solve_for_pfa
from ecdl_tool import crack_baby_giant

STUCK_00 = "00"
STUCK_FF = "FF"
BIT_FLIP = 2

Entry = namedtuple("Entry", ("residue", "modulus"))

def faulty_signature(e, d, curve, fault_model, index):
    # Generate faulty modulus
    if fault_model == STUCK_FF:
        n = curve.order | (255 << (8 * index))
    elif fault_model == STUCK_00:
        n = curve.order - (curve.order & (255 << (8 * index)))
    elif fault_model == BIT_FLIP:
        n = curve.order ^ (1 << index)
    else:
        raise ValueError("Invalid fault model")

    # Ensure k invertible mod n
    k = randrange(curve.order)

    R = k * curve.generator
    r = int(R.x) % n
    s = (pow(k, -1, n) * (e + r * d)) % n
    # Countermeasures test
    # i = pow(d, -1, curve.order)
    # ds = pow(d, 2, curve.order)
    # s = (pow(k, -1, n) * (e + i * r * ds)) % n

    return r, s, n, index, k

def ILP_solver_wrapper(curve, B, fault_model, tau):
    """
    Returns:
        L: indices to fault
        factors: dict index -> {p: max exponent}
    """
    lists = build_factor_lists(curve, fault_model, B)

    if not lists:
        return [], {}

    res = solve_for_pfa(lists, tau, Delta=10)
    if res is None:
        return [], {}

    selected_indices, expo, entropy, faults = res

    factors = {}
    for i in selected_indices:
        factor_list = list()
        for k in lists[i].keys():
            factor = k**lists[i][k]
            if all([factor not in factors[key] for key in factors.keys()]):
                factor_list.append(factor)
        factors[i] = factor_list



    return selected_indices, factors, entropy, faults


def recover_mod_pk(r, s, e, nf, p):
    """
    Recover d mod p^k by lifting (same logic as ecdsa_mod_fa.py)
    """
    if r % p == 0:
        return None
        
    if s % p != 0:
        return None

    residue = (-e * pow(r, -1, p)) % p
    return Entry(residue, p)

def solve_ecdl(curve, Q, entries):
    """
    Combine residues and solve discrete log
    """
    moduli = [e.modulus for e in entries.values()]
    residues = [e.residue for e in entries.values()]

    partial, _ = crt(moduli, residues)
    modulus = math.prod(moduli)

    if modulus >= curve.order:
        return int(partial % curve.order)

    # fallback BSGS
    base = modulus * curve.generator
    target = Q - partial * curve.generator

    x = crack_baby_giant(base, math.ceil(curve.order / modulus), target)
    if x is None:
        return None

    return int(x * modulus + partial)

def optimized_fault_attack_full(curve_name, B, tau, fault_model):
    curve = Curve.get_curve(curve_name)

    n = curve.order
    d = randrange(n)
    Q = d * curve.generator

    print(f"[+] Target key: {d:X}")

    # --- ILP phase ---
    L, factors, entropy, faults = ILP_solver_wrapper(curve, B, fault_model, tau)

    if not L:
        print("[-] ILP failed")
        return None

    print(f"[+] Selected indices: {L}")
    print(f"[+] Expecting: {entropy} bits with {faults} faulty signatures")

    # entries[p] = best Entry (max exponent)
    entries = {}
    faults = 0
    faults_with_error = 0

    for i in L:
        print(f"[+] Processing fault index {i}")
        found_factors = list()

        while len(factors[i]) != len(found_factors):
            try:
                e = randrange(n)

                r, s, nf, _, _ = faulty_signature(e, d, curve, fault_model, i)
                faults += 1
                for p in factors[i]:
                    entry = recover_mod_pk(r, s, e, nf, p)

                    if entry is None:
                        continue

                    if d % p != entry.residue:
                        print("error")

                    if p not in entries:
                        entries[p] = entry
                        print(f"    [+] p={p}")
                        found_factors.append(p)

            # Can't invert ephemeral key with faulty modulus
            except ValueError:
                faults_with_error += 1
                continue

    # --- entropy check ---
    if entries:
        current_prod = math.prod(e.modulus for e in entries.values())
        entropy = math.log2(current_prod)

        if entropy > tau:
            print(f"[+] Entropy reached: {entropy:.2f} > {tau}")

            recovered = solve_ecdl(curve, Q, entries)

            if recovered is None:
                print("[-] ECDL failed")
                return None

            print(f"[+] Recovered key: {recovered:X}")
            print(f"[+] Success: {recovered == d}")
            print(f"[+] Total number of signatures: {faults + faults_with_error} ({faults} used)")


            return recovered
    print("[-] Attack failed")
    return None


def main():
    parser = argparse.ArgumentParser(description="ECDSA Modulo FA Simulation")
    parser.add_argument("--curve", default="secp256r1", help="Curve name (default secp256r1)")
    parser.add_argument("--fault-model", choices=["00", "FF", "ext"], required=True,
                        help="Fault model: 00 (STUCK_00) or FF (STUCK_FF)")
    parser.add_argument("--B", type=int, required=True,
                        help="Factor bound (e.g., 4096)")
    parser.add_argument("--tau", type=int, default=216,
                        help="Target entropy (default: 216)")
    parser.add_argument("--faulty-moduli", type=str, default=None,
                        help="Path to JSON list of faulty moduli")

    args = parser.parse_args()
    fault_model = STUCK_00 if args.fault_model == "00" else STUCK_FF

    optimized_fault_attack_full(
        curve_name=args.curve,
        B=args.B,
        tau=args.tau,
        fault_model=fault_model
    )

if __name__ == "__main__":
    seed(1)
    main()
