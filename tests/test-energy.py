#!/usr/bin/env python

def compare_MD(reference, current, reltol):
    """Compare MD energies.

    Given a reference output and the current output, compare the MD
    energies to within the relative tolerance given by reltol.

    """

    import sys
    #energy = re.compile("MD_data\s+([0-9eE.+-]+)\s+([0-9eE.+-]+)")

    fd = open(reference)
    reference_energies = []
    for line in fd:
        result = line.split()
        reference_energies.append(float(result[0]))
    fd.close()
    
    fd = open(current)
    current_energies = []
    for line in fd:
        result = line.split()
        current_energies.append(float(result[0]))
    fd.close()

 
    if len(reference_energies) != len(current_energies):
        raise Exception("[error] different number of MD steps\n"
                        + ("  reference ran for %4d steps\n" % (len(reference_energies)))
                        + ("  current ran for   %4d steps\n" % (len(current_energies)))
                        + "  can not compare")

    result = True
    for i in range(len(reference_energies)):
        diff = abs(reference_energies[i] - current_energies[i])
        if reference_energies[i] != 0:
            diff = abs(diff/reference_energies[i])
        if diff > reltol:
            print("failure in MD step %d" % (i+1))
            result = False
    if not result:
        raise Exception(("[error] when comparing '%s' with '%s'" % (reference, current))
                        + "energies do not agree")

    print("Energy test passed without failure ...")  

def main():
    """The main function.
    """

    import argparse, os, sys
    
    parser = argparse.ArgumentParser(description="""Script to compare MD results by using the total energy""")
    parser.add_argument("--reference",
                        help="The reference output")
    parser.add_argument("--current",
                        help="The current output")
    parser.add_argument("--reltol",
                        help="Relative tolerance when comparing, default is %(default)s",
                        type=float,
                        default=1e-10)
    options = parser.parse_args()

    compare_MD(options.reference, options.current, options.reltol)

if __name__ == "__main__":
    main()
