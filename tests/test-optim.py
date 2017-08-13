#!/usr/bin/env python

def compare_coordinates(reference, current, reltol):
    """Compare coordinates.

    Given a reference output and the current output, compare the
    coordinates to within the relative tolerance given by reltol.

    """

    import re
    position = re.compile("^\s*([a-zA-Z]+)\s+([0-9eE.+-]+)\s+([0-9eE.+-]+)\s+([0-9eE.+-]+)")

    fd = open(reference)
    lines = fd.readlines()
    fd.close()

    reference_positions = []
    i = 0
    while i < len(lines):
        N = int(lines[i])
        i += 2
        reference_positions.append([])
        for j in range(i, i+N):
            result = position.search(lines[j])
            reference_positions[-1].append({"name": result.group(1),
                                            "position": [float(result.group(2)),
                                                         float(result.group(3)),
                                                         float(result.group(4))]})
        i += N

    fd = open(current)
    lines = fd.readlines()
    fd.close()

    current_positions = []
    i = 0
    while i < len(lines):
        N = int(lines[i])
        i += 2
        current_positions.append([])
        for j in range(i, i+N):
            result = position.search(lines[j])
            current_positions[-1].append({"name": result.group(1),
                                          "position": [float(result.group(2)),
                                                       float(result.group(3)),
                                                       float(result.group(4))]})
        i += N

    if len(reference_positions) != len(current_positions):
        raise Exception("[error] different number of optimization steps\n"
                        + ("  reference ran for %4d steps\n" % (len(reference_positions)))
                        + ("  current ran for   %4d steps\n" % (len(current_positions)))
                        + "  can not compare")

    import math
    result = True
    for i in range(len(reference_positions)):
        rmsd = 0
        for j in range(len(reference_positions[i])):
            ref = reference_positions[i][j]["position"]
            cur = current_positions[i][j]["position"]
            rmsd += (ref[0] - cur[0])**2 + (ref[1] - cur[1])**2 + (ref[2] - cur[2])**2
        rmsd = math.sqrt(rmsd/len(reference_positions[i]))
        if rmsd > reltol:
            print("failure in optimization step %d" % (i+1))
            print("rmsd = %e" % (rmsd))
            result = False
    if not result:
        raise Exception(("[error] when comparing '%s' with '%s'" % (reference, current))
                        + "structures do not agree")
        
    print("optim test passed without failure ...")  

def main():
    """The main function.
    """
    
    import argparse, os, sys
      
    parser = argparse.ArgumentParser(description="""Script to get compare two optimization results""")
    parser.add_argument("--reference",
                        help="The reference output")
    parser.add_argument("--current",
                        help="The current output")
    parser.add_argument("--reltol",
                        help="Relative tolerance when comparing, default is %(default)s",
                        type=float,
                        default=1e-10)
    options = parser.parse_args()

    compare_coordinates(options.reference, options.current, options.reltol)

if __name__ == "__main__":
    main()
