DLtab
-----

Python script to convert DFTB+'s skf format into LATTE's table format.

Usage
-----

Put skf files into a directory (`<skf_DIR>`), then call the python script from terminal:

```bash
python DLtab.py <skf_DIR> <output_DIR>
```

Then three files will be generated into the `<output_DIR>` folder, `electrons.dat`, `bondints.table`, and `ppots.dftb`. The files are ready to used for non spin-polarized calculations with LATTE. For spin-polarized calculations, corresponding spin paramerers for each elements, "Wss", "Wpp", "Wdd", and "Wff", are required in the `electrons.dat` file.

To carry out LATTE calculations with the converted parameters, set following flags in the LATTE input file, `latte.in`:

```latte
PARAMPATH= "<output_DIR>"
SCLTYPE= TABLE
PPOTON= 2
```

Notes
-----

 * The script uses Hubbard U for the s shell as U parameter for the element, which is the same as the default for DFTB+ code.
 * To call LATTE from LAMMPS, make sure masses of elements in `electrons.dat` match values in lammps geometry files.
