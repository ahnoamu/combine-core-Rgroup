# Combine Core and Rgroup

Combine Core and Rgroup with RDKit

## Example

Given a table of the core molecule, associated Rgroups, and reported activity values for QSAR as shown below, the core and Rgroup can be combined as follows:

### Example Table: 
| Compound (Core)                                        | R1      | R2
|--------------------------------------------------------|---------|-------------|
| CCN1C(=CC2=CC=CC=[N+]2CC)SC(=C3N(C4=CC=CC=C4S3)C)C1=O  | -CH2NH  | -CH2CH=CH2  |
|                                                        | -CH2CH3 | -CH2NCH3    |
|                                                        | -CH2NH  | -CH2CH=N    |

### Example use

* Prepare core and Rgroup as shown in [`combineCoreRgroup_procedure.ipynb`](combineCoreRgroup_procedure.ipynb).

* Weld prepared molecules using either methods.

* Method 1
```shell
python3 combineCoreRgroup.py  -c  "[*:1]N1C(=CC2=CC=CC=[N+]2CC)SC(=C3N(C4=CC=CC=C4S3)C)C1=O"  -r  "[*:1]CN"  -m   m1
```

* Method 2
```shell
python3 combineCoreRgroup.py  -c  "[*:2]N1C(=CC2=CC=CC=[N+]2CC)SC(=C3N(C4=CC=CC=C4S3)C)C1=O"  -r  "[*:1]CN"  -m   m2
```

N/B: Method 1: Both molecules MUST have the same pattern of dummy atoms
