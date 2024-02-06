# build_pdb

Simple command-line awk programs for building and manipulating macromolecular structures as PDB files

## Motivation

There are some very powerful model-building tools out there, but sometimes you just want to automate something simple.

## Description

These self-contained awk programs build and modify  take simple keyword-value directives from their standard input, preceding the text of a Protein Data Bank format file, and execute changes or additions to the structure as they come up. These text commands look like this:
```
BUILD MET A51 180 180 180
```
Which means to build a methionine side chain at residue 51 of chain A and give it chi angles of all 180 deg.
 The standard output is a new copy of the input with the edits in place. If you specify `ONLYNEW`, then just the modified atoms are output. The keyword-value directives are not passed through.<br>
 You can also build off the termini of the molecule by specifying phi, psi and omega angles.

## Getting Started

### Dependencies

These are awk scripts, so you will need a working `/bin/awk`.  If your OS has awk installed somewhere else, you will need to either edit the shebang (first line of the script), or run them like this:
```
awk -f build_side.awk
```
There are no other dependencies.

### Installing

* git clone this repo
* copy all the files into somewhere in your shell `$PATH`, and make them executable:
```
    chmod u+x build_*.awk
```
The next thing you might want to do is read about how unix command-line pipes work, and familiarize yourself with concepts like standard input (stdin) and standard output (stdout).

### Executing program

#### peptide synthesis
The first thing you might want do do is try a one-liner, like this:
```
echo "BUILD ALA\nBUILD ALA\nBUILD ALA" | ./build_n2c.awk | tee tripeptide.pdb
```
and you should see:
```
ATOM      1  N   ALA     1      -1.501  -0.000   2.251  1.00 20.00
ATOM      2  CA  ALA     1      -0.467  -0.000   3.268  1.00 20.00
ATOM      3  C   ALA     1      -1.085  -0.000   4.656  1.00 20.00
ATOM      4  O   ALA     1      -2.305  -0.000   4.809  1.00 20.00
ATOM      5  CB  ALA     1       0.412   1.233   3.133  1.00 20.00
ATOM      6  N   ALA     2      -0.222  -0.000   5.669  1.00 20.00
ATOM      7  CA  ALA     2      -0.663  -0.000   7.050  1.00 20.00
ATOM      8  C   ALA     2       0.531  -0.000   7.991  1.00 20.00
ATOM      9  O   ALA     2       1.682   0.000   7.557  1.00 20.00
ATOM     10  CB  ALA     2      -1.503  -1.233   7.341  1.00 20.00
ATOM     11  N   ALA     3       0.241  -0.000   9.289  1.00 20.00
ATOM     12  CA  ALA     3       1.274   0.000  10.306  1.00 20.00
ATOM     13  C   ALA     3       0.656  -0.000  11.694  1.00 20.00
ATOM     14  O   ALA     3      -0.564  -0.000  11.847  1.00 20.00
ATOM     15  CB  ALA     3       2.153   1.233  10.171  1.00 20.00
ATOM     16  OXT ALA     3       1.519  -0.000  12.707  1.00 20.00
```
That vertical bar character `|` is a unix pipe. This tells the shell to take the text that would normally have gone to the screen from the `echo` command and re-direct it as though you typed it yourself on the keyboard into the running `build_n2c.awk` program. The second pipe `|` also takes what would have been printed onto the screen by `build_n2c.awk` and re-directs it into the next program, which is `tee`. This is analogous to a t-junction in a water pipe, where the input is split in two: copied to the screen (stdout), as well as written to the file `tripeptide.pdb`. 
##### debugging
If you only see one ALA, its probably because your `echo` command doesn't know what `\n` means. It is supposed to denote a newline. Try this instead:
```
/bin/echo -e "BUILD ALA\nBUILD ALA\nBUILD ALA" | ./build_n2c.awk | tee tripeptide.pdb
```
If you get an error:
```
./build_n2c.awk: Permission denied.
```
Then you need to run this:
```
chmod a+x ./build_*.awk
```
If you get this error:
```
./build_n2c.awk: Command not found.
```
then your system's awk program is not located in `/bin/awk`.  Try running this:
```
which awk
```
this should tell you where your OS has its copy of `awk`. If, for example, it says `/usr/local/bin/awk`, then you need to edit the first line of the awk program. This line is called the "shebang" because it begins with a hash and a bang `#!`. This is a signal to your shell telling it which program to use to execute the script. So, instead of beginning with `#! /bin/awk -f`, you might want to change it to `#! /usr/local/bin/awk -f`. Alternately, you can run the program like this:
```
echo "BUILD ALA\nBUILD ALA\nBUILD ALA" | awk -f ./build_n2c.awk | tee tripeptide.pdb
```
If the response to `which awk` is `Command not found.`, then you will need to install `awk`.<br>
If you get the error message:
```
tee: Command not found.
``` 
Then you probably have a Mac. So, instead of enjoying the convenience of a pipe tee, you will want to do a two-line command:
```
/bin/echo -e "BUILD ALA\nBUILD ALA\nBUILD ALA" | awk -f ./build_n2c.awk > tripeptide.pdb ; cat tripeptide.pdb
```
that command will run on a Mac.<br>
This may seem like a lot of problems, but it is also all of the problems you might encounter on the various unix distros that are out there. If you have a modern and full-featured Linux installed, then the first command above will work.

#### mutagenesis
Now that you have an alanine tripeptide, lets change one of the residues to methionine:
```
echo "BUILD MET 2 180 180 180" | cat - tripeptide.pdb | ./build_side.awk | tee mutated.pdb
```
and you should see:
```
ATOM      1  N   ALA     1      -1.501  -0.000   2.251  1.00 20.00
ATOM      2  CA  ALA     1      -0.467  -0.000   3.268  1.00 20.00
ATOM      3  C   ALA     1      -1.085  -0.000   4.656  1.00 20.00
ATOM      4  O   ALA     1      -2.305  -0.000   4.809  1.00 20.00
ATOM      5  CB  ALA     1       0.412   1.233   3.133  1.00 20.00
ATOM      1  N   MET     2      -0.222   0.000   5.669  1.00 20.00
ATOM      2  CA  MET     2      -0.663   0.000   7.050  1.00 20.00
ATOM      3  C   MET     2       0.531   0.000   7.991  1.00 20.00
ATOM      4  O   MET     2       1.682   0.000   7.557  1.00 20.00
ATOM      5  CB  MET     2      -1.503  -1.233   7.341  1.00 20.00
ATOM      6  CG  MET     2      -1.959  -1.210   8.811  1.00 20.00
ATOM      7  SD  MET     2      -2.953  -2.618   9.254  1.00 20.00
ATOM      8  CE  MET     2      -3.427  -2.497  10.965  1.00 20.00
ATOM     11  N   ALA     3       0.241  -0.000   9.289  1.00 20.00
ATOM     12  CA  ALA     3       1.274   0.000  10.306  1.00 20.00
ATOM     13  C   ALA     3       0.656  -0.000  11.694  1.00 20.00
ATOM     14  O   ALA     3      -0.564  -0.000  11.847  1.00 20.00
ATOM     15  CB  ALA     3       2.153   1.233  10.171  1.00 20.00
ATOM     16  OXT ALA     3       1.519  -0.000  12.707  1.00 20.00
END
```
The three chi angles are all 180 deg in this case. That `cat - tripeptide.pdb` command means to take the text that comes in from stdin and concatenate it with the `tripeptide.pdb` file, putting the stdin text first.<br>
 If you leave out any of the chi angles then the atom at the end of that dihedral won't be built.  For example, `BUILD LYS 2 60 60` will build a lysine as a norvaline (but call it `LYS`). If you don't care about dihedrals and just want the most popular rotamer, replace the numbers with question marks:
```
echo "BUILD MET 2 ? ? ?" | cat - tripeptide.pdb | ./build_side.awk | tee mutated.pdb
```
which will produce:
```
ATOM      1  N   ALA     1      -1.501  -0.000   2.251  1.00 20.00
ATOM      2  CA  ALA     1      -0.467  -0.000   3.268  1.00 20.00
ATOM      3  C   ALA     1      -1.085  -0.000   4.656  1.00 20.00
ATOM      4  O   ALA     1      -2.305  -0.000   4.809  1.00 20.00
ATOM      5  CB  ALA     1       0.412   1.233   3.133  1.00 20.00
ATOM      1  N   MET     2      -0.222   0.000   5.669  1.00 20.00
ATOM      2  CA  MET     2      -0.663   0.000   7.050  1.00 20.00
ATOM      3  C   MET     2       0.531   0.000   7.991  1.00 20.00
ATOM      4  O   MET     2       1.682   0.000   7.557  1.00 20.00
ATOM      5  CB  MET     2      -1.503  -1.233   7.341  1.00 20.00
ATOM      6  CG  MET     2      -2.794  -1.184   6.503  1.00 20.00
ATOM      7  SD  MET     2      -3.908   0.106   7.015  1.00 20.00
ATOM      8  CE  MET     2      -4.723  -0.375   8.522  1.00 20.00
ATOM     11  N   ALA     3       0.241  -0.000   9.289  1.00 20.00
ATOM     12  CA  ALA     3       1.274   0.000  10.306  1.00 20.00
ATOM     13  C   ALA     3       0.656  -0.000  11.694  1.00 20.00
ATOM     14  O   ALA     3      -0.564  -0.000  11.847  1.00 20.00
ATOM     15  CB  ALA     3       2.153   1.233  10.171  1.00 20.00
ATOM     16  OXT ALA     3       1.519  -0.000  12.707  1.00 20.00
END
```
Which is the `- - -` (minus minus minus) or chi1=-64.5, chi2=-68.5, chi3=-75.6 rotamer for methionine. the most common one reported by Ponder & Richards. You can also use `+`, `t`, or `0` to represent the plus, trans, and zero rotamer angles.<br>
Normally, this script passes-thru all atoms not involved in rebuilding. If you want only the rebuilt residue to be output, put "ONLYNEW" as the first word on a line at the top of the PDB file.

#### adding on to the N terminus
The N terminus can often be disordered, so I wrote a separate program for tacking atoms onto it
```
cat << EOF | cat - mutated.pdb | ./build_c2n.awk | tee addN.pdb
BUILD MET   0 -40
BUILD ALA -64 -40
BUILD PRO -60 -40
BUILD ALA -64 -40
BUILD -64
EOF
```
and you should see:
```
ATOM     16  N   MET    -3      -2.965   3.248  -3.033  1.00 20.00
ATOM     17  CA  MET    -3      -1.774   2.557  -2.578  1.00 20.00
ATOM     18  C   MET    -3      -1.604   2.720  -1.076  1.00 20.00
ATOM     19  O   MET    -3      -1.224   1.786  -0.373  1.00 20.00
ATOM     20  CB  MET    -3      -0.539   3.117  -3.264  1.00 20.00
ATOM     11  N   ALA    -2      -1.891   3.925  -0.593  1.00 20.00
ATOM     12  CA  ALA    -2      -1.777   4.232   0.820  1.00 20.00
ATOM     13  C   ALA    -2      -2.759   3.400   1.628  1.00 20.00
ATOM     14  O   ALA    -2      -2.442   2.913   2.711  1.00 20.00
ATOM     15  CB  ALA    -2      -2.072   5.701   1.072  1.00 20.00
ATOM      6  N   PRO    -1      -3.963   3.242   1.083  1.00 20.00
ATOM      7  CA  PRO    -1      -5.005   2.474   1.736  1.00 20.00
ATOM      8  C   PRO    -1      -4.556   1.037   1.947  1.00 20.00
ATOM      9  O   PRO    -1      -4.821   0.432   2.984  1.00 20.00
ATOM     10  CB  PRO    -1      -6.268   2.463   0.890  1.00 20.00
ATOM      1  N   ALA     0      -3.869   0.497   0.944  1.00 20.00
ATOM      2  CA  ALA     0      -3.375  -0.865   1.000  1.00 20.00
ATOM      3  C   ALA     0      -2.349  -1.015   2.111  1.00 20.00
ATOM      4  O   ALA     0      -2.323  -2.018   2.823  1.00 20.00
ATOM      5  CB  ALA     0      -2.716  -1.247  -0.316  1.00 20.00
ATOM      1  N   ALA     1      -1.501  -0.000   2.251  1.00 20.00
ATOM      2  CA  ALA     1      -0.467  -0.000   3.268  1.00 20.00
ATOM      3  C   ALA     1      -1.085  -0.000   4.656  1.00 20.00
ATOM      4  O   ALA     1      -2.305  -0.000   4.809  1.00 20.00
ATOM      5  CB  ALA     1       0.412   1.233   3.133  1.00 20.00
ATOM      1  N   MET     2      -0.222   0.000   5.669  1.00 20.00
ATOM      2  CA  MET     2      -0.663   0.000   7.050  1.00 20.00
ATOM      3  C   MET     2       0.531   0.000   7.991  1.00 20.00
ATOM      4  O   MET     2       1.682   0.000   7.557  1.00 20.00
ATOM      5  CB  MET     2      -1.503  -1.233   7.341  1.00 20.00
ATOM      6  CG  MET     2      -2.794  -1.184   6.503  1.00 20.00
ATOM      7  SD  MET     2      -3.908   0.106   7.015  1.00 20.00
ATOM      8  CE  MET     2      -4.723  -0.375   8.522  1.00 20.00
ATOM     11  N   ALA     3       0.241  -0.000   9.289  1.00 20.00
ATOM     12  CA  ALA     3       1.274   0.000  10.306  1.00 20.00
ATOM     13  C   ALA     3       0.656  -0.000  11.694  1.00 20.00
ATOM     14  O   ALA     3      -0.564  -0.000  11.847  1.00 20.00
ATOM     15  CB  ALA     3       2.153   1.233  10.171  1.00 20.00
ATOM     16  OXT ALA     3       1.519  -0.000  12.707  1.00 20.00
```
You may want to re-number residues after this.  Note that an empty `BUILD` command comes last because the psi angle of the previous N-terminal residue is not defined. You therefore must specify it.

#### adding an atom in general
If you just want to add an atom with some known bond length, angle and dihedral relative to existing atoms, then this script will do that. Formatting is like this:
```
BUILD A23:N A23:CA A23:C A23:OXT 180
```
Which will add an atom called `OXT` with a dihedral angle of 180 relative to the N-CA-C sequence of atoms of residue 23 in chain A. This is normally where the `OXT` atom should go. Default bond length for the new bond is 1.54 A, the default new bond angle is 109.5 deg, and the default new dihedral angle is 0 deg. You can also specify occupancy and B factor for the new atom like this:
```
BUILD A23:N A23:CA A23:C A23:OXT 180 BOND 1.23 ANGLE 120 CONF A BFAC 15 OCC 0.5
```
Hope that is useful!


## Help

Let me know if you have any questions or find any bugs.

