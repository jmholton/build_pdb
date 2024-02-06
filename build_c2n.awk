#! /bin/awk -f
#
#
#	build a polypeptide chain with given phi-psi angles
#	in C-N direction
#
# specify what residues to tack onto the N terminus like this:
# BUILD RES phi psi
# EG:
# BUILD MET   0 -40
# BUILD ALA -64 -40
# BUILD PRO -60 -40
# BUILD ALA -64 -40
# BUILD -64
# will add the above three residues (with the phi, psi provided)
# to the N-terminus of your pdb.  Residue numbers and chain IDs
# will be the same as for the first complete residue in your pdb.
#
# also, because the phi value of the N-terminal template residue
# is undefined, you must specify it with one more BUILD line:
# BUILD lastphi
#
#
BEGIN{
    if(!resnum) resnum = 1000
    if(!next_phi) next_phi = 180

    # initialize "phantom" atoms to start phi-psi from
    N["X"]=0
    N["Y"]=0
    N["Z"]=0

    CA["X"]=-2.42
    CA["Y"]=0
    CA["Z"]=0

    C["X"]=-1.21
    C["Y"]=0
    C["Z"]=0.953

    O["X"]=0
    O["Y"]=0
    O["Z"]=0
}

/^CRYST/
/^REMARK/
/^SCALE/

/^ATOM/{
    # print these out later
    ++oldatoms
    OldAtom[oldatoms] = $0

    resnum = substr($0, 23, 4)+0
    restyp = substr($0, 18, 3)
    chain  = substr($0, 22, 1)          # O/Brookhaven-style segment ID
    split(substr($0, 13, 4), a)
    atom   = a[1];
    x      = substr($0, 31, 8)+0
    y      = substr($0, 39, 8)+0
    z      = substr($0, 47, 8)+0
   
    # read "seed" pdb as a guide
    if(seed_resnum != "") next
    if(atom == "N")
    {
	N["X"]=x;  N["Y"]=y;  N["Z"]=z;
	seed["N"] = resnum
    }
    if(atom == "CA")
    {
	CA["X"]=x;  CA["Y"]=y;  CA["Z"]=z;
	seed["CA"] = resnum
    }
    if(atom == "C")
    {
	C["X"]=x;  C["Y"]=y;  C["Z"]=z;
	seed["C"] = resnum
    }
    if(atom == "O")
    {
	O["X"]=x;  O["Y"]=y;  O["Z"]=z;
	seed["O"] = resnum
    }
    if((seed["N"] == seed["CA"])&&(seed["CA"] == seed["C"])&&(seed["N"]!=""))
    {
	seed_resnum = resnum
	seed_restyp = restyp
	seed_chain  = chain
	Bfac = substr($0, 61, 6)+0
    }
}

toupper($1) ~ /^BUILD/ && NF>=3{
    ++builds
    restyp = toupper($2); phi=$3; psi=$4;
    if(restyp !~ /^[A-Z][A-Z][A-Z]$/)
    {
	restyp = "ALA"
	phi=$2; psi=$3;
    }
    # cache build commands (since we have to print N->C)
    BUILD[builds] = restyp " " phi " " psi
}

toupper($1) ~ /^BUILD/ && NF==2{
    # short build command for specifying terminal phi value
    next_phi = $2
}

END{
    # see if we were "seeded"
    if(seed_resnum != "")
    {	
	# model residue numbers after this one
	resnum = seed_resnum
	chain  = seed_chain
    }
    else
    {
	resnum = builds+1
    }
    seed_resnum = resnum
    
    # count down the builds from the last one given
    for(build=builds;build>0;--build)
    {
	--resnum
	split(BUILD[build], w)
	restyp = w[1]; phi = w[2]; psi = w[3];

	# position of "this" carbonyl C follows "next" residue's phi angle
	next_atom(C,CA,N,next_phi,121.9,1.33);
	C["X"]=new_atom["X"]; C["Y"]=new_atom["Y"]; C["Z"]=new_atom["Z"];
	
	# position of "this" carbonyl O is determined by peptide bond
	next_atom(CA,N,C,0,123.2,1.23);
	O["X"]=new_atom["X"]; O["Y"]=new_atom["Y"]; O["Z"]=new_atom["Z"];
	
	# position of CA is also determined by peptide bond
	next_atom(CA,N,C,180,115.6,1.52);
	CA["X"]=new_atom["X"]; CA["Y"]=new_atom["Y"]; CA["Z"]=new_atom["Z"];
	
	# position of "this" N follows "this" residue's psi angle
	next_atom(N,C,CA,psi,110.54,1.45);
	N["X"]=new_atom["X"]; N["Y"]=new_atom["Y"]; N["Z"]=new_atom["Z"];
	
	# build up a "block" of atoms (for printing later)
	PDB_block[resnum] = sprint_atom(N,"N");
	PDB_block[resnum] = PDB_block[resnum] sprint_atom(CA,"CA");
	PDB_block[resnum] = PDB_block[resnum] sprint_atom(C,"C");
	PDB_block[resnum] = PDB_block[resnum] sprint_atom(O,"O");
	
	if(restyp != "GLY")
	{
	    # position of "this" CB should be 120 degrees out?
	    next_atom(C,N,CA,-120,110.5,1.52);
	    CB["X"]=new_atom["X"]; CB["Y"]=new_atom["Y"]; CB["Z"]=new_atom["Z"];
	    ++new_atoms;
	    PDB_block[resnum] = PDB_block[resnum] sprint_atom(CB,"CB");
	}
	
	# save this for next residue build
	next_phi = phi;
    }
    
    # now, finally, we can print out the atoms! 
    for(;resnum<seed_resnum;++resnum)
    {
	printf "%s", PDB_block[resnum]
    }
    
    # followed by the input pdb atoms
    for(oldatom=1;oldatom<=oldatoms;++oldatom)
    {
	print OldAtom[oldatom];
    }
    
    if(! atoms_printed) exit
    #print "END"
}



################################################################################
#
#	sprint_atom(atom, name, restyp, resnum, occ, Bfac)
#
#  	Function for creating a standard PDB line
#
################################################################################
function sprint_atom(atom, _name, _restyp, _resnum, _occ, _Bfac) {
# defaults (global)
if(! name)     name = "CA"
if(! restyp)   restyp = "ALA"
if(resnum=="") resnum = 1
if(! chain)    chain = ""
if(occ=="")    occ = 1
if(Bfac=="")   Bfac = 20

if(_name)   name   = toupper(_name)
if(_restyp) restyp = toupper(_restyp)
if(_resnum) resnum = toupper(resnum)
if(_occ)    occ = _occ
if(_Bfac)   Bfac = _Bfac

if(resnum ~ /^[A-Z]/) chain = substr(resnum,1,1)
while(resnum ~ /^[A-Z]/) resnum = substr(resnum,2)
++atoms_printed

entry=sprintf("ATOM %6d  %-3s %3s %1s%4d     %7.3f %7.3f %7.3f %5.2f%6.2f\n",\
atoms_printed,name,restyp,chain,resnum,atom["X"],atom["Y"],atom["Z"],occ,Bfac);

return entry
}


################################################################################
#
#	next_atom(atom1,atom2,atom3, chi, angle, bond)
#
#  	Function for getting "new atom" xyz coordinates using:
#
#	three reference atoms (defining two "previous" bonds)
#	two   angles (the bond angle, and the chi torsion angle)
#	one   distance (length of the "new" bond)
#
#
#     O -atom1                     O - new_atom
#      \                          /
#       \                        /
#        \           angle      / - bond (1.54A)
#         \               \__  /
#          \   chi(= 0)   /   /
#    atom2- O -------------- O -atom3
#
# atom1, atom2, atom3, are 3-membered arrays ["X","Y","Z"]
# new atom returned in new_atom["X","Y","Z"]
#
# note: atom1 serves only to define chi=0, it can really be anywhere
################################################################################
#
#                                   |-> optional
function next_atom(atom1,atom2,atom3, chi, angle, bond) {


    # default bond angle and length of new bond
    if(chi   == "") chi   = 0
    if(angle == "") angle = 109.5
    if(bond == "")
    {
	bond = 1.54
	# assume some double-bond character in non-tetrahedral bonds
	if(angle != 109.5) bond = 1.4
    }
    

    # compute components of "new_atom"-"atom3" vector, relative to chi==0 defined by "bond1"
     axis_component["length"] = bond*sin(3.1415927*(angle-90)/180)
     chi0_component["length"] = -bond*cos(3.1415927*(angle-90)/180)*cos(3.1415927*(chi)/180)
    chi90_component["length"] = -bond*cos(3.1415927*(angle-90)/180)*sin(3.1415927*(chi)/180)

    # now we know the coordinates of the new atom vector in "local" coordinates
    # we need to construct a basis for converting them to "global" coordinates

    # vector subtration of atoms lying on rotation axis
    axis["X"] = atom3["X"]-atom2["X"]
    axis["Y"] = atom3["Y"]-atom2["Y"]
    axis["Z"] = atom3["Z"]-atom2["Z"]

    # vector subtraction of atoms defining "zero" rotation around the axis
    bond1["X"] = atom2["X"]-atom1["X"]
    bond1["Y"] = atom2["Y"]-atom1["Y"]
    bond1["Z"] = atom2["Z"]-atom1["Z"]

    # protect against singular vectors
    if(((axis["X"]^2 + axis["Y"]^2 + axis["Z"]^2) == 0)||((bond1["X"]^2 + bond1["Y"]^2 + bond1["Z"]^2)==0))
    {
	new_atom["X"] = 0;
	new_atom["Y"] = 0;
	new_atom["Z"] = 0; 
	return 0   
    }

    # normalize the "axis" vector
    axis["length"]  = sqrt( (axis["X"]*axis["X"])  + ( axis["Y"]*axis["Y"])  + (axis["Z"]*axis["Z"]));
    axis["X"] = axis["X"]/axis["length"]
    axis["Y"] = axis["Y"]/axis["length"]
    axis["Z"] = axis["Z"]/axis["length"]
    axis["length"] = 1

    # compute amount of "bond1" that needs to be removed in order to orthogonalize it to "axis"
    bond1_dot_axis  =     ( (axis["X"]*bond1["X"]) +  (axis["Y"]*bond1["Y"]) + (axis["Z"]*bond1["Z"]) ) # Dot product

    # subtract the projection from bond1
    bond1["X"] = bond1["X"] - bond1_dot_axis*axis["X"]
    bond1["Y"] = bond1["Y"] - bond1_dot_axis*axis["Y"]
    bond1["Z"] = bond1["Z"] - bond1_dot_axis*axis["Z"]

    # normalize the "bond1" vector to make the "chi0" reference vector
    bond1["length"] = sqrt((bond1["X"]*bond1["X"]) + (bond1["Y"]*bond1["Y"]) +(bond1["Z"]*bond1["Z"]) );
    chi0["X"] = bond1["X"]/bond1["length"]
    chi0["Y"] = bond1["Y"]/bond1["length"]
    chi0["Z"] = bond1["Z"]/bond1["length"]
    chi0["length"] = 1

    # now make the "other" basis vector to complete the right-handed, 
    # orthonormal basis, using a cross-product
    chi90["X"] = axis["Y"] * chi0["Z"] - axis["Z"] * chi0["Y"];
    chi90["Y"] = axis["Z"] * chi0["X"] - axis["X"] * chi0["Z"];
    chi90["Z"] = axis["X"] * chi0["Y"] - axis["Y"] * chi0["X"];
    
    # we now have three unit vectors forming a basis of rotation about the "atom2-atom3" bond

    # the "axis" component of the new atom, offset from atom3, will be 0.514A out along "axis"
    axis_component["X"]  = axis_component["length"]*axis["X"]
    axis_component["Y"]  = axis_component["length"]*axis["Y"]
    axis_component["Z"]  = axis_component["length"]*axis["Z"]

    # apply the "x-y" values from the given chi angle
    chi0_component["X"]  = chi0_component["length"]*chi0["X"]
    chi0_component["Y"]  = chi0_component["length"]*chi0["Y"]
    chi0_component["Z"]  = chi0_component["length"]*chi0["Z"]

    chi90_component["X"] = chi90_component["length"]*chi90["X"]
    chi90_component["Y"] = chi90_component["length"]*chi90["Y"]
    chi90_component["Z"] = chi90_component["length"]*chi90["Z"]

    
    # now generate a "new" atom, in the original coordinate system
    new_atom["X"] = atom3["X"] + axis_component["X"] + chi0_component["X"] + chi90_component["X"]
    new_atom["Y"] = atom3["Y"] + axis_component["Y"] + chi0_component["Y"] + chi90_component["Y"]
    new_atom["Z"] = atom3["Z"] + axis_component["Z"] + chi0_component["Z"] + chi90_component["Z"]

    return 1
}

