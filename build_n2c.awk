#! /bin/awk -f
#
#
#	build a polypeptide chain with given phi-psi-omega angles
#	off the C-terminus of a provided pdb file
#
# command format:
# BUILD RES phi psi omega
#
BEGIN{

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

    last_omega=omega=180
}

/^CRYST/
/^REMARK/
/^SCALE/

toupper($1) ~ /^BUILD/{
    ++resnum
    restyp = toupper($2);
    if(restyp !~ /^[A-Z][A-Z][A-Z]$/)
    {
	restyp = "ALA"
	phi=$2; psi=$3; omega=$4;
	if(phi=="") phi=180
	if(psi=="") psi=180
	if(omega=="") omega=180
    }
    else
    {
	phi=$3; psi=$4; omega=$5;
	if(phi=="") phi=180
	if(psi=="") psi=180
	if(omega=="") omega=180
    }
    # arbitrary number of angles
    for(i=3;i<=NF;++i) {
	if($i=="OCC") {
	    ++i;
	    Occ=$i;
	    continue;
        }
	if($i=="BFAC") {
	    ++i;
	    Bfac=$i;
	    continue;
        }
	if($i=="CONF") {
	    ++i;
	    newconf=$i;
	    continue;
        }
	Build_angles[builds] = Build_angles[builds] " " $i
    }
        
    # print out initial N
    if((lastCA==lastO)&&(lastCA==lastC))
    {
	next_atom(O,CA,C,omega,115.6,1.33);
	N["X"]=new_atom["X"]; N["Y"]=new_atom["Y"]; N["Z"]=new_atom["Z"];
    }
    printf "%s", sprint_atom(N,"N",restyp,resnum);
    lastCA = resnum;
    
    # position of CA follows from last carbonyl carbon
    next_atom(CA,C,N,last_omega,121.9,1.45);
    CA["X"]=new_atom["X"]; CA["Y"]=new_atom["Y"]; CA["Z"]=new_atom["Z"];
    printf "%s", sprint_atom(CA,"CA",restyp,resnum);

    # position of the next carbonyl C follows from phi
    next_atom(C,N,CA,phi,110.54,1.52);
    C["X"]=new_atom["X"]; C["Y"]=new_atom["Y"]; C["Z"]=new_atom["Z"];
    printf "%s", sprint_atom(C,"C",restyp,resnum);

    # carbonyl O follows from psi+180
    next_atom(N,CA,C,psi+180,121.1,1.23);
    O["X"]=new_atom["X"]; O["Y"]=new_atom["Y"]; O["Z"]=new_atom["Z"];
    printf "%s", sprint_atom(O,"O",restyp,resnum);
    
    # now print out CB (and rest of side chain?)
    if(restyp != "GLY")
    {
	# position of CB should be 120 degrees away?
	next_atom(C,N,CA,-120,110.5,1.52);
	CB["X"]=new_atom["X"]; CB["Y"]=new_atom["Y"]; CB["Z"]=new_atom["Z"];
	printf "%s", sprint_atom(CB,"CB",restyp,resnum);
    }    

    # position of next N follows from psi
    next_atom(N,CA,C,psi,115.6,1.33);
    N["X"]=new_atom["X"]; N["Y"]=new_atom["Y"]; N["Z"]=new_atom["Z"];
    # don't print it out yet

    # use this residues omega to build the CA for the next residue
    last_omega=omega;
}


/^ATOM/{
    resnum = substr($0, 23, 4)+0
    restyp = substr($0, 18, 3)
    chain  = substr($0, 22, 1)          # O/Brookhaven-style segment ID
    split(substr($0, 13, 4), a)
    Atom   = a[1];    
    X      = substr($0, 31, 8)+0
    Y      = substr($0, 39, 8)+0
    Z      = substr($0, 47, 8)+0

    # read "seed" pdb as a guide (if available)
    if(Atom == "N")
    {
	N["X"]=X;  N["Y"]=Y;  N["Z"]=Z;
	lastN = resnum;
    }
    if(Atom == "CA") 
    {
	CA["X"]=X; CA["Y"]=Y; CA["Z"]=Z;
	lastCA = resnum;
    }
    if(Atom == "C")   
    {
	C["X"]=X;  C["Y"]=Y;  C["Z"]=Z;
	lastC = resnum;
    }
    if(Atom == "O")   
    {
	O["X"]=X;  O["Y"]=Y;  O["Z"]=Z;
	lastO = resnum;
    }
    
    # detect residue breaks
    if(resnum != last_resnum)
    {
	
    }
    last_resnum = resnum
    
    print;
}

END{

    if(! atoms_printed) exit
    
    if(! noOXT) printf "%s", sprint_atom(N,"OXT",restyp,resnum);
    print "END"
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
if(conf=="")   conf = " "

if(_name)   name   = toupper(_name)
if(_restyp) restyp = toupper(_restyp)
if(_resnum) resnum = toupper(resnum)
if(_occ)    occ = _occ
if(_Bfac)   Bfac = _Bfac

if(resnum ~ /^[A-Z]/) chain = substr(resnum,1,1)
while(resnum ~ /^[A-Z]/) resnum = substr(resnum,2)
++atoms_printed

entry=sprintf("ATOM %6d  %-3s%1s%3s %1s%4d     %7.3f %7.3f %7.3f %5.2f%6.2f\n",\
atoms_printed,name,conf,restyp,chain,resnum,atom["X"],atom["Y"],atom["Z"],occ,Bfac);

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

