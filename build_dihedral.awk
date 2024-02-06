#! /bin/awk -f
#
#
#	(re)build atoms using three "anchor" atoms and a dihedral angle
#
# to rebuild a given side chain, put a line in the PDB (before the ATOM data)
# like this:
# BUILD A23:N A23:CA A23:C A23:OXT 180
#
# I.E.
# echo "BUILD A51:N A51:CA A51:C A51:OXT 180 " |\
# cat - old.pdb |\
# build_dihedral.awk >! new.pdb
#
#  The residue ID is {chain}{number}{conformer}.  I.E. C51A represents residue
#  51 in chain C, conformer A. 
#
#  Normally, this script passes-thru all atoms from the original file.
#  If you want only the rebuilt residue to be output, put "ONLYNEW" as the
#  first word on a line at the top of the PDB
#
#
BEGIN{
    RTD = 45/atan2(1,1);
    Rm["N"]=1;Rm["A"]=2;Rm["C"]=3;Rm["O"]=4;Rm["B"]=5;
    Rm["G"]=6;Rm["D"]=7;Rm["E"]=8;Rm["Z"]=9;Rm["H"]=10;Rm["X"]=11

    alphabet[1]="A"; alphabet[2]="B"; alphabet[3]="C"; alphabet[4]="D"; alphabet[5]="E"; alphabet[6]="F"; alphabet[7]="G"; alphabet[8]="H"; alphabet[9]="I"; alphabet[10]="J"; alphabet[11]="K"; alphabet[12]="L"; alphabet[13]="M"; alphabet[14]="N"; alphabet[15]="O"; alphabet[16]="P"; alphabet[17]="Q"; alphabet[18]="R"; alphabet[19]="S"; alphabet[20]="T"; alphabet[21]="U"; alphabet[22]="V"; alphabet[23]="W"; alphabet[24]="X"; alphabet[25]="Y"; alphabet[26]="Z"
    alphabet["A"]=1; alphabet["B"]=2; alphabet["C"]=3; alphabet["D"]=4; alphabet["E"]=5; alphabet["F"]=6; alphabet["G"]=7; alphabet["H"]=8; alphabet["I"]=9; alphabet["J"]=10; alphabet["K"]=11; alphabet["L"]=12; alphabet["M"]=13; alphabet["N"]=14; alphabet["O"]=15; alphabet["P"]=16; alphabet["Q"]=17; alphabet["R"]=18; alphabet["S"]=19; alphabet["T"]=20; alphabet["U"]=21; alphabet["V"]=22; alphabet["W"]=23; alphabet["X"]=24; alphabet["Y"]=25; alphabet["Z"]=26
    alphabet[0]=" "; alphabet[" "]=0;
}

/^NEWONLY/ || /^ONLYNEW/ {onlynew = 1}

/^OCCUP/ {OVERALL_OCC = $2; next}
/^BFAC/ {OVERALL_BFAC = $2; next}
/^CONF/ {OVERALL_CONF = $2; next}

( /^CRYST/ || /^REMARK/ || /^SCALE/ ) && ! onlynew


toupper($1) ~ /^BUILD/{
    # BUILD A99:CA A99:C A100:N A100A:CA 180
    # store building commands for later
    ++builds
    atoms=0;
    for(i=2;i<=NF;++i)
    {
	# first, check for flags
	if($i=="OCC") {
	    ++i;
	    Build_occ[builds]=$i;
	    continue;
        }
	if($i=="BFAC") {
	    ++i;
	    Build_bfac[builds]=$i;
	    continue;
        }
	if($i=="CONF") {
	    ++i;
	    Build_conf[builds]=$i;
	    continue;
        }
	if($i=="BOND") {
	    ++i;
	    Build_bondlength[builds]=$i;
	    continue;
        }
	if($i=="ANGLE") {
	    ++i;
	    Build_bondangle[builds]=$i;
	    continue;
        }
	if($i ~ /^[A-Z][A-Z][A-Z]$/)
	{
	    # must be a residue type
	    Build_restyp[builds] = $i;
	    continue
	}
	if($i !~ /[^0-9+.-]/)
	{
	    # pure number arguments are angles
	    Build_chi[builds]=$i;
	    continue;
	}
	if($i ~ /:/)
	{
	    # if it contains a colon, it's an atom spec
	    ++atoms;
	
	    # set up associative memory
	    Build_atom[builds,atoms] = $i;
	    if(atoms == 3) Build_num[$i] = builds;
	}
    }
}


/^ATOM|^HETAT/{
    ++n;
    Line[n]   = $0;
    Atom[n]   = substr($0, 13, 4);
    Atom_space[n] = " ";
    if(Atom[n] !~ /^ /) Atom_space[n] = "";
    gsub(" ","",Atom[n]);
    Rem[n]=Rm[substr($0,15,1)];
    if(Rem[n]+0==0)Rem[n]=Rm[substr($0,14,1)];
    if(Rem[n]+0==0)Rem[n]=99;
    Main[n]=(Atom[n] ~ /^[NCO]$/ || Atom[n] ~ /^C[AB]$/);

    Conf[n]   = substr($0, 17, 1);
    Restyp[n] = substr($0, 18, 3);
    Chain[n]  = substr($0, 22, 1);
    Resnum[n] = substr($0, 23, 4)+0;
    X[n]   = substr($0, 31, 8)+0;
    Y[n]   = substr($0, 39, 8)+0;
    Z[n]   = substr($0, 47, 8)+0;

    Occ[n]  = substr($0, 55,6)+0;
    Bfac[n] = substr($0, 61,6)+0;

    cfmr[n] = Chain[n] Resnum[n] Conf[n];
    gsub(" ","",cfmr[n]);

    # may need sequence associations
    Restyp_of[cfmr[n]] = Restyp[n];
    Restyp_of[Chain[n] Resnum[n]] = Restyp[n];

    id = cfmr[n] ":" Atom[n] ;
    gsub(" ","",id);
    id_of[n] = id
    n_of[id] = n;

    # linked list of what to print in what order
    nxt_atom[n-1] = n;
    if(n>1) prv_atom[n] = n-1;
}

END{
    # now that we have all atom positions we can start inserting new atoms
    firstatom = 1;
    firstnew = n+1;
    for(a=1;a<=n;++a)
    {
        # see if we're ready to build something
	build = Build_num[id_of[a]];
	if(build == "" || built[build])
	{
	    # this is not a build point.  Move on.
	    continue;
	}

	# do we have all three anchor atoms?
	i = n_of[Build_atom[build,1]];
	j = n_of[Build_atom[build,2]];
	k = n_of[Build_atom[build,3]];
	if(k != a) print "WHAT THE FRAK? "

	if(i == "" || j == "" || k == "")
	{
	    # one of the atoms we need is missing!
	    printf("ERROR: cannot build %s because",Build_atom[build,4]);
	    if(i == "") printf("%s is missing. "),Build_atom[build,1];
	    if(j == "") printf("%s is missing. "),Build_atom[build,2];
	    if(k == "") printf("%s is missing. "),Build_atom[build,3];
	    print "";
	    continue;
	}

	# construct the new atom position with sensible defaults for missing values
	n = next_atom(i,j,k,Build_chi[build],Build_bondangle[build],Build_bondlength[build]);
		
	# get atom name, chain, resnum and conformer from the build command
	split(Build_atom[build,4],w,":")
	Atom[n] = w[2];
	# two-letter element name if atom is long, or contains "spaces"
        Atom_space[n]=" ";
	if( length(Atom[n])>3 || Atom[n] ~ /_/ ) Atom_space[n]="";

	Resnum[n] = w[1];
	Conf[n] = Chain[n] = ""
	if(Resnum[n] ~ /^[A-Z]/)
	{
	    Chain[n] = substr(Resnum[n],1,1)
	    Resnum[n] = substr(Resnum[n],2)
	}
	if(Resnum[n] ~ /[^0-9]$/)
	{
	    Conf[n] = substr(Resnum[n],length(Resnum[n]))
	    Resnum[n] = substr(Resnum[n],1,length(Resnum[n])-1)
	}
	while(Resnum[n] ~ /^[A-Z]/)  Resnum[n] = substr(Resnum[n],2)
	while(Resnum[n] ~ /[^0-9]$/) Resnum[n] = substr(Resnum[n],1,length(Resnum[n])-1)

	# apply any user-specified conformers
	if(OVERALL_CONF != "") Conf[n] = OVERALL_CONF;
        if(Build_conf[build] != ""){
	    Conf[n] = Build_conf[build];
	}


	# by default, inherit x-ray characteristic of fulcrum atom
	Occ[n]  = Occ[k];
	Bfac[n] = Bfac[k];
	# unless a global property has been specified
	if(OVERALL_OCC  != "") Occ[n]  = OVERALL_OCC;
	if(OVERALL_BFAC != "") Bfac[n] = OVERALL_BFAC;
	# or a specific property was specified in the build command
        if(Build_occ[build] != ""){
	    Occ[n] = Build_occ[build];
	}
        if(Build_bfac[build] != ""){
	    Bfac[n] = Build_bfac[build];
	}

	# fill in with other defaults, such as inheritance from fulcrum atom
	if(Resnum[n]=="")
	{
	    Resnum[n] = Resnum[k];
	    # these must have been missing too?
	    if(Chain[n]=="") Chain[n] = Chain[k];
	    if(Conf[n]=="") Conf[n]  = Conf[k];
	    # maybe we have moved on to the next residue?
	    if(Atom[k] == "C" && Atom[n] == "N")
	    {
	        Resnum[n]=Resnum[k]+1;
	    }
	    if(Atom[k] == "N" && Atom[n] == "C")
	    {
		# or the previous
	        Resnum[n]=Resnum[k]-1;
	    }
	}
	if(Chain[n] !~ /^[A-Z]$/) Chain[n] = " "
	if(Conf[n]  == "") Conf[n] = " "

	# maintain for future builds
	cfmr[n] = Chain[n] Resnum[n] Conf[n];
	gsub(" ","",cfmr[n]);
	# should be identical format to construction from ATOM records
	id = cfmr[n] ":" Atom[n];
	id_of[n] = id;
	n_of[id] = n;

	# all that remains now is the three-letter code
	Restyp[n] = Build_restyp[build];
	if(Restyp[n]=="") Restyp[n] = Restyp_of[cfmr[n]];
	if(Restyp[n]=="") Restyp[n] = Restyp_of[Chain[n] Resnum[n]];
	if(Restyp[n]=="") Restyp[n] = "UNK"
	Restyp_of[cfmr[n]] = Restyp[n];
	Restyp_of[Chain[n] Resnum[n]] = Restyp[n];

	# now we finally have enough information for the whole line
	Line[n]=sprintf("ATOM %6d %-4s%1s%3s %1s%4d     %7.3f %7.3f %7.3f %5.2f%6.2f",\
	  n,Atom_space[n] Atom[n],Conf[n],Restyp[n],Chain[n],Resnum[n],X[n],Y[n],Z[n],Occ[n],Bfac[n]);

	# OK. Now where do we put it?
        Rem[n]=Rm[substr(Line[n],15,1)];
        if(Rem[n]+0==0)Rem[n]=Rm[substr(Line[n],14,1)];
        if(Rem[n]+0==0)Rem[n]=99;
        Main[n]=(Atom[n] ~ /^[NCO]$/ || Atom[n] ~ /^C[AB]$/);

	ins = a;
	while(alphabet[Chain[n]]<alphabet[Chain[ins]] && prv_atom[ins] != "")
	{
	    # advance to end of insertion-point residue
	    ins=prv_atom[ins];
	}
	while(alphabet[Chain[n]]>alphabet[Chain[ins]] && nxt_atom[ins] != "")
	{
	    # advance to end of insertion-point residue
	    ins=nxt_atom[ins];
	}
	while(Resnum[n]<Resnum[ins] && Chain[n]==Chain[ins] && prv_atom[ins] != "")
	{
	    # advance to end of insertion-point residue
	    ins=prv_atom[ins];
	}
	while(Rem[n]>Rem[ins] && Resnum[n]==Resnum[ins] && Chain[n]==Chain[ins] && nxt_atom[ins] != "")
	{
	    # advance to end of insertion-point residue
	    ins=nxt_atom[ins];
	}
	while(Rem[n]<Rem[ins] && Resnum[n]==Resnum[ins] && Chain[n]==Chain[ins] && prv_atom[ins] != "")
	{
	    # advance to end of insertion-point residue
	    ins=prv_atom[ins];
	}

	prv = prv_atom[ins];
	nxt = nxt_atom[ins];
	if(Resnum[n]<Resnum[ins] || Rem[n]<Rem[ins])
	{
	    # we are building "up", "n" must be printed before "ins"
	    nxt_atom[n]   = ins;
	    nxt_atom[prv] = n;
	    
	    prv_atom[n]   = prv;
	    prv_atom[ins] = n;
	    if(prv=="") firstatom = n;
	}
	else
	{
	    # we are building "down", "n" must be printed after "ins"
	    nxt_atom[ins] = n;
	    nxt_atom[n]   = nxt;
	    
	    prv_atom[n]   = ins;
	    prv_atom[nxt] = n;
	}

        # we are done! don't do this again.
	built[build] = 1;	
    }

    # now, finally, print it all out
    maxprint = n;
    if(onlynew) firstatom = firstnew;
    for(i=firstatom;i<=n && maxprint>0;i=nxt_atom[i])
    {
	print Line[i],i;
	--maxprint;
    }
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
# atom1, atom2, atom3, are indices of the arrays X[],Y[],Z[]
# return value is index of new atom in these same arrays
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

    # find an empty slot in the "X,Y,Z" arrays.  The max index is usually "n"
    atom4 = n;
    while(X[atom4] != ""){
	++atom4;
	++n;
    }

    # compute components of "new_atom"-"atom3" vector, relative to chi==0 defined by "bond1"
     axis_component["length"] = bond*sin((angle-90)/RTD)
     chi0_component["length"] = -bond*cos((angle-90)/RTD)*cos(chi/RTD)
    chi90_component["length"] = -bond*cos((angle-90)/RTD)*sin(chi/RTD)

    # now we know the coordinates of the new atom vector in "local" coordinates
    # we need to construct a basis for converting them to "global" coordinates

    # vector subtration of atoms lying on rotation axis
    axis["X"]  = X[atom3]-X[atom2]
    axis["Y"]  = Y[atom3]-Y[atom2]
    axis["Z"]  = Z[atom3]-Z[atom2]

    # vector subtraction of atoms defining "zero" rotation around the axis
    bond1["X"] = X[atom1]-X[atom2]
    bond1["Y"] = Y[atom1]-Y[atom2]
    bond1["Z"] = Z[atom1]-Z[atom2]

    # protect against singular vectors
    if(((axis["X"]^2 + axis["Y"]^2 + axis["Z"]^2) == 0)||((bond1["X"]^2 + bond1["Y"]^2 + bond1["Z"]^2)==0))
    {
	X[atom4] = 0;
	Y[atom4] = 0;
	Z[atom4] = 0; 
	print "ERROR: trying to build from exactly overlapping atoms"
	return 0;
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
    X[atom4] = X[atom3] + axis_component["X"] + chi0_component["X"] + chi90_component["X"]
    Y[atom4] = Y[atom3] + axis_component["Y"] + chi0_component["Y"] + chi90_component["Y"]
    Z[atom4] = Z[atom3] + axis_component["Z"] + chi0_component["Z"] + chi90_component["Z"]

    return atom4;
}


