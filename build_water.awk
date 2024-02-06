#! /bin/awk -f
#
#
#	build amorphous water using the 
#
#
BEGIN{

    # initialize "phantom" atoms to start phi-psi from
    beforelast["X"]=1
    beforelast["Y"]=0
    beforelast["Z"]=-3.8

    last["X"]=0
    last["Y"]=0
    last["Z"]=-3.8

    this["X"]=0
    this["Y"]=0
    this["Z"]=0

    hbond = 2.976
    angle = 110
}

{
  maxrad = $1
  maxbox = $2
  maxatoms = $3
  if(maxrad < 3) maxrad = 99999
  if(maxbox < 3) maxbox = 99999
  if(maxatoms < 3) maxatoms = 999999

  if(maxbox < 999) printf "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f\n",maxbox*2,maxbox*2,maxbox*2,90,90,90;

  while ( resnum+0 < maxatoms )
  {
    ++resnum

    clash = trials = 1
    while ( clash )
    {
	++trials;
#	next_atom(beforelast,last,this,360*rand(),104.474+0*(rand()-0.5),2.976+0.0*(rand()-0.5));
#	next_atom(beforelast,last,this,360*rand(),angle+30*(2*(rand()-0.5)),hbond+0.0*(2*(rand()-0.5)));
	next_atom(beforelast,last,this,360*rand(),angle,hbond);

	# check for clashes
	clash = buddies = 0
	for(i=1;i<resnum-2;++i) {
	    if(new_atom["X"]<-maxbox) { clash=1; break }
	    if(new_atom["Y"]<-maxbox) { clash=1; break }
	    if(new_atom["Z"]<-maxbox) { clash=1; break }
	    if(new_atom["X"]>maxbox) { clash=1; break }
	    if(new_atom["Y"]>maxbox) { clash=1; break }
	    if(new_atom["Z"]>maxbox) { clash=1; break }
	    dist=sqrt((Xhistory[i]-new_atom["X"])^2+(Yhistory[i]-new_atom["Y"])^2+(Zhistory[i]-new_atom["Z"])^2);
	    radius=sqrt((new_atom["X"])^2+(new_atom["Y"])^2+(new_atom["Z"])^2);
	    if(radius > maxrad) {
#		print "clash with outer rim"
		clash=1;
		break;
	    }
	    if(dist < hbond-0.1) {
#		print "clash with atom" i
		clash=1;
		break;
	    }
	    if(dist < 1.5*hbond) {
		++buddies;
		buddy[buddies]=resnum;
	    }
	}
	if(buddies==0 && resnum > 4 && trials < 100) {
	    # don't let chain wander into empty space
	    clash=1
	}

	if(clash && trials > 100) {
	    # chain has run into itself, find another seed pair
	    ++seed;
	    trials = 0
	    print "REMARK changing seed to atoms", seed, pivotseed
	    this["X"]=Xhistory[seed];
	    this["Y"]=Yhistory[seed];
	    this["Z"]=Zhistory[seed];
	}
	if(! clash) seed = 1
	if(clash && ( seed > resnum ) ) {
	    print "REMARK out of seeds"
	    seed=1
	    ++pivotseed
	    last["X"]=Xhistory[pivotseed];
	    last["Y"]=Yhistory[pivotseed];
	    last["Z"]=Zhistory[pivotseed];
	}
	if(clash && ( seed > resnum ) && ( pivotseed > resnum ) ) {
	    print "REMARK out of options"
	    exit
	}
    }

    beforelast["X"]=last["X"]; beforelast["Y"]=last["Y"]; beforelast["Z"]=last["Z"];
    last["X"]=this["X"]; last["Y"]=this["Y"]; last["Z"]=this["Z"];
    this["X"]=new_atom["X"]; this["Y"]=new_atom["Y"]; this["Z"]=new_atom["Z"];

    Xhistory[resnum]=this["X"];Yhistory[resnum]=this["Y"];Zhistory[resnum]=this["Z"];
    printf "%s", sprint_atom(this,"OW","WAT",resnum);
  }
  exit
}



END{
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
    if(axis["length"] == 0)
    {
	new_atom["X"] = 0;
	new_atom["Y"] = 0;
	new_atom["Z"] = 0; 
	return 0   
    }
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
    if(bond1["length"] == 0)
    {
	new_atom["X"] = 0;
	new_atom["Y"] = 0;
	new_atom["Z"] = 0; 
	return 0   
    }
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

