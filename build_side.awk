#! /bin/awk -f
#
#
#        (re)build amino acid side chains with given chi angles                      -James Holton  3-21-23
#
# to rebuild a given side chain, put a line in the PDB (before the ATOM data)
# like this:
# BUILD RES Chain|Number chi1 chi2 chi3 ...
#
# I.E.
# echo "BUILD MET A51 180 180 180" |\
# cat - old.pdb |\
# nawk -f build_side.awk >! new.pdb
#
#  The residue ID is {chain}{number}{conformer}.  I.E. C51A represents residue
#  51A in chain C.  (o differentiates residue 51 from 51A)
#
#  chi angles can also be expressed as one of + - t or 0, which correspond to
#  the Ponder and Richards rotamer types.  Missing chi angles will not be built.
#
#  If you don't care what a particular chi angle is, and want to build the
#  most probable Ponder and Richards value indicate the indeterminate chi 
#  chi value with "?" 
#
#  for example:
# BUILD LYS A5 60 60  
#  would build Lysine a5 as a norvaline (but call it LYS)
#
# BUILD LYS A5 ? ? ? 
#  would build Lysine a5 as the - t 180 rotamer
#
#  Normally, this script passes-thru all atoms not involved in rebuilding.
#  If you want only the rebuilt residue to be output, put "ONLYNEW" as the
#  first word on a line at the top of the PDB
#
#
BEGIN{
    # user can set default Occupancy, B-factor (otherwise taken from main chain)
    if(! Occ)  Occ  = ""
    if(! Bfac) Bfac = ""
}

/^NEWONLY/ || /^ONLYNEW/ {onlynew = 1}

/^OCCUP/ {Occ = $2; next}
/^BFAC/ {Bfac = $2; next}



toupper($1) ~ /^BUILD/ && NF>=3{
    # BUILD RES A99 chi1 chi2 ...
    # store building commands for later
    ++builds
    restyp = toupper($2);
    resnum = toupper($3);
    conf = chain = ""
    while(resnum ~ /^[A-Z]/)
    {
        chain = substr(resnum,1,1)
        resnum = substr(resnum,2)
    }
    while(resnum ~ /[A-Z]$/)
    {
        conf = substr(resnum,length(resnum))
        resnum = substr(resnum,1,length(resnum)-1)
    }
    if(chain !~ /^[A-Z]$/) chain = " "
    if(conf  !~ /^[A-Z]$/) conf = " "

    # set up associative memory
    Build_typ[builds]    = restyp
    Build_num[builds]    = resnum
    Build_chain[builds]  = chain
    Build_conf[builds]   = conf
    
    # reverse association
    BUILD[chain,resnum,conf]  = builds
    BUILD[chain,resnum]  = builds
    # arbitrary number of angles
    for(i=4;i<=NF;++i) {
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
            Build_newconf[builds]=$i;
            continue;
        }
        if($i=="p") $i="+"
        if($i=="m") $i="-"
        Build_angles[builds] = Build_angles[builds] " " $i
    }
    next;
}

( ! /^ATOM|^HETAT/ ) && ! onlynew {print;next}


/^ATOM|^HETAT/{
    resnum = substr($0, 23, 4)+0
    conf   = substr($0, 17, 1)
    restyp = substr($0, 18, 3)
    chain  = substr($0, 22, 1)          # O/Brookhaven-style segment ID
    split(substr($0, 13, 4), a)
    atom   = a[1];
    x      = substr($0, 31, 8)+0
    y      = substr($0, 39, 8)+0
    z      = substr($0, 47, 8)+0
    occ    = substr($0, 55, 6)+0
    bfac   = substr($0, 61, 6)+0
    
    OCC[atom]=occ;
    BFAC[atom]=bfac;
    
    if(! BUILD[chain,resnum])
    {
        # nothing interesting to do here
        if(! onlynew) print;
        next
    }
    
    # this residue is supposed to be rebuilt, but maybe not this conf
    build  = BUILD[chain,resnum]
    restyp = Build_typ[build]
    if(atom == "N")
    {
        N["X"]=x;  N["Y"]=y;  N["Z"]=z; gotN=build;
    }
    if(atom == "CA") 
    {
        CA["X"]=x; CA["Y"]=y; CA["Z"]=z; gotCA=build;
    }
    if(atom == "C")   
    {
        C["X"]=x;  C["Y"]=y;  C["Z"]=z; gotC=build;
    }
    if((atom == "O")||(atom == "OXT"))   
    {
        O["X"]=x;  O["Y"]=y;  O["Z"]=z; gotO=build;
    }

    if((atom == "CB")&&(gotN==build)&&(gotCA==build)&&(gotC!=build))   
    {
        # was carbonyl C missing? 
        CB["X"]=x;  CB["Y"]=y;  CB["Z"]=Z;
        next_atom(CB,N,CA,120,110.54,1.52)
        C["X"]=new_atom["X"]; C["Y"]=new_atom["Y"]; C["Z"]=new_atom["Z"];
        gotC=build;
    }
    

    if(! BUILD[chain,resnum,conf])
    {
        # this is the wrong conformer, but maybe backbone will be useful
        if(! onlynew) print;
        next
    }
    else
    {
        if(atom=="OXT") print;
    }
    
    # see if we're ready to build the side chain
    if((gotN==build)&&(gotCA==build)&&(gotC==build)&&(gotO==build))
    {
      for(build=1;build<=builds;++build){
        if(Build_num[build]   != resnum) continue;
        if(Build_chain[build] != chain) continue;

        if(Build_occ[build] != ""){
            Occ = Build_occ[build];
            OCC["N"]=OCC["CA"]=OCC["C"]=OCC["O"]=Occ;
        }
        if(Build_bfac[build] != ""){
            Bfac = Build_bfac[build];
            BFAC["N"]=BFAC["CA"]=BFAC["C"]=BFAC["O"]=Bfac;
        }
        if(Build_newconf[build] != ""){
            conf = Build_newconf[build];
        }

        # print out main chain (in right order)
        main_chain =            sprint_atom(N,"N", "", "", OCC["N"], BFAC["N"]) 
        main_chain = main_chain sprint_atom(CA,"CA", "", "", OCC["CA"], BFAC["CA"]) 
        main_chain = main_chain sprint_atom(C,"C", "", "", OCC["C"], BFAC["C"]) 
        main_chain = main_chain sprint_atom(O,"O", "", "", OCC["O"], BFAC["O"]) 
        printf "%s", main_chain
        
        # we have everything we need to build the side chain
        split(Build_angles[build], chi)

        if(restyp == "GLY") side_chain = ""
        if(restyp == "ALA") side_chain = build_ALA(N, CA, C);

        if(restyp == "SER") side_chain = build_SER(N,CA,C, chi[1]);
        if(restyp == "CYS") side_chain = build_CYS(N,CA,C, chi[1]);
        if(restyp == "VAL") side_chain = build_VAL(N,CA,C, chi[1]);
        if(restyp == "THR") side_chain = build_THR(N,CA,C, chi[1]);
        if(restyp == "ABU") side_chain = build_ABU(N,CA,C, chi[1]);

        if(restyp == "LEU") side_chain = build_LEU(N,CA,C, chi[1], chi[2]);
        if(restyp == "ILE") side_chain = build_ILE(N,CA,C, chi[1], chi[2]);
        if(restyp == "IIL") side_chain = build_IIL(N,CA,C, chi[1], chi[2]);
        if(restyp == "PRO") side_chain = build_PRO(N,CA,C, chi[1], chi[2]);
        if(restyp == "ASP") side_chain = build_ASP(N,CA,C, chi[1], chi[2]);
        if(restyp == "ASN") side_chain = build_ASN(N,CA,C, chi[1], chi[2]);
        if(restyp == "NRV") side_chain = build_NRV(N,CA,C, chi[1], chi[2]);

        if(restyp == "HIS") side_chain = build_HIS(N,CA,C, chi[1], chi[2]);
        if(restyp == "PHE") side_chain = build_PHE(N,CA,C, chi[1], chi[2]);
        if(restyp == "TYR") side_chain = build_TYR(N,CA,C, chi[1], chi[2]);
        if(restyp == "TRP") side_chain = build_TRP(N,CA,C, chi[1], chi[2]);

        if(restyp == "MET") side_chain = build_MET(N,CA,C, chi[1], chi[2], chi[3]);
        if(restyp == "NRL") side_chain = build_NRL(N,CA,C, chi[1], chi[2], chi[3]);
        if(restyp == "GLU") side_chain = build_GLU(N,CA,C, chi[1], chi[2], chi[3]);
        if(restyp == "GLN") side_chain = build_GLN(N,CA,C, chi[1], chi[2], chi[3]);

        if(restyp == "LYS") side_chain = build_LYS(N,CA,C, chi[1], chi[2], chi[3], chi[4]);
        if(restyp == "GLA") side_chain = build_GLA(N,CA,C, chi[1], chi[2], chi[3], chi[4]);

        if(restyp == "ARG") side_chain = build_ARG(N,CA,C, chi[1], chi[2], chi[3], chi[4], chi[5]);
        
        # print it out
        printf "%s", side_chain
      }
        # reset finder flags
        gotN=gotCA=gotC=gotO=0;
        
    }
}



END{

    if(! atoms_printed) exit
    
    if(build=="N2C") print sprint_atom(N,"OXT",Restyp,Resnum,Conf);
    #print "END"
}



################################################################################
#
#        build_ALA(N,CA,C)
#
#          Function for building an alanine side chain
#
################################################################################
function build_ALA(N,CA,C) {
    # position of CB should be 120 degrees away?
    next_atom(C,N,CA,-120,110.5,1.52);
    CB["X"]=new_atom["X"]; CB["Y"]=new_atom["Y"]; CB["Z"]=new_atom["Z"];

    side_chain = sprint_atom(CB,"CB");
    return side_chain
}


################################################################################
#
#        build_SER(N,CA,C,chi1)
#
#          Function for building a serine side chain
#
################################################################################
function build_SER(N,CA,C,chi1) {
    if(chi1=="?") chi1 = "+"
    
    if(chi1=="+") chi1 = 64.7
    if(chi1=="-") chi1 = -69.7
    if(chi1=="t") chi1 = -176.1

    # use alanine builder
    side_chain = build_ALA(N,CA,C);

    # honor request to not build
    if(chi1=="") return side_chain

    # add oxygen with chi angle
    next_atom(N,CA,CB,chi1);
    OG["X"]=new_atom["X"]; OG["Y"]=new_atom["Y"]; OG["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(OG,"OG");
    return side_chain
}


################################################################################
#
#        build_CYS(N,CA,C,chi1)
#
#          Function for building a cystine side chain
#
################################################################################
function build_CYS(N,CA,C,chi1) {
    if(chi1=="?") chi1 = "-"
    
    if(chi1=="+") chi1 = 63.5
    if(chi1=="-") chi1 = -65.2
    if(chi1=="t") chi1 = -179.6

    # use alanine builder
    side_chain = build_ALA(N,CA,C);

    # honor request to not build
    if(chi1=="") return side_chain

    # add sulphur with chi angle
    next_atom(N,CA,CB,chi1,113,1.82);
    SG["X"]=new_atom["X"]; SG["Y"]=new_atom["Y"]; SG["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(SG,"SG");
    return side_chain
}
################################################################################
#
#        build_ABU(N,CA,C,chi1)
#
#          Function for building an aminobutric acid side chain
#
################################################################################
function build_ABU(N,CA,C,chi1) {
    if(chi1=="?") chi1 = "t"
    
    if(chi1=="+") chi1 = 60
    if(chi1=="-") chi1 = -60
    if(chi1=="t") chi1 = 180

    # use alanine builder for CB
    side_chain = build_ALA(N,CA,C);

    # honor request to not build
    if(chi1=="") return side_chain

    # add CG with chi angle
    next_atom(N,CA,CB,chi1);
    CG["X"]=new_atom["X"]; CG["Y"]=new_atom["Y"]; CG["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CG,"CG");
    return side_chain
}

################################################################################
#
#        build_NRV(N,CA,C,chi1,chi2)
#
#          Function for building a norvaline side chain
#
################################################################################
function build_NRV(N,CA,C,chi1,chi2) {
    if(chi1=="?") chi1 = "t"
    if(chi2=="?") chi2 = "t"

    if(chi1=="+") chi1 = 60
    if(chi1=="-") chi1 = -60
    if(chi1=="t") chi1 = 180
    
    if(chi2=="+") chi2 = 60
    if(chi2=="-") chi2 = -60
    if(chi2=="t") chi2 = 180
    
    # use abu builder for CB & CG
    side_chain = build_ABU(N,CA,C,chi1);

    # honor request to not build
    if(chi2=="") return side_chain

    # add CD with 2nd chi angle
    next_atom(CA,CB,CG,chi2);
    CD["X"]=new_atom["X"]; CD["Y"]=new_atom["Y"]; CD["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CD,"CD");
    return side_chain
}


################################################################################
#
#        build_PRO(N,CA,C,chi1,chi2)
#
#          Function for building a proline side chain
#
################################################################################
function build_PRO(N,CA,C,chi1,chi2) {
    if(chi1=="?") chi1 = "+"
    if(chi1=="+") {chi1 =  26.9; chi2= -29.4}
    if(chi1=="-") {chi1 = -21.8; chi2=  31.2}
    if(chi1=="0") {chi1 =   0.3; chi2=  -0.8}
    if(chi2=="?") chi2 = -29.4
    
    # position of CB should be 120 degrees away?
    next_atom(C,N,CA,-120,103.5,1.5);
    CB["X"]=new_atom["X"]; CB["Y"]=new_atom["Y"]; CB["Z"]=new_atom["Z"];
    side_chain = sprint_atom(CB,"CB");

    # honor request to not build
    if(chi1=="") return side_chain

    # add CG with 1st chi angle
    next_atom(N,CA,CB,chi1,105,1.5);
    CG["X"]=new_atom["X"]; CG["Y"]=new_atom["Y"]; CG["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CG,"CG");

    # honor request to not build
    if(chi2=="") return side_chain

    # add CD with 2nd chi angle
    next_atom(CA,CB,CG,chi2,105,1.5);
    CD["X"]=new_atom["X"]; CD["Y"]=new_atom["Y"]; CD["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CD,"CD");

    return side_chain
}


################################################################################
#
#        build_NRL(N,CA,C,chi1,chi2,chi3)
#
#          Function for building a norleucine side chain
#
################################################################################
function build_NRL(N,CA,C,chi1,chi2,chi3) {
    if(chi1=="?") chi1 = "t"
    if(chi2=="?") chi2 = "t"
    if(chi3=="?") chi3 = "t"

    if(chi1=="+") chi1 = 60
    if(chi1=="-") chi1 = -60
    if(chi1=="t") chi1 = 180
    
    if(chi2=="+") chi2 = 60
    if(chi2=="-") chi2 = -60
    if(chi2=="t") chi2 = 180
    
    if(chi3=="+") chi3 = 60
    if(chi3=="-") chi3 = -60
    if(chi3=="t") chi3 = 180
    
    # use norvaline builder for CB, CG &CD
    side_chain = build_NRV(N,CA,C,chi1,chi2);

    # honor request to not build
    if(chi3=="") return side_chain

    # add CE with 3rd chi angle
    next_atom(CB,CG,CD,chi3);
    CE["X"]=new_atom["X"]; CE["Y"]=new_atom["Y"]; CE["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CE,"CE");
    return side_chain
}


################################################################################
#
#        build_MET(N,CA,C,chi1,chi2,chi3)
#
#          Function for building a methionine side chain
#
################################################################################
function build_MET(N,CA,C,chi1,chi2,chi3) {
    if(chi1=="?") chi1 = "-"
    if(chi2=="?") chi2 = "-"
    if(chi3=="?") chi3 = "-"
    if((chi1=="-")&&(chi2=="t")) {chi1= -78.3; chi2= -174.7}
    if((chi1=="t")&&(chi2=="t")) {chi1= 178.9; chi2= 179.0}
    
    if(chi1=="+") chi1 = 60
    if(chi1=="-") chi1 = -64.5
    if(chi1=="t") chi1 = 178.9
    
    if(chi2=="+") chi2 = 60
    if(chi2=="-") chi2 = -68.5
    if(chi2=="t") chi2 = 180
    
    if(chi3=="+") chi3 = 60
    if(chi3=="-") chi3 = -75.6
    if(chi3=="t") chi3 = 180
    
    # use abu builder for CB & CG
    side_chain = build_ABU(N,CA,C,chi1);

    # honor request to not build
    if(chi2=="") return side_chain

    # add SD with 2nd chi angle
    next_atom(CA,CB,CG,chi2,113,1.78);
    SD["X"]=new_atom["X"]; SD["Y"]=new_atom["Y"]; SD["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(SD,"SD");

    # honor request to not build
    if(chi3=="") return side_chain

    # add CE with 3rd chi angle
    next_atom(CB,CG,SD,chi3,109.5,1.78);
    CE["X"]=new_atom["X"]; CE["Y"]=new_atom["Y"]; CE["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CE,"CE");
    
    return side_chain
}




################################################################################
#
#        build_THR(N,CA,C,chi1)
#
#          Function for building a threonine side chain
#
################################################################################
function build_THR(N,CA,C,chi1) {
    if(chi1=="?") chi1 = "+"
    
    if(chi1=="+") chi1 = 62.7
    if(chi1=="-") chi1 = -59.7
    if(chi1=="t") chi1 = -169.5
    
    # use alanine builder
    side_chain = build_ALA(N,CA,C);

    # honor request to not build rest of side chain
    if(chi1=="") return side_chain

    # add CG and OG with chi angle
    next_atom(N,CA,CB,chi1);
    OG1["X"]=new_atom["X"]; OG1["Y"]=new_atom["Y"]; OG1["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(OG1,"OG1");
    
    next_atom(N,CA,CB,chi1-120);
    CG2["X"]=new_atom["X"]; CG2["Y"]=new_atom["Y"]; CG2["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CG2,"CG2");

    return side_chain
}


################################################################################
#
#        build_VAL(N,CA,C,chi1)
#
#          Function for building a valine side chain
#
################################################################################
function build_VAL(N,CA,C,chi1) {
    if(chi1=="?") chi1 = "t"
    
    if(chi1=="+") chi1 = 69.3
    if(chi1=="-") chi1 = -63.4
    if(chi1=="t") chi1 = 173.5
    
    # use alanine builder
    side_chain = build_ALA(N,CA,C);

    # honor request to not build rest of side chain
    if(chi1=="") return side_chain

    # add CGs with chi angle
    next_atom(N,CA,CB,chi1);
    CG1["X"]=new_atom["X"]; CG1["Y"]=new_atom["Y"]; CG1["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CG1,"CG1");
    
    next_atom(N,CA,CB,chi1+120);
    CG2["X"]=new_atom["X"]; CG2["Y"]=new_atom["Y"]; CG2["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CG2,"CG2");

    return side_chain
}


################################################################################
#
#        build_ILE(N,CA,C,chi1,chi2)
#
#          Function for building an isoleucie side chain
#
################################################################################
function build_ILE(N,CA,C,chi1,chi2) {
    if(chi1=="?") chi1 = "-"
    if(chi2=="?") chi2 = "t"
    if((chi1=="-")&&(chi2=="t")) {chi1= -60.9; chi2= 168.7}
    if((chi1=="-")&&(chi2=="-")) {chi1= -59.6; chi2= -64.1}
    if((chi1=="+")&&(chi2=="t")) {chi1=  61.7; chi2= 163.8}
    if((chi1=="t")&&(chi2=="t")) {chi1=-166.6; chi2= 166.0}
    if((chi1=="t")&&(chi2=="+")) {chi1=-174.8; chi2=  72.1}

    if(chi1=="+") chi1 = 61.7
    if(chi1=="-") chi1 = -60.9
    if(chi1=="t") chi1 = -166.6
    
    if(chi2=="+") chi2 = 72.1
    if(chi2=="-") chi2 = -64.1
    if(chi2=="t") chi2 = 168.7

    # use alanine builder
    side_chain = build_ALA(N,CA,C);

    # honor request to not build rest of side chain
    if(chi1=="") return side_chain

    # add CGs with 1st chi angle
    next_atom(N,CA,CB,chi1);
    CG1["X"]=new_atom["X"]; CG1["Y"]=new_atom["Y"]; CG1["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CG1,"CG1");
    
    next_atom(N,CA,CB,chi1-120);
    CG2["X"]=new_atom["X"]; CG2["Y"]=new_atom["Y"]; CG2["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CG2,"CG2");

    # honor request to not build rest of side chain
    if(chi2=="") return side_chain

    # add CD1 with 2nd chi angle
    next_atom(CA,CB,CG1,chi2);
    CD1["X"]=new_atom["X"]; CD1["Y"]=new_atom["Y"]; CD1["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CD1,"CD1");

    return side_chain
}


################################################################################
#
#        build_IIL(N,CA,C,chi1,chi2)
#
#          Function for building an isoisoleucine side chain
#
################################################################################
function build_IIL(N,CA,C,chi1,chi2) {
    if(chi1=="?") chi1 = "-"
    if(chi2=="?") chi2 = "t"
    if((chi1=="-")&&(chi2=="t")) {chi1= -60.9; chi2= 168.7}
    if((chi1=="-")&&(chi2=="-")) {chi1= -59.6; chi2= -64.1}
    if((chi1=="+")&&(chi2=="t")) {chi1=  61.7; chi2= 163.8}
    if((chi1=="t")&&(chi2=="t")) {chi1=-166.6; chi2= 166.0}
    if((chi1=="t")&&(chi2=="+")) {chi1=-174.8; chi2=  72.1}

    if(chi1=="+") chi1 = 61.7
    if(chi1=="-") chi1 = -60.9
    if(chi1=="t") chi1 = -166.6
    
    if(chi2=="+") chi2 = 72.1
    if(chi2=="-") chi2 = -64.1
    if(chi2=="t") chi2 = 168.7

    # use valine builder
    side_chain = build_VAL(N,CA,C,chi1);

    # honor request to not build rest of side chain
    if(chi2=="") return side_chain

    # add CD1 with 2nd chi angle
    next_atom(CA,CB,CG1,chi2);
    CD1["X"]=new_atom["X"]; CD1["Y"]=new_atom["Y"]; CD1["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CD1,"CD1");

    return side_chain
}


################################################################################
#
#        build_LEU(N,CA,C,chi1,chi2)
#
#          Function for building a valine side chain
#
################################################################################
function build_LEU(N,CA,C,chi1,chi2) {
    if(chi1=="?") chi1 = "-"
    if(chi2=="?") chi2 = "t"
    if((chi1=="-")&&(chi2=="t")) {chi1= -64.9; chi2= 176.0}
    if((chi1=="t")&&(chi2=="+")) {chi1=-176.4; chi2=  63.1}
    if((chi1=="t")&&(chi2=="t")) {chi1=-165.3; chi2= 168.2}
    if((chi1=="+")&&(chi2=="+")) {chi1=  44.3; chi2=  60.4}

    if(chi1=="+") chi1 = 44.3
    if(chi1=="-") chi1 = -64.9
    if(chi1=="t") chi1 = -176.4
    
    if(chi2=="+") chi2 = 63.1
    if(chi2=="-") chi2 = -60
    if(chi2=="t") chi2 = 176.0
    
    # use ABU builder for CB & CG
    side_chain = build_ABU(N,CA,C,chi1);

    # honor request to not build rest of side chain
    if(chi2=="") return side_chain

    # add CDs with 2nd chi angle
    next_atom(CA,CB,CG,chi2,113);
    CD1["X"]=new_atom["X"]; CD1["Y"]=new_atom["Y"]; CD1["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CD1,"CD1");
    next_atom(CA,CB,CG,chi2+120,113);
    CD2["X"]=new_atom["X"]; CD2["Y"]=new_atom["Y"]; CD2["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CD2,"CD2");

    return side_chain
}



################################################################################
#
#        build_ASP(N,CA,C,chi1,chi2)
#
#          Function for building an aspartate side chain
#
################################################################################
function build_ASP(N,CA,C,chi1,chi2) {
    if(chi1=="?") chi1 = "-"
    if(chi2=="?") chi2 = -25.7

    if(chi1=="-") {chi1= -68.3; chi2= -25.7}
    if(chi1=="t") {chi1=-169.1; chi2=   3.9}
    if(chi1=="+") {chi1=  63.7; chi2=   2.4}
    if(chi1=="0") {chi1=   0.0; chi2=   0.0}
    
    # use ABU builder for CB & CG
    side_chain = build_ABU(N,CA,C,chi1);

    # honor request to not build rest of side chain
    if(chi2=="") return side_chain

    # add ODs with 2nd chi angle
    next_atom(CA,CB,CG,chi2, 120,1.29);
    OD1["X"]=new_atom["X"]; OD1["Y"]=new_atom["Y"]; OD1["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(OD1,"OD1");
    next_atom(CA,CB,CG,chi2+180, 120,1.29);
    OD2["X"]=new_atom["X"]; OD2["Y"]=new_atom["Y"]; OD2["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(OD2,"OD2");

    return side_chain
}


################################################################################
#
#        build_ASN(N,CA,C,chi1,chi2)
#
#          Function for building an asparagine side chain
#
################################################################################
function build_ASN(N,CA,C,chi1,chi2) {
    if(chi1=="?") chi1 = "-"
    if(chi2=="?") chi2 = "-"
    if((chi1=="-")&&(chi2=="-")) {chi1= -68.3; chi2= -36.0}
    if((chi1=="t")&&(chi2=="0")) {chi1=-177.1; chi2=   1.3}
    if((chi1=="-")&&(chi2=="+")) {chi1= -67.2; chi2= 128.8}
    if((chi1=="+")&&(chi2=="0")) {chi1=  63.9; chi2=  -6.8}
    if((chi1=="t")&&(chi2=="t")) {chi1=-174.9; chi2=-156.8}
    if((chi1=="+")&&(chi2=="+")) {chi1=  63.6; chi2=  53.8}

    if(chi1=="+") chi1 =  63.9
    if(chi1=="-") chi1 = -68.3
    if(chi1=="t") chi1 = -177.1
    
    if(chi2=="+") chi2 = 128.8
    if(chi2=="-") chi2 = -36.0
    if(chi2=="t") chi2 =-156.8
    
    # use ABU builder for CB & CG
    side_chain = build_ABU(N,CA,C,chi1);

    # honor request to not build rest of side chain
    if(chi2=="") return side_chain

    # add O/NDs with 2nd chi angle
    next_atom(CA,CB,CG,chi2, 120,1.23);
    OD1["X"]=new_atom["X"]; OD1["Y"]=new_atom["Y"]; OD1["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(OD1,"OD1");
    next_atom(CA,CB,CG,chi2+180, 120,1.32);
    ND2["X"]=new_atom["X"]; ND2["Y"]=new_atom["Y"]; ND2["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(ND2,"ND2");

    return side_chain
}

################################################################################
#
#        build_GLU(N,CA,C,chi1,chi2,chi3)
#
#          Function for building a glutamate side chain
#
################################################################################
function build_GLU(N,CA,C,chi1,chi2,chi3) {
    if(chi1=="?") chi1 = "-"
    if(chi2=="?") chi2 = "t"
    if(chi3=="?") chi3 = "0"
    if((chi1=="-")&&(chi2=="t")) {chi1= -69.6; chi2=-177.2; chi3= -11.4}
    if((chi1=="t")&&(chi2=="t")) {chi1=-176.2; chi2= 175.4; chi3=  -6.7}
    if((chi1=="-")&&(chi2=="-")) {chi1= -64.6; chi2= -69.1; chi3= -33.4}
    if((chi1=="-")&&(chi2=="+")) {chi1= -55.6; chi2=  77.0; chi3=  25.3}
    if((chi1=="+")&&(chi2=="t")) {chi1=  69.8; chi2=-179.0; chi3=   6.6}
    if((chi1=="t")&&(chi2=="+")) {chi1=-173.6; chi2=  70.6; chi3=  14.0}
    if((chi1=="+")&&(chi2=="-")) {chi1=  63.0; chi2= -80.4; chi3=  16.3}

    if(chi1=="+") chi1 =  69.8
    if(chi1=="-") chi1 = -69.6
    if(chi1=="t") chi1 = -176.2
    
    if(chi2=="+") chi2 = 77.0
    if(chi2=="-") chi2 = -69.1
    if(chi2=="t") chi2 =-177.2

    if(chi3=="+") chi3 =  25.3
    if(chi3=="-") chi3 = -11.4
    if(chi3=="t") chi3 = 180

    # use norvaline builder for CB, CG & CD
    side_chain = build_NRV(N,CA,C,chi1,chi2);

    # honor request to not build rest of side chain
    if(chi3=="") return side_chain

    # add OEs with 3rd chi angle
    next_atom(CB,CG,CD,chi3,120,1.29);
    OE1["X"]=new_atom["X"]; OE1["Y"]=new_atom["Y"]; OE1["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(OE1,"OE1");
    next_atom(CB,CG,CD,chi3+180,120,1.29);
    OE2["X"]=new_atom["X"]; OE2["Y"]=new_atom["Y"]; OE2["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(OE2,"OE2");

    return side_chain
}

################################################################################
#
#        build_GLA(N,CA,C,chi1,chi2,chi3,chi4)
#
#          Function for building a gamma-carboxy-glutamate side chain
#
################################################################################
function build_GLA(N,CA,C,chi1,chi2,chi3,chi4) {
    if(chi1=="?") chi1 = "-"
    if(chi2=="?") chi2 = "t"
    if(chi3=="?") chi3 = "0"
    if((chi1=="-")&&(chi2=="t")) {chi1= -69.6; chi2=-177.2; chi3= -11.4}
    if((chi1=="t")&&(chi2=="t")) {chi1=-176.2; chi2= 175.4; chi3=  -6.7}
    if((chi1=="-")&&(chi2=="-")) {chi1= -64.6; chi2= -69.1; chi3= -33.4}
    if((chi1=="-")&&(chi2=="+")) {chi1= -55.6; chi2=  77.0; chi3=  25.3}
    if((chi1=="+")&&(chi2=="t")) {chi1=  69.8; chi2=-179.0; chi3=   6.6}
    if((chi1=="t")&&(chi2=="+")) {chi1=-173.6; chi2=  70.6; chi3=  14.0}
    if((chi1=="+")&&(chi2=="-")) {chi1=  63.0; chi2= -80.4; chi3=  16.3}

    if(chi1=="+") chi1 =  69.8
    if(chi1=="-") chi1 = -69.6
    if(chi1=="t") chi1 = -176.2
    
    if(chi2=="+") chi2 = 77.0
    if(chi2=="-") chi2 = -69.1
    if(chi2=="t") chi2 =-177.2

    if(chi3=="+") chi3 =  25.3
    if(chi3=="-") chi3 = -11.4
    if(chi3=="t") chi3 = 180

    if(chi4=="") chi4=chi3;

    # use glutamate builder for first branch
    side_chain = build_GLU(N,CA,C,chi1,chi2,chi3);
    
    # edit the side chain so that CD becomes CD1
    sub(" CD ", " CD1", side_chain)
    
    # add 2nd CD with 2nd chi angle
    next_atom(CA,CB,CG,chi2+120);
    CD2["X"]=new_atom["X"]; CD2["Y"]=new_atom["Y"]; CD2["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CD2,"CD2");
    
    # honor request to not build rest of side chain
    if(chi4=="") return side_chain

    # add OEs with 4th chi angle
    next_atom(CB,CG,CD2,chi4,120,1.29);
    OE3["X"]=new_atom["X"]; OE3["Y"]=new_atom["Y"]; OE3["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(OE3,"OE3");
    next_atom(CB,CG,CD2,chi4+180,120,1.29);
    OE4["X"]=new_atom["X"]; OE4["Y"]=new_atom["Y"]; OE4["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(OE4,"OE4");

    return side_chain
}




################################################################################
#
#        build_GLN(N,CA,C,chi1,chi2,chi3)
#
#          Function for building a glutamine side chain
#
################################################################################
function build_GLN(N,CA,C,chi1,chi2,chi3) {
    if(chi1=="?") chi1 = "-"
    if(chi2=="?") chi2 = "t"
    if(chi3=="?") chi3 = "0"
    if((chi1=="-")&&(chi2=="t"))              {chi1= -66.7; chi2=-178.5}
    if((chi1=="t")&&(chi2=="t"))              {chi1=-174.6; chi2= 177.7}
    if((chi1=="-")&&(chi2=="-")&&(chi3=="0")) {chi1= -58.7; chi2= -63.8; chi3= -46.3}
    if((chi1=="t")&&(chi2=="+")&&(chi3=="0")) {chi1=-179.4; chi2=  67.3; chi3=  26.8}
    if((chi1=="+")&&(chi2=="t"))              {chi1=  70.8; chi2=-165.6}
    if((chi1=="-")&&(chi2=="-")&&(chi3=="t")) {chi1= -51.3; chi2= -90.4; chi3= 165.0}
    if((chi1=="t")&&(chi2=="+")&&(chi3=="t")) {chi1= 167.5; chi2=  70.9; chi3= 174.2}

    if(chi1=="+") chi1 =  70.8
    if(chi1=="-") chi1 = -66.7
    if(chi1=="t") chi1 = -174.6
    
    if(chi2=="+") chi2 = 67.3
    if(chi2=="-") chi2 = -63.8
    if(chi2=="t") chi2 = -178.5

    if(chi3=="+") chi3 =  26.8
    if(chi3=="-") chi3 = -46.3
    if(chi3=="t") chi3 = 165.0

    # use norvaline builder for CB, CG & CD
    side_chain = build_NRV(N,CA,C,chi1,chi2);

    # honor request to not build rest of side chain
    if(chi3=="") return side_chain

    # add O/NEs with 3rd chi angle
    next_atom(CB,CG,CD,chi3,120,1.23);
    OE1["X"]=new_atom["X"]; OE1["Y"]=new_atom["Y"]; OE1["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(OE1,"OE1");
    next_atom(CB,CG,CD,chi3+180,120,1.32);
    NE2["X"]=new_atom["X"]; NE2["Y"]=new_atom["Y"]; NE2["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(NE2,"NE2");

    return side_chain
}




################################################################################
#
#        build_LYS(N,CA,C,chi1,chi2,chi3,chi4)
#
#          Function for building a lysine side chain
#
################################################################################
function build_LYS(N,CA,C,chi1,chi2,chi3,chi4) {
    if(chi1=="?") chi1 = "-"
    if(chi2=="?") chi2 = "t"
    if(chi3=="?") chi3 = "t"
    if(chi4=="?") chi4 = "t"
    if((chi1=="-")&&(chi2=="t")) {chi1= -68.9; chi2=-178.4}
    if((chi1=="t")&&(chi2=="t")) {chi1=-172.1; chi2= 175.3}
    if((chi1=="-")&&(chi2=="-")) {chi1= -58.1; chi2= -74.9}
    if((chi1=="t")&&(chi2=="+")) {chi1= 173.4; chi2=  83.4}
    if((chi1=="+")&&(chi2=="t")) {chi1=  71.5; chi2=-174.3}
    if((chi1=="t")&&(chi2=="-")) {chi1=-175.8; chi2= -63.9}
    if((chi1=="-")&&(chi2=="+")) {chi1=-104.0; chi2=  74.6}

    if(chi1=="+") chi1 =  63.9
    if(chi1=="-") chi1 = -68.3
    if(chi1=="t") chi1 = -177.1
    
    if(chi2=="+") chi2 = 128.8
    if(chi2=="-") chi2 = -36.0
    if(chi2=="t") chi2 =-156.8
    
    if(chi3=="+") chi3 =  60
    if(chi3=="-") chi3 = -60
    if(chi3=="t") chi3 = 180
    
    if(chi4=="+") chi4 =  60
    if(chi4=="-") chi4 = -60
    if(chi4=="t") chi4 = 180
    
    # use norleuine builder for CB, CG & CD
    side_chain = build_NRL(N,CA,C,chi1,chi2,chi3);

    # honor request to not build rest of side chain
    if(chi4=="") return side_chain

    # add NZ with 4th chi angle
    next_atom(CG,CD,CE,chi4);
    NZ["X"]=new_atom["X"]; NZ["Y"]=new_atom["Y"]; NZ["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(NZ,"NZ");

    return side_chain
}





################################################################################
#
#        build_ARG(N,CA,C,chi1,chi2,chi3,chi4,chi5)
#
#          Function for building an arginine side chain
#
################################################################################
function build_ARG(N,CA,C,chi1,chi2,chi3,chi4,chi5) {
    if(chi1=="?") chi1 = "-"
    if(chi2=="?") chi2 = "t"
    if(chi3=="?") chi3 = "t"
    if(chi4=="?") chi4 = "t"
    if(chi5=="") chi5 = "0"
    if((chi1=="-")&&(chi2=="t")) {chi1= -67.6; chi2= 176.9}
    if((chi1=="t")&&(chi2=="t")) {chi1=-174.1; chi2=-178.6}
    if((chi1=="+")&&(chi2=="t")) {chi1=  80.0; chi2= 175.6}
    if((chi1=="-")&&(chi2=="-")) {chi1= -67.0; chi2= -71.7}
    if((chi1=="t")&&(chi2=="+")) {chi1= 178.2; chi2=  69.5}
    if((chi1=="+")&&(chi2=="+")) {chi1=  57.1; chi2=  82.8}
    if((chi1=="-")&&(chi2=="+")) {chi1= -76.9; chi2=  54.2}
    
    if(chi1=="+") chi1 =  80.0
    if(chi1=="-") chi1 = -67.6
    if(chi1=="t") chi1 = -174.1
    
    if(chi2=="+") chi2 =  69.5
    if(chi2=="-") chi2 = -71.7
    if(chi2=="t") chi2 = 176.9
    
    if(chi3=="+") chi3 =  60
    if(chi3=="-") chi3 = -60
    if(chi3=="t") chi3 = 180
    
    if(chi4=="+") chi4 =  60
    if(chi4=="-") chi4 = -60
    if(chi4=="t") chi4 = 180
    
    # use norvaline builder for CB, CG & CD
    side_chain = build_NRV(N,CA,C,chi1,chi2);

    # honor request to not build rest of side chain
    if(chi3=="") return side_chain

    # add NE with 3rd chi angle
    next_atom(CB,CG,CD,chi3, 109.5,1.47);
    NE["X"]=new_atom["X"]; NE["Y"]=new_atom["Y"]; NE["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(NE,"NE");

    # honor request to not build rest of side chain
    if(chi4=="") return side_chain

    # add CZ with 4th chi angle
    next_atom(CG,CD,NE,chi4, 124,1.33);
    CZ["X"]=new_atom["X"]; CZ["Y"]=new_atom["Y"]; CZ["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CZ,"CZ");

    # honor request to not build rest of side chain
    if(chi5=="") return side_chain

    # add NHs with 5th chi angle
    next_atom(CD,NE,CZ,chi5, 122,1.33);
    NH1["X"]=new_atom["X"]; NH1["Y"]=new_atom["Y"]; NH1["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(NH1,"NH1");
    next_atom(CD,NE,CZ,chi5+180, 122,1.33);
    NH2["X"]=new_atom["X"]; NH2["Y"]=new_atom["Y"]; NH2["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(NH2,"NH2");

    return side_chain
}





################################################################################
#
#        build_HIS(N,CA,C,chi1,chi2)
#
#          Function for building a histidine side chain
#
################################################################################
function build_HIS(N,CA,C,chi1,chi2) {
    if(chi1=="?") chi1 = "-"
    if(chi2=="?") chi2 = "-"
    if((chi1=="-")&&(chi2=="-")) {chi1= -62.8; chi2= -74.3}
    if((chi1=="t")&&(chi2=="+")) {chi1=-175.2; chi2= -87.7}
    if((chi1=="-")&&(chi2=="+")) {chi1=  69.8; chi2=  96.1}
    if((chi1=="+")&&(chi2=="-")) {chi1=  67.9; chi2= -80.5}
    if((chi1=="t")&&(chi2=="-")) {chi1=-177.3; chi2= 100.5}
    if((chi1=="+")&&(chi2=="+")) {chi1=  48.0; chi2=  85.9}

    if(chi1=="+") chi1 = 69.8
    if(chi1=="-") chi1 = -62.8
    if(chi1=="t") chi1 = -175.2
    
    if(chi2=="+") chi2 = 96.1
    if(chi2=="-") chi2 = -74.3
    if(chi2=="t") chi2 = 180
    
    # use ABU builder for CB & CG
    side_chain = build_ABU(N,CA,C,chi1);

    # honor request to not build rest of side chain
    if(chi2=="") return side_chain

    # add ND1 with 2nd chi angle
    next_atom(CA,CB,CG,chi2, 123,1.34);
    ND1["X"]=new_atom["X"]; ND1["Y"]=new_atom["Y"]; ND1["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(ND1,"ND1");
    # add CD2 180 degrees away
    next_atom(CA,CB,CG,chi2+180, 131,1.34);
    CD2["X"]=new_atom["X"]; CD2["Y"]=new_atom["Y"]; CD2["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CD2,"CD2");

    # finish the imidazole ring
    next_atom(CB,CG,ND1,180, 106.5,1.34);
    CE1["X"]=new_atom["X"]; CE1["Y"]=new_atom["Y"]; CE1["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CE1,"CE1");
    next_atom(CB,CG,CD2,180, 109.5,1.34);
    NE2["X"]=new_atom["X"]; NE2["Y"]=new_atom["Y"]; NE2["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(NE2,"NE2");

    return side_chain
}



################################################################################
#
#        build_PHE(N,CA,C,chi1,chi2)
#
#          Function for building a phenylalanine side chain
#
################################################################################
function build_PHE(N,CA,C,chi1,chi2) {
    if(chi1=="?") chi1 = "-"
    if(chi2=="?") chi2 = "90"
    if((chi1=="-")&&(chi2=="90")) {chi1= -66.3; chi2=  94.3}
    if((chi1=="t")&&(chi2=="90")) {chi1=-179.2; chi2=  78.9}
    if((chi1=="+")&&(chi2=="90")) {chi1=  66.0; chi2=  90.7}
    if((chi1=="-")&&(chi2=="0" )) {chi1= -71.9; chi2=  -0.4}

    if(chi1=="+") chi1 = 66.0
    if(chi1=="-") chi1 = -66.3
    if(chi1=="t") chi1 = -179.2
    
    if(chi2=="0")  chi2 = 0
    if(chi2=="90") chi2 = 94.3
    
    # use ABU builder for CB & CG
    side_chain = build_ABU(N,CA,C,chi1);

    # honor request to not build rest of side chain
    if(chi2=="") return side_chain

    # add CDs with 2nd chi angle
    next_atom(CA,CB,CG,chi2,120);
    CD1["X"]=new_atom["X"]; CD1["Y"]=new_atom["Y"]; CD1["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CD1,"CD1");
    next_atom(CA,CB,CG,chi2+180,120);
    CD2["X"]=new_atom["X"]; CD2["Y"]=new_atom["Y"]; CD2["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CD2,"CD2");

    # finish off the ring
    next_atom(CB,CG,CD1,180,120);
    CE1["X"]=new_atom["X"]; CE1["Y"]=new_atom["Y"]; CE1["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CE1,"CE1");
    next_atom(CB,CG,CD2,180,120);
    CE2["X"]=new_atom["X"]; CE2["Y"]=new_atom["Y"]; CE2["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CE2,"CE2");
    next_atom(CG,CD1,CE1,0,120);
    CZ["X"]=new_atom["X"]; CZ["Y"]=new_atom["Y"]; CZ["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CZ,"CZ");

    return side_chain
}


################################################################################
#
#        build_TYR(N,CA,C,chi1,chi2)
#
#          Function for building a tyrosine side chain
#
################################################################################
function build_TYR(N,CA,C,chi1,chi2) {
    if(chi1=="?") chi1 = "-"
    if(chi2=="?") chi2 = "90"
    if((chi1=="-")&&(chi2=="90")) {chi1= -66.5; chi2=  96.6}
    if((chi1=="t")&&(chi2=="90")) {chi1=-179.7; chi2=  71.9}
    if((chi1=="+")&&(chi2=="90")) {chi1=  63.3; chi2=  99.1}
    if((chi1=="-")&&(chi2=="0" )) {chi1= -67.2; chi2=  -1.0}

    if(chi1=="+") chi1 = 63.3
    if(chi1=="-") chi1 = -66.5
    if(chi1=="t") chi1 = -179.7
    
    if(chi2=="0")  chi2 = 0
    if(chi2=="90") chi2 = 96.6
    
    # use PHE builder for CB & CG
    side_chain = build_PHE(N,CA,C,chi1,chi2);

    # honor request to not build rest of side chain
    if(chi2=="") return side_chain

    next_atom(CD1,CE1,CZ,180,120);
    OH["X"]=new_atom["X"]; OH["Y"]=new_atom["Y"]; OH["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(OH,"OH");

    return side_chain
}


################################################################################
#
#        build_TRP(N,CA,C,chi1,chi2)
#
#          Function for building a tryptophan side chain
#
################################################################################
function build_TRP(N,CA,C,chi1,chi2) {
    if(chi1=="?") chi1 = "-"
    if(chi2=="?") chi2 = "+"
    if((chi1=="-")&&(chi2=="+")) {chi1= -70.4; chi2= 100.5}
    if((chi1=="+")&&(chi2=="-")) {chi1=  64.8; chi2= -88.9}
    if((chi1=="t")&&(chi2=="-")) {chi1=-177.3; chi2= -95.1}
    if((chi1=="t")&&(chi2=="+")) {chi1=-179.5; chi2=  87.5}
    if((chi1=="-")&&(chi2=="-")) {chi1= -73.3; chi2= -87.7}
    if((chi1=="+")&&(chi2=="+")) {chi1=  62.2; chi2= 112.5}

    if(chi1=="+") chi1 = 64.8
    if(chi1=="-") chi1 = -70.4
    if(chi1=="t") chi1 = -177.3
    
    if(chi2=="+") chi2 = 100.5
    if(chi2=="-") chi2 = -88.9
    
    # use ABU builder for CB & CG
    side_chain = build_ABU(N,CA,C,chi1);

    # honor request to not build rest of side chain
    if(chi2=="") return side_chain

    # add CD1 with 2nd chi angle
    next_atom(CA,CB,CG,chi2, 127,1.4);
    CD1["X"]=new_atom["X"]; CD1["Y"]=new_atom["Y"]; CD1["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CD1,"CD1");
    # add CD2 180 degrees away
    next_atom(CA,CB,CG,chi2+180, 127,1.4);
    CD2["X"]=new_atom["X"]; CD2["Y"]=new_atom["Y"]; CD2["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CD2,"CD2");

    # finish the rings...
    next_atom(CB,CG,CD1,180, 110.5,1.4);
    NE1["X"]=new_atom["X"]; NE1["Y"]=new_atom["Y"]; NE1["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(NE1,"NE1");

    next_atom(CB,CG,CD2,180, 108.5,1.4);
    CE2["X"]=new_atom["X"]; CE2["Y"]=new_atom["Y"]; CE2["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CE2,"CE2");

    next_atom(CB,CG,CD2,0, 133.3,1.4);
    CE3["X"]=new_atom["X"]; CE3["Y"]=new_atom["Y"]; CE3["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CE3,"CE3");

    next_atom(CG,CD2,CE2,180, 120,1.4);
    CZ2["X"]=new_atom["X"]; CZ2["Y"]=new_atom["Y"]; CZ2["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CZ2,"CZ2");

    next_atom(CG,CD2,CE3,180, 120,1.4);
    CZ3["X"]=new_atom["X"]; CZ3["Y"]=new_atom["Y"]; CZ3["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CZ3,"CZ3");

    next_atom(CD2,CE2,CZ2,0, 120,1.4);
    CH2["X"]=new_atom["X"]; CH2["Y"]=new_atom["Y"]; CH2["Z"]=new_atom["Z"];
    side_chain = side_chain sprint_atom(CH2,"CH2");

    return side_chain
}



















################################################################################
#
#        sprint_atom(atom, name, restyp, resnum, occ, bfac)
#
#          Function for creating a standard PDB line
#
################################################################################
function sprint_atom(atom, _name, _restyp, _resnum, _occ, _bfac, _conf) {
if(! name)     name = "CA"
if(! conf)     conf = " "
if(! restyp)   restyp = "ALA"
if(resnum=="") resnum = 1
if(! chain)    chain = ""

# default parameter values
if(_occ=="")   _occ = Occ
if(_bfac=="")  _bfac = Bfac
if(_name)   name   = toupper(_name)
if(_restyp) restyp = toupper(_restyp)
if(_resnum) resnum = toupper(resnum)
if(_conf)   conf = _conf
if(_occ)    occ  = _occ
if(_bfac)   bfac = _bfac

if(resnum ~ /^[A-Z]/) chain = substr(resnum,1,1)
while(resnum ~ /^[A-Z]/) resnum = substr(resnum,2)
if(resnum ~ /[A-Z]$/) resletter = substr(resnum,length(resnum),1)

++atoms_printed

entry=sprintf("ATOM %6d  %-3s%1s%3s %1s%4d     %7.3f %7.3f %7.3f %5.2f%6.2f\n",\
atoms_printed,name,conf,restyp,chain,resnum,atom["X"],atom["Y"],atom["Z"],occ,bfac);


return entry
}














################################################################################
#
#        next_atom(atom1,atom2,atom3, chi, angle, bond)
#
#          Function for getting "new atom" xyz coordinates using:
#
#        three reference atoms (defining two "previous" bonds)
#        two   angles (the bond angle, and the chi torsion angle)
#        one   distance (length of the "new" bond)
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







function dihedral(atom1,atom2,atom3,atom4) {
#
#     O -atom1                     O - atom4
#      \                          /
#       \                        /
#        \                      /
#         \                    /
#          \      chi = 0     /
#    atom2- O -------------- O -atom3
#
# atom1, atom2, atom3, atom4, are 3-membered arrays ["X","Y","Z"]
# return value is chi
#
    
    # we need to construct a basis for converting them to "global" coordinates
    
    # get vector of first dihedral bond
    bond1["X"] = atom1["X"]-atom2["X"]
    bond1["Y"] = atom1["Y"]-atom2["Y"]
    bond1["Z"] = atom1["Z"]-atom2["Z"]
    
    # get vector of "second" (rotating axis) bond
    axis["X"]  = atom3["X"]-atom2["X"]
    axis["Y"]  = atom3["Y"]-atom2["Y"]
    axis["Z"]  = atom3["Z"]-atom2["Z"]

    # get vector of "third" dihedral bond
    bond3["X"] = atom4["X"]-atom3["X"]
    bond3["Y"] = atom4["Y"]-atom3["Y"]
    bond3["Z"] = atom4["Z"]-atom3["Z"]
    
    # normalize the "axis" to unit length
    norm = sqrt(axis["X"]^2 + axis["Y"]^2 + axis["Z"]^2)
    if(norm == 0) return "axis error"
    axis["X"] = axis["X"]/norm
    axis["Y"] = axis["Y"]/norm
    axis["Z"] = axis["Z"]/norm
    
    # reduce "bond" to their components perpendicular to the "axis"
    component  = bond1["X"]*axis["X"] + bond1["Y"]*axis["Y"] + bond1["Z"]*axis["Z"]
    bond1["X"] = bond1["X"]-component*axis["X"]
    bond1["Y"] = bond1["Y"]-component*axis["Y"]
    bond1["Z"] = bond1["Z"]-component*axis["Z"]
    
    component  = bond3["X"]*axis["X"] + bond3["Y"]*axis["Y"] + bond3["Z"]*axis["Z"]
    bond3["X"] = bond3["X"]-component*axis["X"]
    bond3["Y"] = bond3["Y"]-component*axis["Y"]
    bond3["Z"] = bond3["Z"]-component*axis["Z"]
    
    
    # now the angle between bond1 and bond3 is the dihedral angle
    

    # normalize the first and last bond vectors
    norm = sqrt(bond1["X"]^2 + bond1["Y"]^2 + bond1["Z"]^2)
    if(norm == 0) return "bond1 error"
    bond1["X"] = bond1["X"]/norm
    bond1["Y"] = bond1["Y"]/norm
    bond1["Z"] = bond1["Z"]/norm
    
    norm = sqrt(bond3["X"]^2 + bond3["Y"]^2 + bond3["Z"]^2)
    if(norm == 0) return "bond3 error"
    bond3["X"] = bond3["X"]/norm
    bond3["Y"] = bond3["Y"]/norm
    bond3["Z"] = bond3["Z"]/norm


    # construct a vector perpendicular to both the axis and bond1
    # (this differentiates "sides" of the dihedral)
    chi90["X"] = axis["Y"] * bond1["Z"] - axis["Z"] * bond1["Y"];
    chi90["Y"] = axis["Z"] * bond1["X"] - axis["X"] * bond1["Z"];
    chi90["Z"] = axis["X"] * bond1["Y"] - axis["Y"] * bond1["X"];
    
    
    # get the component of bond3 along bond1
    adjacent = bond1["X"]*bond3["X"] + bond1["Y"]*bond3["Y"] + bond1["Z"]*bond3["Z"]
    
    # get the component of bond3 along bond1
    opposite = chi90["X"]*bond3["X"] + chi90["Y"]*bond3["Y"] + chi90["Z"]*bond3["Z"]
    
    # use ArcTan to get the angle
    angle = atan2(opposite, adjacent)*180/3.1415927;
    
    return angle
}

