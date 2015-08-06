#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#  
#  This file is part of ESPResSo++.
#  
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. 


# -*- coding: utf-8 -*-
import math
from collections import namedtuple, defaultdict
import espressopp
from topology_helper import *
from operator import itemgetter # for sorting a dict


"""This Python module allows one to use GROMACS data files as the
   input to an ESPResSo++ simulation, set interactions for given
   particle types and convert GROMACS potential tables into
   ESPResSo++ tables.
   It containts functions: read(), setInteractions(), convertTable()
"""

GromacsTopology = namedtuple(
    'GromacsTopolog', [
        'defaults',
        'types',
        'masses',
        'charges',
        'res_ids',
        'atomtypeparams',
        'bondtypes',
        'bondtypeparams',
        'angletypes',
        'angletypeparams',
        'dihedraltypes',
        'dihedraltypeparams',
        'impropertypes',
        'impropertypeparams',
        'pairtypes',
        'pairtypeparams',
        'nonbond_params',
        'exclusions',
        'resname',
        'x', 'y', 'z',
        'vx', 'vy', 'vz',
        'Lx', 'Ly', 'Lz'
    ])


def read(gro_file, top_file="", doRegularExcl=True, defines=None, return_tuple=True):
    """ Read GROMACS data files.

    Keyword arguments:
    gro_file -- contains coordinates of all particles, the number of particles, velocities and box size.
    top_file -- contains topology information. Included topology files (.itp) are also read
    doRegularExcl -- if True, exclusions are generated automatically based on the nregxcl parameter (see gromacs manual)
    defines: Dictionary with #define <key> <value>.
    return_tuple: If set to True then output is in the old-fashion way. (default: True)
    """
    
    if defines is None:
        defines = {}

    x, y, z = [], [], []
    vx, vy, vz = [], [], []
    resname = []
    # read gro file
    if gro_file != "":
        f = open(gro_file)
        f.readline() # skip comment line
        total_num_particles = int(f.readline())
        
        # store coordinates and velocities
        for i in range(total_num_particles):
            line = f.readline()
            s = line[20:69]
            resname.append(line[5:8])
            # coordinates
            x.append(float(s[0:8]))
            y.append(float(s[8:16]))
            z.append(float(s[16:24]))
            
            if len(s.split()) > 3:
                # velocities
                vx.append(float(s[24:32]))
                vy.append(float(s[32:40]))
                vz.append(float(s[40:49]))
        
        # store box size
        Lx, Ly, Lz = map(float, f.readline().split()) # read last line, convert to float
        f.close()

    # read top and itp files
    masses, charges = [], [] # mases and charges of the whole configuration
    types=[] # tuple: atomindex(int) to atomtypeid(int)
    bonds={} # dict: key bondtypeid value: tuple of bond pairs
    angles={} # dict: key angletype value: tuple of triples
    dihedrals = {} #dict: key is tuple of dtypeid, value: tuple of quadruples
    impropers = {} #dict: key is tuple of dtypeid, value: tuple of quadruples
    exclusions=[] #list of atom pairs no considered in non-bonded interactions
    pairs_1_4 = {} #list of atom pairs with 1-4 interaction (scaled non-bonded interaction)
    atomtypes={} # a dict: key atomtypename(str) value: atomtypeid(int)
    
    defaults={} # gromacs default values
    atomtypeparams={} # a dict: key atomtypeid , value : class storing actual parameters of each type e.g. c6, c12, etc..
    use_atomtypeparams = {}  # dict with the atomtypes that are use in the topology
    nonbond_params = {}
    use_nonbond_params = {}
    bondtypeparams={} # same for bonds
    angletypeparams={} # same for angles
    dihedraltypeparams={} # same for dihedrals
    impropertypeparams={} # same for dihedrals
    pairtypeparams={}
    use_pairtypeparams={}

    if top_file != "":
        #f = open(top_file)
        # FileBuffer: a class which behaves like a file, but all lines are in memory
        # we use this for emulating a 'preprocessor' which handles the #include
        # statements in the .top and .itp files
        fb=FileBuffer()
        
        FillFileBuffer(top_file, fb, defines=defines)
        f = PostProcessFileBuffer(fb, defines)

        print "Reading top file: "+top_file
        line = ''
        a = 0
        defaults={} # gromacs default values
        bondtypes={} # a dict: key atomindex(int),atomindex(int)  value: bondtypeid(int)
        angletypes={} # a dict: key atomindex(int), atomindex(int),atomindex(int) value: angletypeid(int)
        dihedraltypes={} # a dict: key atomtindex(int), atomindex(int), atomindex(int),atomindex(int) value: tuple of dihedraltypeid(int)
        impropertypes={} # a dict: key atomtindex(int), atomindex(int), atomindex(int),atomindex(int) value: tuple of dihedraltypeid(int)
        atnum_attype = {}
        
        molecules=[]
        skip_section = False

        current_section = None
        previous_section = None
       
        lineindex=-1
        for line in f.lines:
            lineindex+=1
            line = line.strip()
           
            if line[0] == ";":  # skip comment line
                continue

            if skip_section and line.startswith('#end'):
                skip_section = False
                continue

            if skip_section:
                continue

            if line.startswith('#ifdef'):
                define_tmp = line.split()
                if len(define_tmp) > 1:
                    skip_section = defines.get(define_tmp[1], False)
                else:
                    skip_section = True
                continue

            if line.startswith('#else'):
                skip_section = True
                continue

            if line.startswith("#define"):
                define_tmp = line.split()
                defines[define_tmp[1]] = True
                continue

            if line.startswith('['):  # section part
                previous_section = current_section
                current_section = line.replace('[', '').replace(']', '').strip()
                continue
                
            if current_section == 'defaults':
                fields=line.split()
                # TODO: support for genpairs
                if len(fields)==5: 
                    defaults={"nbtype": int(fields[0]), "combinationrule":int(fields[1]),
                    "genpairs":fields[2] == 'yes', "fudgeLJ": float(fields[3]),
                    "fudgeQQ":float(fields[4])}
                else: 
                    defaults={"nbtype": int(fields[0]),
                              "combinationrule":int(fields[1])}

            if previous_section == 'atomtypes' and current_section != 'atomtypes':
                # add dihedral wildcard atomtyp
                atomtypes.update({'X':a})
                atomtype_wildcard = a
                # prints gromacs type and esp++ type
                for t in sorted(atomtypes.items(), key=itemgetter(1)):
                    print " %s: %d"%(t[0],t[1])

            if current_section == 'atomtypes':
                fields=line.split()
                attypename = fields[0]
                
                #make a map containing the properties
                # sig, eps may be c6 and c12: this is specified in the defaults
                # and converted later

                if fields[0].startswith('opls'):
                    tmpprop = {
                        'atnum': fields[1],
                        'mass': float(fields[3]),
                        'charge': float(fields[4]),
                        'particletype': fields[5],
                        'sig': float(fields[6]),
                        'eps': float(fields[7])
                        }
                    atnum_attype[fields[1]] = attypename
                elif len(fields)==7:
                    tmpprop={
                            "atnum":int(fields[1]),
                            "atname": fields[0], 
                            "mass": float(fields[2]),
                            "charge":float(fields[3]),
                            "particletype":fields[4],
                            "sig":float(fields[5]),
                            "eps":float(fields[6])}
                elif len(fields) == 8:
                    tmpprop = {
                        'atnum': fields[1],
                        'mass': float(fields[3]),
                        'charge': float(fields[4]),
                        'sig': float(fields[6]),
                        'eps': float(fields[7])
                        }
                else:
                    print('atomtypes other')
                    tmpprop={
                    "atnum": fields[0], "mass":float(fields[1]),
                    "charge":float(fields[2]), "particletype":fields[3],
                    "sig":float(fields[4]), "eps":float(fields[5])}
                    
                if attypename not in atomtypes:
                    atomtypes.update({attypename:a}) # atomtypes is used when reading the "atoms" section
                    atomtypeparams.update({a:tmpprop})
                    a += 1

            if current_section == 'nonbond_params':
                fields = l.split(';')[0].split()
                if len(fields) > 0:
                    a1, a2, fn, c6, c12 = fields[:5]
                    if int(fn) != 1:
                        continue
                    at1, at2 = sorted([atomtypes.get(a1), atomtypes.get(a2)])
                    if (at1, at2) not in nonbond_params:
                        nonbond_params[(at1, at2)] = {
                            'sig': float(c6),
                            'eps': float(c12)
                            }

            if current_section == 'pairtypes':
                fields = l.split(';')[0].split()
                if len(fields) > 0:
                    a1, a2, fn, c6, c12 = fields[:5]
                    if int(fn) != 1:
                        continue
                    at1, at2 = sorted([atomtypes.get(a1), atomtypes.get(a2)])
                    if at1 and at2:
                        if (at1, at2) not in pairtypeparams:
                            pairtypeparams[(at1, at2)] = {
                                'sig': float(c6),
                                'eps': float(c12)
                                }

            if current_section == 'bondtypes':
                tmp = line.split()
                    # i: i-atomname  i, j: j-atomname
                pairs = map(atomtypes.get, tmp[:2])
                if None in pairs:
                    pairs = [atomtypes.get(atnum_attype.get(i, i)) for i in tmp[:2]]
                    if None in pairs:
                        continue
                i, j = pairs
                p=ParseBondTypeParam(line)
                #check if this type has been defined before
                bdtypeid=FindType(p, bondtypeparams)
                if bdtypeid==None:
                    bdtypeid=len(bondtypeparams)
                    bondtypeparams.update({bdtypeid:p})

                if i in bondtypes:
                    bondtypes[i].update({j:bdtypeid})
                else:
                    bondtypes.update({i:{j:bdtypeid}})

            if current_section == 'angletypes':
                tmp = line.split()
                triplet = map(atomtypes.get, tmp[:3])
                if None in triplet:
                    triplet = [atomtypes.get(atnum_attype.get(i, i)) for i in tmp[:3]]
                    if None in triplet:
                        continue
                i, j, k = triplet
                p=ParseAngleTypeParam(line)

                atypeid=FindType(p, angletypeparams)
                if atypeid==None:
                    atypeid=len(angletypeparams)
                    angletypeparams.update({atypeid:p})

                if i in angletypes:
                    if j in angletypes[i]:
                        angletypes[i][j].update({k:atypeid})
                    else:
                        angletypes[i].update({j:{k:atypeid}})
                else:
                    angletypes.update({i:{j:{k:atypeid}}})

            if current_section == 'dihedraltypes':
                tmp = line.split()
                try:
                    fn = int(tmp[4])
                except ValueError:
                    continue
                quadruplet = map(atomtypes.get, tmp[:4])
                if None in quadruplet:
                    quadruplet = [atomtypes.get(atnum_attype.get(xx, xx)) for xx in tmp[:4]]
                    if None in quadruplet:
                        continue
                i, j, k, l = quadruplet
                if fn == 4:  # imporper dihedrals
                    parse_fn = ParseImproperTypeParam
                    d_types = impropertypes
                    d_type_params = impropertypeparams
                else:
                    parse_fn = ParseDihedralTypeParam
                    d_types = dihedraltypes
                    d_type_params = dihedraltypeparams

                p=parse_fn(line)

                dtypeid=FindType(p, d_type_params)
                if dtypeid==None:
                    dtypeid=len(d_type_params)
                    d_type_params.update({dtypeid:p})
                if i in d_types:
                    if j in d_types[i]:
                        if k in d_types[i][j]:
                            if l in d_types[i][j][k]:
                                d_types[i][j][k][l] += dtypeid
                            else:
                                d_types[i][j][k].update({l: dtypeid})
                        else:
                            d_types[i][j].update({k:{l:dtypeid}})
                    else:
                        d_types[i].update({j:{k:{l:dtypeid}}})
                else:
                    d_types.update({i:{j:{k:{l:dtypeid}}}})
            
            if current_section == 'molecules':
                print " "+line.strip('\n')
                mol, nrmol = line.split()
                #we have to check if the same molecules comes multiple times in the molecules section
                if len(molecules) == 0:
                    molecules.append({'name':mol, 'count':int(nrmol)})
                elif molecules[-1]['name'] == mol: #check if mol was added earlier already
                    molecules[-1]['count'] = molecules[-1]['count'] + int(nrmol) #update count
                else: molecules.append({'name':mol, 'count':int(nrmol)}) #if mol newly added
              
        
        molstartindex=0 #this is the index of the first atom in the molecule being parsed
        res_idx = 0 # index of molecule like single polymer chain or protein
        
        f.seek(0) # Now we search for bonds, angles definitions and start from the beginning of the file buffer

        for mol in molecules: ### this loop not modified for 1-4 pairs
            print "Preparing %d %s molecules... " %(mol['count'], mol['name']) 
            print "-----------------------------"

            # find and store number of molecules
            num_molecule_copies=mol['count']
            # this does not what the name suggests....
            nrexcl = storeMolecules(f, molecules, mol)
            # find and store atom types
            types, masses, charges, num_atoms_molecule, res_ids = \
                storeAtoms(f, defaults, types, atomtypes,
                           atomtypeparams, use_atomtypeparams, nonbond_params, use_nonbond_params,
                           masses, charges, num_molecule_copies, res_idx)

            # find and store bonds
            bonds = storeBonds(f, types, bondtypes, bondtypeparams, 
                               bonds, num_atoms_molecule, num_molecule_copies, molstartindex, 
                               exclusions, nrexcl, doRegularExcl)

            # find and store angles
            angles = storeAngles(f, types, angletypes, angletypeparams, angles,
                                 num_atoms_molecule, num_molecule_copies, molstartindex)

            # find and store dihedrals
            dihedrals = storeDihedrals(f, types, dihedraltypes, dihedraltypeparams, dihedrals,
                                       num_atoms_molecule, num_molecule_copies, molstartindex, atomtype_wildcard)

            # find and store impropers
            impropers = storeImpropers(f, types, impropertypes, impropertypeparams, impropers,
                                       num_atoms_molecule, num_molecule_copies, molstartindex, 
                                       atomtype_wildcard)

            pairs_1_4 = storePairs(f, defaults, types, pairtypeparams, use_pairtypeparams,
                                   atomtypeparams, pairs_1_4,
                                   num_atoms_molecule, num_molecule_copies,
                                   molstartindex)

            molstartindex+=num_molecule_copies*num_atoms_molecule
            res_idx += num_molecule_copies
            
            
    # Update typeparams. Store only those one which are really necessary. 
    bondtypeparams = {k: v for k, v in bondtypeparams.iteritems() if k in bonds}
    angletypeparams = {k: v for k, v in angletypeparams.iteritems() if k in angles}
    dihedraltypeparams = {k: v for k, v in dihedraltypeparams.iteritems() if k in dihedrals}

    try:
      del atomtypes['X'] #don't export wildcard atomtype
    except KeyError:
      pass

    # The data is packed into a touple, unpackvars contains a string which
    # tells the user which kind of data was read.
    print('Found default values {}'.format(defaults))
    print('Found {} particles'.format(len(x)))
    print('Found {} types'.format(len(types)))
    print('Found {} nonbonded_pairs'.format(len(use_nonbond_params)))
    print('Found {} masses'.format(len(masses)))
    print('Found {} charges'.format(len(charges)))
    print('Found {} atom type parameters'.format(len(use_atomtypeparams)))
    print('Found {} bonds'.format(sum(map(len, bonds.values()))))
    print('Found {} bond type parameters'.format(len(bondtypeparams)))
    print('Found {} angles'.format(sum(map(len, angles.values()))))
    print('Found {} angle type parameters'.format(len(angletypeparams)))
    print('Found {} dihedrals'.format(sum(map(len, dihedrals.values()))))
    print('Found {} dihedral type parameters'.format(len(dihedraltypeparams)))
    print('Found {} 1-4 pair type parameters'.format(len(use_pairtypeparams)))
    print('Found {} 1-4 pairs'.format(sum(map(len, pairs_1_4.values()))))
    print('Found {} bond exclusions'.format(len(exclusions)))
    print('Found box: {}'.format([Lx, Ly, Lz]))

    if return_tuple:
        params = []
        unpackvars=[]

        try:
            del atomtypes['X'] #don't export wildcard atomtype
        except KeyError:
            pass

        # The data is packed into a touple, unpackvars contains a string which
        # tells the user which kind of data was read.

        if len(defaults) != 0:
            unpackvars.append("defaults")
            params.append(defaults)
        if len(types) != 0:
            unpackvars.append("types")
            params.append(types)
            unpackvars.append("atomtypes")
            params.append(atomtypes)
        if len(masses) != 0:
            unpackvars.append("masses")
            params.append(masses)
        if len(charges) != 0:
            unpackvars.append("charges")
            params.append(charges)
        if len(atomtypeparams) !=0:
            unpackvars.append("atomtypeparameters")
            params.append(use_atomtypeparams)
        if len(bonds) != 0:
            unpackvars.append("bondtypes")
            params.append(bonds)
        if len(bondtypeparams) !=0:
            unpackvars.append("bondtypeparams")
            params.append(bondtypeparams)        
        if len(angles) != 0:
            unpackvars.append("angletypes")
            params.append(angles)
        if len(angletypeparams) != 0:
            unpackvars.append("angletypeparams")
            params.append(angletypeparams)      
        if len(dihedrals) != 0:
            unpackvars.append("dihedraltypes")
            params.append(dihedrals)
        if len(dihedraltypeparams) != 0:
            unpackvars.append("dihedraltypeparams")
            params.append(dihedraltypeparams)       
        if len(impropers) != 0:
            unpackvars.append("impropertypes")
            params.append(impropers)
        if len(impropertypeparams) != 0:
            unpackvars.append("impropertypeparams")
            params.append(impropertypeparams)
        if len(exclusions) != 0:
            unpackvars.append("exclusions")
            params.append(exclusions)  
        if len(pairs_1_4) != 0:
            unpackvars.append("pairs_1_4")
            params.append(pairs_1_4)
            unpackvars.append("pairstypeparams")
            params.append(use_pairtypeparams)
        if len(res_ids) != 0:
            unpackvars.append('res_ids')
            params.append(res_ids)
        unpackvars.append("x, y, z")
        params.extend([x, y, z])
        if len(vx) != 0:
            params.extend([vx, vy, vz])
            unpackvars.append("vx, vy, vz")
        
        unpackvars.append("resname")
        params.append(resname)

        params.extend([Lx, Ly, Lz])
        unpackvars.append("Lx, Ly, Lz")

        print "USAGE: unpack as"
        s=""
        for i in range(len(unpackvars)):
            s+=str(unpackvars[i])
            if (i< len(unpackvars)-1): s+=", "
        print s, "=gromacs.read( ... )"

        return tuple(params)
    else:
        return GromacsTopology(
            defaults, 
            types, 
            masses, 
            charges, 
            res_ids, 
            use_atomtypeparams,
            bonds, 
            bondtypeparams, 
            angles, 
            angletypeparams,
            dihedrals, 
            dihedraltypeparams, 
            impropers, 
            impropertypeparams,
            pairs_1_4, 
            use_pairtypeparams,
            use_nonbond_params, 
            exclusions, 
            resname,
            x, y, z, 
            vx, vy, vz, 
            Lx, Ly, Lz)


def storeMolecules(f, molecules, mol=""):
    nrexcl=0
    line = ''
    line=f.readlastline()
    while not 'moleculetype' in line:
        line = f.readline()
        if not line: break # break out of while if EOF
    line = f.readline()
    while(len(line) > 1 and not '[' in line):
        if line[0] == ";":   # skip comment lines
            #print "skipping line: "+line.strip("\n")
            line = f.readline()
            continue
        fields=line.split()
        #mol = fields[0]
        nrexcl=int(fields[1])
        line = f.readline()
    return nrexcl


def storeAtoms(f, defaults, types, atomtypes,
               atomtypeparams,
               use_atomtypeparams,
               nonbondedparams,
               use_nonbond_params,
               masses,
               charges,
               num_molecule_copies,
               res_idx):
    types_tmp = []
    charge_tmp =[]
    mass_tmp=[]
    molecule_index = []
    pos = f.tell()

    combinationrule = defaults['combinationrule']

    line=f.readlastline()
    while not 'atoms' in line:
        line = f.readline()
        if not line: break # break out of while if EOF
    line = f.readline()
    while(len(line) > 1 and not '[' in line):
        if line[0] == ";":   # skip comment lines
            line = f.readline()
            continue
        fields=line.split()
        attypeid=atomtypes[fields[1]] # map str type to int type
        types_tmp.append(attypeid)
        if len(fields) > 6:
            # this atom has a charge different from its atomtype
            charge_tmp.append(float(fields[6]))
        else:
            #look up default values for this atom type
            charge_tmp.append(atomtypeparams[attypeid]['charge'])
        if len(fields) > 7:
            # also has a special mass
            mass_tmp.append(float(fields[7]))
        else:
            mass_tmp.append(atomtypeparams[attypeid]['mass'])

        use_atomtypeparams.update({attypeid: atomtypeparams[attypeid]})

        line = f.readline()

    # Convert to sigma/epsilon
    if combinationrule == 1:
        for k, v in use_atomtypeparams.iteritems():
            c6, c12 = float(v['sig']), float(v['eps'])
            sig, eps = convertc6c12(c6, c12)
            print '{}, Convert C6({}), C12({}) to sig({}), eps({})'.format(
                k, c6, c12, sig, eps)
            use_atomtypeparams[k]['sig'] = sig
            use_atomtypeparams[k]['eps'] = eps

    # Prepare nonbonded params to contains only that store in atomtypes
    for k, v in nonbondedparams.iteritems():
        if k[0] in use_atomtypeparams and k[1] in use_atomtypeparams:
            use_nonbond_params.update({k: v})
            if combinationrule == 1:
                c6, c12 = float(v['sig']), float(v['eps'])
                sig, eps = convertc6c12(c6, c12)
                print '{}, Convert C6({}), C12({}) to sig({}), eps({})'.format(
                    k, c6, c12, sig, eps)
                use_nonbond_params[k]['sig'] = sig
                use_nonbond_params[k]['eps'] = eps

    # Reset position
    f.seek(pos)

    # extend copies of this molecule
    num_atoms_molecule = len(types_tmp)
    for i in range(num_molecule_copies):
        types.extend(types_tmp)
        charges.extend(charge_tmp)
        masses.extend(mass_tmp)
        molecule_index.extend([res_idx+i]*num_atoms_molecule)
   
    return types, masses, charges, num_atoms_molecule, molecule_index


def storePairs(f, defaults, types, pairtypeparams,
               use_pairtypeparams,
               atomtypeparams, pairs, num_atoms_molecule, num_molecule_copies, molstartindex):
    pairs_tmp = []
    file_pos = f.tell()
    line = f.readlastline()
    fudgeLJ = float(defaults.get('fudgeLJ', 1.0))
    print('Using fudgeLJ: {}'.format(fudgeLJ))
    combinationrule = defaults['combinationrule']
    types_pairtypeid = {}

    line = f.readline().strip()
    in_section = False
    while line:
        if line.startswith('['):
            if 'pairs' in line:
                in_section = True
            else:
                in_section = False
            line = f.readline().strip()
            continue
        elif line.startswith(';'):
            line = f.readline().strip()
            continue
        elif in_section:
            tmp = line.split(';')[0].split()
            lookup = len(tmp) <= 3
            pid1, pid2 = sorted(map(int, tmp[0:2]))
            t1, t2 = sorted([types[pid1-1], types[pid2-1]])
            pairtypeid = max(use_pairtypeparams) + 1 if use_pairtypeparams else 0
            if lookup:
                at1 = atomtypeparams[t1]
                at2 = atomtypeparams[t2]
                if (t1, t2) in pairtypeparams:
                    if types_pairtypeid:
                        pairtypeid = types_pairtypeid.setdefault(
                            (t1, t2),
                            max(types_pairtypeid.values())+1)
                    else:
                        pairtypeid = 0
                        types_pairtypeid[(t1, t2)] = 0
                    use_pairtypeparams[pairtypeid] = pairtypeparams[(t1, t2)]
                else:
                    sig_1, eps_1 = at1['sig'], at1['eps']
                    sig_2, eps_2 = at2['sig'], at2['eps']
                    eps = fudgeLJ*(eps_1*eps_2)**(1.0/2.0)
                    if combinationrule == 2:
                        sig = 0.5*(sig_1 + sig_2)
                    else:
                        sig = (sig_1*sig_2)**(1.0/2.0)
                    pairtypeid = max(use_pairtypeparams) + 1 if use_pairtypeparams else 0
                    use_pairtypeparams[pairtypeid] = {'sig': sig, 'eps': eps}
                    pairtypeparams[(t1, t2)] = use_pairtypeparams[pairtypeid]
                    types_pairtypeid[(t1, t2)] = pairtypeid
                pairs_tmp.append((pid1, pid2, pairtypeid))
            else:  # Params provided
                if int(tmp[2]) != 1:
                    print('Warning! Supported only pair with type 1, given: {}'.format(
                        tmp[2]))
                    line = f.readline()
                    continue
                sig = float(tmp[3])
                eps = float(tmp[4])
                use_pairtypeparams.update({
                    pairtypeid: {
                        'sig': sig,
                        'eps': eps
                        }})
                pairs_tmp.append((pid1, pid2, pairtypeid))
                pairtypeid += 1
            line = f.readline()
        else:
            line = f.readline()

    f.seek(file_pos)

    if combinationrule == 1:
        for k, v in use_pairtypeparams.iteritems():
            c6, c12 = float(v['sig']), float(v['eps'])
            sig, eps = convertc6c12(c6, c12)
            use_pairtypeparams[k]['sig'] = sig
            use_pairtypeparams[k]['eps'] = eps

    # Extend pairs to copies of molecule
    pairs_per_molecule = len(pairs_tmp)
    for i in range(num_molecule_copies):
        for j in range(pairs_per_molecule):
            pid1, pid2, pairtypeid = pairs_tmp[j]
            ia = molstartindex + pid1 + (i * num_atoms_molecule)
            ib = molstartindex + pid2 + (i * num_atoms_molecule)
            if pairtypeid in pairs:
                pairs[pairtypeid].append((ia, ib))
            else:
                pairs.update({pairtypeid: [(ia, ib)]})
    return pairs

def storeBonds(f, types, bondtypes, bondtypeparams, bonds, num_atoms_molecule,\
    num_molecule_copies, molstartindex, exclusions, nregxcl, doRegularExcl=True):
    line = ''
    bonds_tmp = []
    top = False
    pos = f.tell()
    line=f.readlastline()
    local_exclusions=[] # excluded pairs of atoms within this mol (local ids)

    while not 'bonds' in line:
        line = f.readline()
        if 'moleculetype' in line or not line:
            f.seek(pos)
            return bonds

    line = f.readline()
    while(len(line) > 1 and not '[' in line):
        if line[0] == ";":   # skip comment lines
            line = f.readline()
            continue
        tmp = line.split()
        lookup=(len(tmp)<=3) # if the bond has < 3 arguments, it is defined in the bondtypes section and we have to look it up
        pid1, pid2 = map(int, tmp[0:2])
        if lookup:
            # based on atom names: potential has to be defined in bondtypes already
            # this is for tabulated bond potentials specified based on type
            t1, t2 = types[pid1-1], types[pid2-1]
            if t1 > t2: # interactions in the other way
                t1, t2 = t2, t1
            bdtypeid = bondtypes[t1][t2] #bondtypes[t1][t2]
        else:
            # this one is specific for this pair of atoms: check if we need to make a new type
            temptype=ParseBondTypeParam(line)
            bdtypeid=FindType(temptype, bondtypeparams)
            if bdtypeid==None:
                bdtypeid=len(bondtypeparams)
                bondtypeparams.update({bdtypeid:temptype})
            
        bonds_tmp.append((pid1, pid2, bdtypeid)) # store bondtypes for this molecule
        if bondtypeparams[bdtypeid].automaticExclusion():
             # this bond type generates an exclusion as defined by the
             # function type (see gromacs manual)
            local_exclusions.append((pid1, pid2))
        line = f.readline()
    

    if doRegularExcl:
        # generate exclusions for atoms up to a number of nregxcl bonds away
        # see gromacs manual, section 5.4
        exclusions_bonds=[]
        for b in bonds_tmp:
            pid1, pid2, bdtypeid = b[0:3]
            exclusions_bonds.append((pid1, pid2))   
        print "Generating Regular exclusions nregxcl=", nregxcl
        print "Warning: this doesn't work for systems containing (1,5)-bonds, i.e. additional bonds between"
        print "particles which are 4 bonds apart along a chain, e.g. some CG polymer models (see gromacs.py for solution)"
        #for systems with (1,5) bonds, use local_exclusions=GenerateRegularExclusions(local_exclusions, nregxcl,local_exclusions)
        local_exclusions=GenerateRegularExclusions(exclusions_bonds, nregxcl,local_exclusions)
    # extend bonds to copies of this molecule
    bonds_per_mol = len(bonds_tmp)
    for i in range(num_molecule_copies):
        for j in range(bonds_per_mol):
            pid1, pid2, bdtypeid = bonds_tmp[j][0:3]
            ia=molstartindex+pid1 + (i * num_atoms_molecule) # index of copy atom i
            ib=molstartindex+pid2 + (i * num_atoms_molecule) # index of copy atom j
            
            if bdtypeid in bonds:
                bonds[bdtypeid].append((ia, ib))
            else:
                bonds.update({bdtypeid:[(ia, ib)]})
                
                
    # now, extend also the regular exclusions
    for i in range(num_molecule_copies):
        for exclpair in local_exclusions:
            pid1, pid2 = exclpair[0:2]
            ia=molstartindex+pid1 + (i * num_atoms_molecule) # index of copy atom i
            ib=molstartindex+pid2 + (i * num_atoms_molecule) # index of copy atom j
            exclusions.append((ia,ib))
    return bonds
        
def storeAngles(f, types, angletypes, angletypeparams, angles, num_atoms_molecule, num_molecule_copies, molstartindex):
    line = ''
    angles_tmp = []
    pos = f.tell()
    line=f.readlastline()
    while not 'angles' in line:
        line = f.readline()
        if 'moleculetype' in line or not line:
            f.seek(pos)
            return angles
        
    line = f.readline()
    while(len(line) > 1 and not '[' in line):
        if line[0] == ";": # skip comment lines
            line = f.readline()
            continue
        tmp = line.split()
        lookup=(len(tmp)<=4)
        pid1, pid2, pid3 = map(int, tmp[0:3])
        if lookup:
            t1, t2, t3 = types[pid1-1], types[pid2-1], types[pid3-1]
            try:
                antypeid = angletypes[t1][t2][t3]
            except KeyError:
                #todo: is this good style?
                t1, t3 = t3, t1
                antypeid = angletypes[t1][t2][t3]
        else:
            #check if we need to make new type
            temptype=ParseAngleTypeParam(line) 
            antypeid=FindType(temptype, angletypeparams)
            if antypeid==None:
                antypeid=len(angletypeparams)
                angletypeparams.update({antypeid:temptype})
                
        angles_tmp.append((pid1, pid2, pid3, antypeid)) # store angletypes for this molecule
        
        line = f.readline()
        
    # extend angles to copies of this molecule
    angles_per_mol = len(angles_tmp)
    for i in range(num_molecule_copies):
        for j in range(angles_per_mol):
            pid1, pid2, pid3, antypeid = angles_tmp[j][0:4]
            ia=molstartindex+pid1 + (i * num_atoms_molecule) # index of copy atom i
            ib=molstartindex+pid2 + (i * num_atoms_molecule) # index of copy atom j
            ic=molstartindex+pid3 + (i * num_atoms_molecule) # index of copy atom k
            if antypeid in angles:
                angles[antypeid].append((ia, ib, ic))
            else:
                angles.update({antypeid:[(ia, ib, ic)]})
    return angles  


def storeDihedrals(f, types, dihedraltypes, dihedraltypeparams, dihedrals, num_atoms_molecule, num_molecule_copies, molstartindex, atomtype_wildcard):
    line = ''
    dihedrals_tmp = []
    pos = f.tell()
    line=f.readlastline()
    while not 'dihedrals' in line:
        line = f.readline()
        if 'moleculetype' in line or not line:
            f.seek(pos)
            return dihedrals

    line = f.readline()
    while(len(line) > 1 and not '[' in line):
        if line[0] == ";": # skip comment lines
            line = f.readline()
            continue
        tmp = line.split()
        lookup=(len(tmp)<=5)
        pid1, pid2, pid3, pid4 = map(int, tmp[0:4])
        if lookup:
            t1, t2, t3, t4 = types[pid1-1], types[pid2-1], types[pid3-1], types[pid4-1] # get types of particles
            try:
                dihtypeid = dihedraltypes[t1][t2][t3][t4] #dihtypeid is now a tuple
            #if t1 not in dihedraltypes: # interactions in the other way
            except KeyError:
                t1, t2, t3, t4 = t4, t3, t2, t1
                try: 
                    dihtypeid = dihedraltypes[t1][t2][t3][t4]
                except KeyError:
                    t1, t2, t3, t4 = atomtype_wildcard, t2, t3, atomtype_wildcard
                    try:
                        dihtypeid = dihedraltypes[t1][t2][t3][t4]
                    except KeyError:
                        t1, t2, t3, t4 = t1, t3, t2, t4
                        dihtypeid = dihedraltypes[t1][t2][t3][t4]
                #t1, t2, t3, t4 = t4, t1, t2, t3
                #dihtypeid = dihedraltypes[t1][t2][t3][t4]
        else:
            #check if we need to make new type
            temptype=ParseDihedralTypeParam(line)
            dihtypeid=FindType(temptype, dihedraltypeparams) #here,dihtypeid is an int, not a tuple
            if dihtypeid==None:
                dihtypeid=len(dihedraltypeparams)
                dihedraltypeparams.update({dihtypeid:temptype})

            dihtypeid=(dihtypeid,) #convert to tuple for putting in dihedrals_tmp
        
        dihedrals_tmp.append((pid1, pid2, pid3,pid4, dihtypeid)) # 
        line = f.readline()
        
    # extend angles to copies of this molecule
    dihedrals_per_mol = len(dihedrals_tmp)
    for i in range(num_molecule_copies):
        for j in range(dihedrals_per_mol):
            pid1, pid2, pid3, pid4, dihtypeid = dihedrals_tmp[j][0:5]
            ia=molstartindex+pid1 + (i * num_atoms_molecule) # index of copy atom i
            ib=molstartindex+pid2 + (i * num_atoms_molecule) # index of copy atom j
            ic=molstartindex+pid3 + (i * num_atoms_molecule) # index of copy atom k
            id=molstartindex+pid4 + (i * num_atoms_molecule) # index of copy atom l
            if dihtypeid in dihedrals:
                dihedrals[dihtypeid].append((ia, ib, ic, id)) # ###what happens now that it's a tuple? this has only been briefly tested for the case of more than one molecule containing dihedrals
            else:
                dihedrals.update({dihtypeid:[(ia, ib, ic, id)]})
    return dihedrals
    
def storeImpropers(f, types, impropertypes, impropertypeparams, impropers, num_atoms_molecule, num_molecule_copies, molstartindex, atomtype_wildcard):
    line = ''                          
    impropers_tmp = []
    pos = f.tell()
    line=f.readlastline()
    while not 'impropers' in line:     
        line = f.readline()
        if 'moleculetype' in line or not line:
            f.seek(pos)
            return impropers
    
    line = f.readline()
    while(len(line) > 1 and not '[' in line):
        if line[0] == ";": # skip comment lines
            line = f.readline()
            continue
        tmp = line.split()
        lookup=(len(tmp)<=5)
        pid1, pid2, pid3, pid4 = map(int, tmp[0:4])
        if lookup:                                        
            t1, t2, t3, t4 = types[pid1-1], types[pid2-1], types[pid3-1], types[pid4-1] # get types of particles
            try:
                dihtypeid = impropertypes[t1][t2][t3][t4] #dihtypeid is now a tuple
#                print t1, t2, t3, t4, 'found'
            except KeyError:                    
#                print t1, t2, t3, t4, 'not yet found'
                t1, t2, t3, t4 = atomtype_wildcard, t2, t3, t4
                try:
                    dihtypeid = impropertypes[t1][t2][t3][t4]
                except KeyError:
#                    print t1, t2, t3, t4, 'not yet found'
                    t1, t2, t3, t4 = atomtype_wildcard, atomtype_wildcard, t3, t4
                    try:
                        dihtypeid = impropertypes[t1][t2][t3][t4]
                    except KeyError:
#                        print t1, t2, t3, t4, 'not yet found'
                        t1, t2, t3, t4 = atomtype_wildcard, atomtype_wildcard,types[pid1-1], types[pid2-1]
                        try:
                            dihtypeid = impropertypes[t1][t2][t3][t4]
                        except KeyError:
#                            print t1, t2, t3, t4, 'not yet found'
                            t1, t2, t3, t4 = types[pid4-1], types[pid2-1], types[pid3-1], types[pid1-1]
                            try:
                                dihtypeid = impropertypes[t1][t2][t3][t4]
                            except KeyError:
#                                print t1, t2, t3, t4, 'not yet found'
                                t1, t2, t3, t4 = types[pid4-1], types[pid3-1], types[pid2-1], types[pid1-1]
                                try:
                                    dihtypeid = impropertypes[t1][t2][t3][t4]
                                except KeyError:
                                    print t1, t2, t3, t4, 'not yet found in impropers'
                                    quit()
        else:
            #check if we need to make new type
            temptype=ParseImproperTypeParam(line)
            dihtypeid=FindType(temptype, impropertypeparams) #here,dihtypeid is an int, not a tuple
            if dihtypeid==None:
                dihtypeid=len(impropertypeparams)
                impropertypeparams.update({dihtypeid:temptype})

        impropers_tmp.append((pid1, pid2, pid3,pid4, dihtypeid)) # 
        line = f.readline()

    # extend angles to copies of this molecule
    impropers_per_mol = len(impropers_tmp)
    for i in range(num_molecule_copies):
        for j in range(impropers_per_mol):
            pid1, pid2, pid3, pid4, dihtypeid = impropers_tmp[j][0:5]
            ia=molstartindex+pid1 + (i * num_atoms_molecule) # index of copy atom i
            ib=molstartindex+pid2 + (i * num_atoms_molecule) # index of copy atom j
            ic=molstartindex+pid3 + (i * num_atoms_molecule) # index of copy atom k
            id=molstartindex+pid4 + (i * num_atoms_molecule) # index of copy atom l
            if dihtypeid in impropers:
                impropers[dihtypeid].append((ia, ib, ic, id)) # ###what happens now that it's a tuple?
            else:
                impropers.update({dihtypeid:[(ia, ib, ic, id)]})
    return impropers
### adapt for impropers

def setBondedInteractions(system, bonds, bondtypeparams):
    list={}
    bc=0
    for id, bondlist in bonds.iteritems():
        fpl = espressopp.FixedPairList(system.storage)
        fpl.addBonds(bondlist)
        bc+=len(bondlist) 
        bdinteraction=bondtypeparams[id].createEspressoInteraction(system, fpl)
        if bdinteraction:
            system.addInteraction(bdinteraction)
            list.update({id: bdinteraction})
    return list

def setAngleInteractions(system, angles, angletypeparams):
    list={}
    
    for id, anglelist in angles.iteritems():
        fpl = espressopp.FixedTripleList(system.storage)
        fpl.addTriples(anglelist)
        angleinteraction=angletypeparams[id].createEspressoInteraction(system, fpl)
        if angleinteraction:
            system.addInteraction(angleinteraction)
            list.update({id: angleinteraction})
    return list

def setDihedralInteractions(system, dihedrals, dihedraltypeparams):
    d_list={}
    
    for dihid, dihedrallist in dihedrals.iteritems():
        fpl = espressopp.FixedQuadrupleList(system.storage)
        fpl.addQuadruples(dihedrallist)
        dihedralinteraction=dihedraltypeparams[dihid].createEspressoInteraction(system, fpl)
        if dihedralinteraction:
            system.addInteraction(dihedralinteraction)
            ii = len(d_list)
            d_list.update({ii: dihedralinteraction}) #ii instead of id bcs same id may already have been encountered in another idlist (tuple of id's)
    return d_list

def setImproperInteractions(system, impropers, impropertypeparams):
    list={}
    
    for idlist, improperlist in impropers.iteritems():
        fpl = espressopp.FixedQuadrupleList(system.storage)
        fpl.addQuadruples(improperlist)
        for i in range(len(idlist)):
          id=idlist[i]
          improperinteraction=impropertypeparams[id].createEspressoInteraction(system, fpl)
          if improperinteraction:
              system.addInteraction(improperinteraction)
              ii = len(list)
              list.update({ii: improperinteraction}) #ii instead of id bcs same id may already have been encountered in another idlist (tuple of id's)
    return list

def setBondedInteractionsAdress(system, bonds, bondtypeparams,ftpl):
    list={}
    bc=0
    for id, bondlist in bonds.iteritems():
        fpl = espressopp.FixedPairListAdress(system.storage,ftpl)
        fpl.addBonds(bondlist)
        bc+=len(bondlist) 
        bdinteraction=bondtypeparams[id].createEspressoInteraction(system, fpl)
        if bdinteraction:
            system.addInteraction(bdinteraction)
            list.update({id: bdinteraction})
    return list

def setAngleInteractionsAdress(system, angles, angletypeparams,ftpl):
    list={}
    
    for id, anglelist in angles.iteritems():
        fpl = espressopp.FixedTripleListAdress(system.storage,ftpl)
        fpl.addTriples(anglelist)
        angleinteraction=angletypeparams[id].createEspressoInteraction(system, fpl)
        if angleinteraction:
            system.addInteraction(angleinteraction)
            list.update({id: angleinteraction})
    return list

def setDihedralInteractionsAdress(system, dihedrals, dihedraltypeparams,ftpl):
    list={}
    
    for idlist, dihedrallist in dihedrals.iteritems():
        fpl = espressopp.FixedQuadrupleListAdress(system.storage,ftpl)
        fpl.addQuadruples(dihedrallist)
        for i in range(len(idlist)):
          id=idlist[i]
          dihedralinteraction=dihedraltypeparams[id].createEspressoInteraction(system, fpl)
          if dihedralinteraction:
              system.addInteraction(dihedralinteraction)
              ii = len(list)
              list.update({ii: dihedralinteraction}) #ii instead of id bcs same id may already have been encountered in another idlist (tuple of id's)
    return list

def setImproperInteractionsAdress(system, impropers, impropertypeparams,ftpl):
    list={}
    
    for idlist, improperlist in impropers.iteritems():
        fpl = espressopp.FixedQuadrupleListAdress(system.storage,ftpl)
        fpl.addQuadruples(improperlist)
        for i in range(len(idlist)):
          id=idlist[i]
          improperinteraction=impropertypeparams[id].createEspressoInteraction(system, fpl)
          if improperinteraction:
              system.addInteraction(improperinteraction)
              ii = len(list)
              list.update({ii: improperinteraction}) #ii instead of id bcs same id may already have been encountered in another idlist (tuple of id's)
    return list

def setLennardJonesInteractions(system, defaults, atomtypeparams, verletlist, cutoff, 
                                nonbonded_params=None, hadress=False, ftpl=None, 
                                table_groups=None):
    """ Set lennard jones interactions which were read from gromacs based on the atomypes
        
        Args:
            system: The system object.
            defaults: The dictionary with defaults from topology.
            atomtypeparams: The dictionary with atom types parameters.
            verletlist: The VerletList object.
            cutoff: The cutoff value for Lennard Jones potential.
            nonbonded_params: The dictionary with non-bonded parameters.
            hadress: If set to true then H-AdResS is enabled. (default: False)
            ftpl: The FixedTupleList object used by AdResS.
            table_groups: List of atom types that interaction is tabulated.
        Returns:
            interaction object.
    """
    if table_groups is None:
        table_groups = []

    if ftpl: # In AdResS mode
        if hadress:
            interaction=espressopp.interaction.VerletListHadressLennardJones(verletlist, ftpl)
        else:
            interaction=espressopp.interaction.VerletListAdressLennardJones(verletlist, ftpl)
    else:
        interaction = espressopp.interaction.VerletListLennardJones(verletlist)
   
    if nonbonded_params is None:
        nonbonded_params = {}

    combinationrule = defaults['combinationrule']

    type_pairs = sorted({
        tuple(sorted([type_1, type_2])) 
            for type_1, pi in atomtypeparams.iteritems()
            for type_2, pj in atomtypeparams.iteritems()
            if pi['atnum'] not in table_groups and pj['atnum'] not in table_groups 
        })
    print('Number of pairs: {}'.format(len(type_pairs)))
    for type_1, type_2 in type_pairs:
        pi = atomtypeparams[type_1]
        pj = atomtypeparams[type_2]
        if pi['particletype'] == 'V' or pj['particletype'] == 'V':
            print('Skip {}-{}'.format(type_1, type_2))
            continue
        param = nonbonded_params.get((type_1, type_2))
        if param:
            print 'Using defined non-bonded cross params', param
            sig, eps = param['sig'], param['eps']
        else:
            sig_1, eps_1 = float(pi['sig']), float(pi['eps'])
            sig_2, eps_2 = float(pj['sig']), float(pj['eps'])
            if combinationrule == 2:
                sig = 0.5*(sig_1 + sig_2)
                eps = (eps_1*eps_2)**(1.0/2.0)
            else:
                sig = (sig_1*sig_2)**(1.0/2.0)
                eps = (eps_1*eps_2)**(1.0/2.0)
            print 'Combination rule: {}, sig={}, eps={}'.format(combinationrule, sig, eps)
        if sig > 0.0 and eps > 0.0:
            print "Setting LJ interaction for", type_1, type_2, "to sig ", sig, "eps", eps, "cutoff", cutoff
            ljpot = espressopp.interaction.LennardJones(epsilon=eps, sigma=sig, shift='auto', cutoff=cutoff)
            if ftpl:
                interaction.setPotentialAT(type1=type_1, type2=type_2, potential=ljpot)
            else:
                interaction.setPotential(type1=type_1, type2=type_2, potential=ljpot)

    system.addInteraction(interaction)
    return interaction

def setLennardJones14Interactions(system, defaults, atomtypeparams, onefourlist, cutoff):
    """ Set lennard jones interactions which were read from gromacs based on the atomypes"""
    interaction = espressopp.interaction.FixedPairListTypesLennardJones(system,onefourlist)
    
    #for i in range(len(atomtypeparams)):
    #    for j in range(i, len(atomtypeparams)):
    for i in atomtypeparams.keys():
         for j in atomtypeparams.keys():
            pi=atomtypeparams[i]
            pj=atomtypeparams[j]
            if pi!=pj:
                sig=0.5*(float(pi['sig'])+float(pj['sig']))
                eps=math.sqrt(float(pi['eps'])*float(pj['eps']))
            else:
                sig=float(pi['sig'])
                eps=float(pi['eps'])
            if (sig>0 and eps >0):
                eps = eps*fudge
                #print "Setting 1-4 LJ interaction for", i, j, "to sig ", sig, "eps", eps, "cutoff", cutoff
                ljpot= espressopp.interaction.LennardJones(epsilon=eps, sigma=sig, cutoff=cutoff, shift=0)
                #ljpot= espressopp.interaction.LennardJonesGromacs(epsilon=eps, sigma=sig, cutoff=cutoff, shift=0)
                interaction.setPotential(type1=i, type2=j, potential=ljpot)
    system.addInteraction(interaction)
    return interaction

def setCoulombInteractions(system, verletlist, rc, types, epsilon1, epsilon2,kappa, hadress=False, adress=False, ftpl=None):
 
    print "# Setting up Coulomb reaction field interactions"

    pref=138.935485 # we want gromacs units, so this is 1/(4 pi eps_0) ins units of kJ mol^-1 e^-2
    
    pot = espressopp.interaction.ReactionFieldGeneralized(prefactor=pref, kappa=kappa, epsilon1=epsilon1, epsilon2=epsilon2, cutoff=rc)
    #pot = espressopp.interaction.CoulombTruncated(prefactor=pref, cutoff=rc)
    if (hadress and adress):
      print "Error! In gromacs.setCoulombInteractions, you cannot use adress and hadress at the same time"
      return
    if (hadress):
        interaction=espressopp.interaction.VerletListHadressReactionFieldGeneralized(verletlist, ftpl)
    elif (adress):
        interaction=espressopp.interaction.VerletListAdressReactionFieldGeneralized(verletlist, ftpl)
    else: 
	interaction=espressopp.interaction.VerletListReactionFieldGeneralized(verletlist)
   # interaction=espressopp.interaction.VerletListCoulombTruncated(verletlist)
            
    for i in range(max(types)+1):
        for k in range(i, max(types)+1):
	    if (hadress or adress):
		interaction.setPotentialAT(type1=i, type2=k, potential=pot)
	    else:
		interaction.setPotential(type1=i, type2=k, potential=pot)

    system.addInteraction(interaction)
    return interaction

def setCoulombInteractionsProtein(system, verletlist, rc, types, epsilon1, epsilonprot,epsilonwat,kappa,otype,htype, hadress=False, adress=False, ftpl=None):
 
    print "# Setting up Coulomb reaction field interactions"
    print "# Using ",epsilonwat," for water and wat-prot and ",epsilonprot," for protein"

    pref=138.935485 # we want gromacs units, so this is 1/(4 pi eps_0) ins units of kJ mol^-1 e^-2
    
    potwat = espressopp.interaction.ReactionFieldGeneralized(prefactor=pref, kappa=kappa, epsilon1=epsilon1, epsilon2=epsilonwat, cutoff=rc)
    potprot = espressopp.interaction.ReactionFieldGeneralized(prefactor=pref, kappa=kappa, epsilon1=epsilon1, epsilon2=epsilonprot, cutoff=rc)
    #pot = espressopp.interaction.CoulombTruncated(prefactor=pref, cutoff=rc)
    if (hadress and adress):
      print "Error! In gromacs.setCoulombInteractions, you cannot use adress and hadress at the same time"
      return
    if (hadress):
        interaction=espressopp.interaction.VerletListHadressReactionFieldGeneralized(verletlist, ftpl)
    elif (adress):
        interaction=espressopp.interaction.VerletListAdressReactionFieldGeneralized(verletlist, ftpl)
    else: 
	interaction=espressopp.interaction.VerletListReactionFieldGeneralized(verletlist)
   # interaction=espressopp.interaction.VerletListCoulombTruncated(verletlist)
            
    for i in range(max(types)+1):
        for k in range(i, max(types)+1):
            if (i==otype or i==htype or k==otype or k==htype):
	        if (hadress or adress):
	            interaction.setPotentialAT(type1=i, type2=k, potential=potwat)
	        else:
	            interaction.setPotential(type1=i, type2=k, potential=potwat)
            else: 
	        if (hadress or adress):
	            interaction.setPotentialAT(type1=i, type2=k, potential=potprot)
	        else:
	            interaction.setPotential(type1=i, type2=k, potential=potprot)

    system.addInteraction(interaction)
    return interaction

def setCoulomb14Interactions(system, defaults, onefourlist, rc, types): 

    #in Gromas, 1-4 interactions don't have reaction field correction
    print "# Setting up 1-4 Coulomb interactions"

    if defaults:
        fudge=float(defaults['fudgeQQ'])
        print "# Using electrostatics 1-4 fudge factor ",fudge

    pref=138.935485*fudge # we want gromacs units, so this is 1/(4 pi eps_0) ins units of kJ mol^-1 e^-2, scaled by fudge factor

    #pot = espressopp.interaction.CoulombRSpace(prefactor=pref, alpha=0.0, cutoff=rc)
    pot = espressopp.interaction.CoulombTruncated(prefactor=pref, cutoff=rc)

    #interaction=espressopp.interaction.FixedPairListTypesCoulombRSpace(system,onefourlist)
    interaction=espressopp.interaction.FixedPairListTypesCoulombTruncated(system,onefourlist)

    for i in range(max(types)+1):
        for k in range(i, max(types)+1):
            interaction.setPotential(type1=i, type2=k, potential=pot)

    system.addInteraction(interaction)
    return interaction

def setTabulatedInteractions(system, atomtypeparams, vl, cutoff=None, interaction=None, ftpl=None,
                             table_groups=None, adress_is_cg=True, hadress=False,
                             atomtypes_group=None,
                             spline_type=3):
    """Sets tabulated potentials.

    Args:
        system: The system object.
        atomtypeparams: Dictionary with key as atom type and values with atomtypes from topology.
        vl: The VerletList object.
        cutoff: The cutoff value for potentials.
        interaction: The interaction object.
        ftpl: The fixed tuple list object for AdResS.
        table_groups: The list of atom type names that should be use for tabulated potentials. If
            empty then all atom types will be consider to use tabulated potential.
        atomtypes_group: The dictionary that groups number of atomtypes into single type. The
            key is an original atom type and the value is corresponding group.
        spline_types: The type of sp-line. (default: 3)
    """
    if table_groups is None:
        table_groups = []
    if atomtypes_group is None:
        atomtypes_group = {}

    if interaction is None:
        if ftpl:
            interaction = espressopp.interaction.VerletListAdressTabulated(vl, ftpl)
        else:
            interaction = espressopp.interaction.VerletListTabulated(vl)
    
    name_pairs = set()
    name_types = defaultdict(set)
    for type_1, v1 in atomtypeparams.iteritems():
        for type_2, v2 in atomtypeparams.iteritems():
            t1 = atomtypes_group.get(v1['atnum'], v1['atnum'])
            t2 = atomtypes_group.get(v2['atnum'], v2['atnum'])
            if not table_groups or (t1 in table_groups and t2 in table_groups):
                kk = tuple(sorted([t1, t2]))
                tt = tuple(sorted([type_1, type_2]))
                name_pairs.add(kk)
                name_types[kk].add(tt)
    
    for name_1, name_2 in name_pairs:
        print('Set tabulated potential {}-{}'.format(name_1, name_2))
        table_name = '{}-{}.tab'.format(name_1, name_2)
        orig_table_name = 'table_{}_{}.xvg'.format(name_1, name_2)
        print('Converting table_{name1}_{name2}.xvg to {name1}-{name2}.tab'.format(
            name1=name_1, name2=name_2))
        convertTable(orig_table_name, table_name)
        if ftpl:
            if adress_is_cg:
                set_potential_fn = interaction.setPotentialCG
            else:
                set_potential_fn = interaction.setPotentialAT
        else:
            set_potential_fn = interaction.setPotential
        for type_1, type_2 in name_types[(name_1, name_2)]:
            print('Defining potential for types: {}-{}'.format(type_1, type_2))
            set_potential_fn(
                type1=type_1,
                type2=type_2,
                potential=espressopp.interaction.Tabulated(
                    itype=spline_type,
                    filename=table_name,
                    cutoff=cutoff))

    system.addInteraction(interaction)

    return interaction


def convertTable(gro_in_file, esp_out_file, sigma=1.0, epsilon=1.0, c6=1.0, c12=1.0):
    """Convert GROMACS tabulated file into ESPResSo++ tabulated file (new file
    is created). First column of input file can be either distance or angle.
    For non-bonded files, c6 and c12 can be provided. Default value for sigma, epsilon,
    c6 and c12 is 1.0. Electrostatics are not taken into account (f and fd columns).
    
    Keyword arguments:
    gro_in_file -- the GROMACS tabulated file name (bonded, nonbonded, angle
    or dihedral).
    esp_out_file -- filename of the ESPResSo++ tabulated file to be written.
    sigma -- optional, depending on whether you want to convert units or not.
    epsilon -- optional, depending on whether you want to convert units or not.
    c6 -- optional
    c12 -- optional
    """

    

    # determine file type
    bonded, angle, dihedral = False, False, False
    if gro_in_file[6] == "b":
        bonded = True
    if gro_in_file[6] == "a":
        angle  = True
        bonded = True
    if gro_in_file[6] == "d":
        dihedral = True
        bonded = True

    fin = open(gro_in_file, 'r')
    fout= open(esp_out_file, 'w')

    if bonded: # bonded has 3 columns
        for line in fin:
            if line[0] == "#": # skip comment lines
                continue
            
            columns = line.split()
            r = float(columns[0])
            f = float(columns[1]) # energy
            fd= float(columns[2]) # force
            
            # convert units
            if angle or dihedral: # degrees to radians
                r = math.radians(r)
                fd=fd*180/math.pi
            else:
                r = r / sigma
            e = f / epsilon
            f = fd*sigma / epsilon
            
            if (not angle and not dihedral and r != 0) or \
                 (angle and r <= math.pi and r >= 0) or \
                  (dihedral and r >= -math.pi and r <= math.pi):
                fout.write("%15.8g %15.8g %15.8g\n" % (r, e, f))
    
    else: # non-bonded has 7 columns
        for line in fin:
            if line[0] == "#": # skip comment lines
                continue
            
            columns = line.split()
            r = float(columns[0])
            #f = float(columns[1]) # electrostatics not implemented yet
            #fd= float(columns[2]) # electrostatics not implemented yet
            g = float(columns[3]) # dispersion
            gd= float(columns[4])
            h = float(columns[5]) # repulsion
            hd= float(columns[6])
            
            e = c6*g + c12*h
            f = c6*gd+ c12*hd
            
            # convert units
            r = r / sigma
            e = e / epsilon
            f = f*sigma / epsilon
            
            if r != 0: # skip 0
                fout.write("%15.8g %15.8g %15.8g\n" % (r, e, f))
    
    fin.close()
    fout.close()
