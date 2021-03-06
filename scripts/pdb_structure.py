import math
import numpy as np

class Atom:
    def __init__(self, atom_string):
        self.id = int(atom_string[6:11])
        self.name = atom_string[11:16].strip()
        self.alt = atom_string[16]
        self.resn = atom_string[17:20].strip()
        self.chain = atom_string[21]
        self.resi = int(atom_string[22:26])
        self.x = float(atom_string[30:38])
        self.y = float(atom_string[38:46])
        self.z = float(atom_string[46:54])
        self.pos = np.array([self.x,self.y,self.z])
        self.occ = float(atom_string[54:60])
        self.temp_factor = float(atom_string[60:66])
        if len(atom_string)>=78:
            self.elem = atom_string[76:78].strip()

    def __str__(self):
        return self.name+"_"+str(self.resi)

    def __repr__(self):
        return 'Atom("ATOM  %5d%5s %3s %c%4d    %8.3f%8.3f%8.3f")' % (self.id, self.name, self.resn, self.chain, self.resi, self.x, self.y, self.z)

    def pdbString(self):
        return 'ATOM  %5d%5s %3s %c%4d    %8.3f%8.3f%8.3f' % (self.id, self.name, self.resn, self.chain, self.resi, self.x, self.y, self.z)

    def __sub__(self, other):
        return np.array([self.x-other.x, self.y-other.y, self.z-other.z])

    def __add__(self, other):
        return np.array([self.x+other.x, self.y+other.y, self.z+other.z])

    def distance(self, a):
        return math.sqrt((self.x-a.x)**2 + (self.y-a.y)**2 + (self.z-a.z)**2)

    def distanceSquared(self, a):
        return (self.x-a.x)**2 + (self.y-a.y)**2 + (self.z-a.z)**2


class PDBFile:
    """ A representation of a PDB-file """

    def pdbDownload(self,pdb_id):
        hostname="ftp.wwpdb.org"
        directory="/pub/pdb/data/structures/all/pdb/"
        prefix="pdb"
        suffix=".ent.gz"

        import os, sys, ftplib, shutil, gzip

        # Log into server
        #print "Downloading %s from %s ..." % (pdb_id, hostname)
        ftp = ftplib.FTP()
        ftp.connect(hostname)
        ftp.login()

        # Download all files in file_list
        to_get = "%s/%s%s%s" % (directory,prefix,pdb_id.lower(),suffix)
        to_write = "%s%s" % (pdb_id,suffix)
        final_name = "%s.pdb" % to_write[:to_write.index(".")]
        try:
            ftp.retrbinary("RETR %s" % to_get,open(to_write,"wb").write)
            f = gzip.open(to_write,'r')
            g = open(final_name,'w')
            g.writelines(f.readlines())
            f.close()
            g.close()
            os.remove(to_write)
        except ftplib.error_perm:
            os.remove(to_write)
            print("ERROR! %s could not be retrieved from PDB!" % to_get)
            ftp.quit()
            return None

        # Log out
        ftp.quit()
        return final_name


    def __init__(self, pdb):
        """ Initialize a PDBFile object using a PDB-file or PDB-id. If pdb is 4 characters long 
        its assumed to be a PDB-id and the corresponding PDB-file will be downloaded and used. """

        self.file_name = pdb
        self.models = []
        cur_model = None

        if len(pdb)==4: 
            self.file_name = self.pdbDownload(pdb)

        if self.file_name.endswith(".gz"):
            import gzip
            f = gzip.open(self.file_name, "r")
            lines = map(lambda l: l.decode('ascii'), f.readlines())
        else:
            f = open(self.file_name,'r')
            lines = f.readlines()

        for line in lines:
            if line[0:4] == "ATOM":
                    if cur_model==None: cur_model = []
                    cur_model.append(Atom(line))
            if (line[0:6] == "ENDMDL" or line[0:5] == "MODEL") and cur_model!=None:
                    self.models.append(cur_model)
                    cur_model = None
        if cur_model!=None:
                self.models.append(cur_model)

        f.close()

    def removeResidues(self, residues):
        for model in range(len(self.models)):
            self.models[model] = [ a for a in self.models[model] if not a.resi in residues ]


    def getAtom(self, res_number, atom_name, model_number = 0):
        for atom in self.models[model_number]:
            if atom.resi==res_number and atom.name==atom_name:
                return atom

    def getAtomById(self, atom_id):
        for model in self.models:
            for atom in model:
                if atom.id==atom_id:
                    return atom

    def getAtoms(self, model_number = 0):
        return self.models[model_number]

    def getAtomsInResi(self, resi, model_number = 0):
        ret = []
        for atom in self.models[model_number]:
            if atom.resi==resi:
                ret.append(atom)
        return ret

    def getResidues(self, model_number = 0):
        '''
        Return a sorted list of all residue numbers in this structure
        '''
        ret = set()
        for atom in self.models[model_number]:
            ret.add(atom.resi)
        return sorted(list(ret))

    def getResidueIDsandNames(self, model_number = 0):
        '''
        Return a sorted list of all residue numbers and names in this structure
        '''
        ret = set()
        for atom in self.models[model_number]:
            ret.add((atom.resi,atom.resn))
        return sorted(list(ret))

    def getChains(self, model_number = 0):
        '''
        Return a set of unique chain identifiers
        '''
        return set(map(lambda a: a.chain, self.models[model_number]))
    
    def getSequence(self, model_number = 0):
        '''
        Get the sequence of this structure. Currently only works for RNA (single-char resn)
        '''

        protresnmap={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','ASX':'B','CYS':'C','GLU':'E','GLN':'Q','GLX':'Z','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

        ret = ''
        prevResi = -1000
        for atom in self.models[model_number]:
            if prevResi==-1000 or prevResi+1==atom.resi:
                if len(atom.resn)==1: #Its RNA/DNA: Just use the resn
                    ret+=atom.resn
                else: #Its probably protein. Use the resn map
                    if atom.resn in protresnmap:
                        ret+=protresnmap[atom.resn]
                    else:
                        ret+='?'
                prevResi=atom.resi
            else:
                while prevResi<atom.resi:
                    ret+='_'
                    prevResi+=1
        return ret


    def bFactorList(self, model_number = 0, names=["C4'","CA"],resis=None):
        """
        Get a list of b-factors number of atoms with one of the specified names, an optional list of residues
        that can limit b factors to specified residues only, useful for
        comparison across non-identical sequences or with missing loops.
        """
        ret =[]
        for atom in self.models[model_number]:
            if names and atom.name not in names: continue
            if resis and atom.resi not in resis: continue
            ret.append(atom.temp_factor)
        return ret

    def coordMatrix(self, model_number = 0, names=["C4'","CA"],resis=None):
        """
        Get a coordinate-matrix of shape (a, 3) where a is the number of atoms with one of the specified names.
        New: an optional list of residues can limit the coordinate Matrix to specified residues only, useful for
        comparison across non-identical sequences or with missing loops.
        """
        ret = np.zeros(  shape=( len(self.models[model_number]) , 3 )  )
        a = 0
        for atom in self.models[model_number]:
            if names and atom.name not in names: continue
            if resis and atom.resi not in resis: continue
            ret[a][0] = atom.x
            ret[a][1] = atom.y
            ret[a][2] = atom.z
            a+=1
        ret.resize(a,3)
        return ret


    def rmsd(self, pdbFile, model1=0, model2=0, names=None):
        """
        Return the smallest root-mean-square-deviation between coordinates in self and pdbFile.
        If names is None, then all atoms are used.
        """

        crds1 = self.coordMatrix(model1, names=names)
        crds2 = pdbFile.coordMatrix(model2, names=names)
        assert(crds1.shape[1] == 3)
        if crds1.shape[0] != crds2.shape[0]:
            print("Structure 1 size does not match structure 2 (", crds1.shape[0], "vs", crds2.shape[0], ")")
            assert(crds1.shape == crds2.shape)
        n = np.shape(crds1)[0]

        # Move crds1 to origo
        avg1 = np.zeros(3)
        for c1 in crds1:
            avg1 += c1
        avg1 /= n
        for c1 in crds1:
            c1 -= avg1

        # Move crds2 to origo
        avg2 = np.zeros(3)
        for c2 in crds2:
            avg2 += c2
        avg2 /= n
        for c2 in crds2:
            c2 -= avg2

        # Get optimal rotation
        # From http://boscoh.com/protein/rmsd-root-mean-square-deviation.html
        correlation_matrix = np.dot(np.transpose(crds1), crds2)
        u, s, v_tr = np.linalg.svd(correlation_matrix)
        r = np.dot(u, v_tr)
        r_det = np.linalg.det(r)
        if r_det < 0:
            # print 'WARNING: MIRRORING'
            u[:,-1] = -u[:,-1]
            r = np.dot(u, v_tr)
        # is_reflection = (np.linalg.det(u) * np.linalg.det(v_tr))
        # if is_reflection:
        # 	u[:,-1] = -u[:,-1]

        # Apply rotation and find rmsd
        import itertools
        rms = 0.
        for c1, c2 in itertools.izip(crds1, crds2):
            c2_r = np.dot(r, c2) #rotate centroid-shifted coordinates
            #compute optimal RMSD
            tmp = c1-c2_r
            rms += np.dot(tmp, tmp)
        
        # Return RMSD
        return math.sqrt(rms/n)


    def alignTo(self, pdbFile, model1=0, model2=0):
        """
        Return the smallest root-mean-square-deviation between coordinates in self and pdbFile
        and moves self to best align with pdbFile
        """

        crds2 = self.coordMatrix(model1, names=None)
        crds1 = pdbFile.coordMatrix(model2, names=None)

        assert(crds1.shape[1] == 3)

        if crds1.shape[0]!=crds2.shape[0]:
            print("Structure 1 size does not match structure 2 (",crds1.shape[0],"vs",crds2.shape[0],")")
            #assert(crds1.shape == crds2.shape)
        n = np.shape(crds1)[0]


        #Move crds1 to origo
        avg1 = np.zeros(3)
        for c1 in crds1:
            avg1 += c1
        avg1 /= n
        for c1 in crds1:
            c1-=avg1

        #Move crds2 to origo
        avg2 = np.zeros(3)
        for c2 in crds2:
            avg2 += c2
        avg2 /= n
        for c2 in crds2:
            c2-=avg2


        #Get optimal rotation
        #From http://boscoh.com/protein/rmsd-root-mean-square-deviation.html
        correlation_matrix = np.matmul(np.transpose(crds1), crds2)
        u, s, v_tr = np.linalg.svd(correlation_matrix)
        # r = np.matmul(u,v_tr)
        # r_det = np.linalg.det(r)
        r_det = (np.linalg.det(u) * np.linalg.det(v_tr))
        if r_det<0.0:
            #print 'WARNING: MIRRORING'
            u[:,-1] = -u[:,-1]
        r = np.dot(u,v_tr)
        #is_reflection = (np.linalg.det(u) * np.linalg.det(v_tr)) 
        #if is_reflection:
        #	u[:,-1] = -u[:,-1]

        #Apply rotation and find rmsd
        a = 0
        import itertools
        rms = 0.
        for c1, c2 in itertools.izip(crds1, crds2):
            c2_r = np.dot(r, c2) #rotate coordinates
            #RMSD
            tmp = c1-c2_r
            rms += np.dot(tmp, tmp)
            
            c2_r+= avg1 #apply centroid translation of first point cloud to have best overlap
            #Save new positions
            self.models[model1][a].x = c2_r[0]
            self.models[model1][a].y = c2_r[1]
            self.models[model1][a].z = c2_r[2]
            self.models[model1][a].pos = np.array([c2_r[0],c2_r[1],c2_r[2]])
            a+=1
        
        #Return RMSD    
        return math.sqrt(rms/n)

    def save(self, fileName):
        f = open(fileName, "w")
        for model in range(len(self.models)):
            if len(self.models)>1:
                f.write("MODEL "+str(model+1)+"\n")

            for atom in self.models[model]:
                f.write(atom.pdbString()+"\n")

            if len(self.models)>1:
                f.write("ENDMDL\n")
        f.close()

    def __repr__(self):
        return "PDBFile('"+self.file_name+"')"

