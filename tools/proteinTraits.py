'''
Created on 18/02/2014

@author: victor
'''
import prody
import numpy
import Bio.pairwise2
import Bio.SubsMat
import Levenshtein # pip install python-Levenshtein
import scipy.linalg
import math

def getPDB(pdbId):
    """
    Downloads a pdb from the Protein Data Bank (if necessary) and removes all models so that it only has one
    model.

    @param pdbId: A 4 letter pdb id

    @return: The downloaded pdb data structure.
    """
    # Download pdb
    path = prody.fetchPDB(pdbId, compressed = False)
    # Get pdb data structure
    pdb = prody.parsePDB(path)
    number_of_models = pdb.select("protein").numCoordsets()
    # Delete all coordsets but coordset 0
    [pdb.delCoordset(1) for i in range(1,number_of_models)]
    return pdb

def savePDB(pdb, filename):
    """
    Saves a prody pdb data structure in a pdb format file (simple wrapper).

    @param pdb: A prody pdb data structure.
    """
    prody.writePDB(filename, pdb)

def getNumberOfAtoms(pdb):
    """
    Prody wrapper, returns the number of atoms of a pdb.

    @param pdb: A prody pdb data structure.

    @return: The number of atoms of that protein.
    """
    return pdb.numAtoms

def squared_distance(point1,point2):
    """
    Calculates the distance between 2 points.

    @param point1: Tridimensional point.
    @param point2: Another tridimensional point.

    @return: The squared distance between this two points.
    """
    return (point1[0]-point2[0])**2+(point1[1]-point2[1])**2+(point1[2]-point2[2])**2

def getAverageNeighbourCount(pdb, cutoffDistance):
    """
    Calculates the average value of the number of neighbours of each atom for a given cutoff radius.

    @param pdb: A prody pdb data structure.

    @param cutoffDistance: the cutoff radius.

    @return: The average value of the number of neighbours of each atom.

    """
    sq_cutoff_distance = cutoffDistance* cutoffDistance
    coordsets = pdb.getCoordsets()
    # Must have only one coordset
    num_neigs = [0]*pdb.numAtoms
    for i in range(pdb.numAtoms -1):
        atom_i = coordsets[i]
        for j in range(i+1, pdb.numAtoms):
            if squared_distance(atom_i, coordsets[j]) < sq_cutoff_distance:
                num_neigs[i] += 1
                num_neigs[j] += 1
    return numpy.mean(num_neigs), numpy.std(num_neigs)

def getProteinSequence(pdb):
    """
    Generates the 1 letter per residue sequence for a protein. Uses a dictionary that maps the 3 letter naming with the 1 letter naming convention
    Source:
        - Biskit (http://biskit.pasteur.fr/)

    @param pdb: A prody pdb data structure.

    @return: A string with the sequence of this protein.
    """
    aaDicStandard ={   'asp':'D', 'glu':'E', 'lys':'K', 'his':'H', 'arg':'R',
                       'gln':'Q', 'asn':'N', 'ser':'S', 'asx':'B', 'glx':'Z',
                       'phe':'F', 'trp':'W', 'tyr':'Y','gly':'G', 'ala':'A',
                       'ile':'I', 'leu':'L', 'cys':'C', 'met':'M', 'thr':'T',
                       'val':'V', 'pro':'P' ,'cyx':'C', 'hid':'H', 'hie':'H',
                       'hip':'H', 'unk':'X', 'ace':'X', 'nme':'X'}

    # One-liner just for the sake of the challenge
    return "".join([aaDicStandard[resname] if resname in aaDicStandard else "X" for resname in [residue.getResname().lower() for residue in prody.HierView(pdb).iterResidues()]])

def getSequenceIdentity(seq1,seq2):
    """
    Calculates the best identity score of the aligned sequences.

    @param seq1: A string representing one protein sequence.
    @param seq2: A string representing the other protein sequence.

    @return: The identity score [0,1]
    """
    alignments = Bio.pairwise2.align.globalds(seq1, seq2, Bio.SubsMat.MatrixInfo.blosum62, -10, -0.5)
    scores = [Levenshtein.ratio(s1,s2) for (s1,s2,sc1,sc2,sc3) in alignments]
    return numpy.max(scores)

def getOBBAxisLengths(pdb):
    """
    Calculates the length of the sides of the Oriented Bounding Box that encloses the atoms of the conformation. Only CA,CB and CG
    atoms are used in order to decrease the number of variables.

    @param pdb: A prody pdb data structure.

    @return: The sorted lengths of the OBB edges.
    """
    coordset = pdb.select("protein and not hydrogen").getCoordsets()[0] # first model
    axis = pca(coordset)
    return project_points_and_get_lengths_in_axis(coordset, axis)

def normalize3D(vector):
    """
    Normalizes a tridimensional vector.

    @param vector: The 3D vector we want to normalize.

    @results: The normalized vector.
    """
    return numpy.array(vector) / math.sqrt(numpy.dot(vector,vector))

def project_points_and_get_lengths_in_axis(points, axis):
    """
    Calculates the length of the projection area of the studied points over each of the axis.

    @param points: An array of points from which we want to extract the PCA.
    @param points: The points of the dataset we want to project in that axis.

    @return: The lengths of the projected range for each axis.
    """
    lengths = []

    # Normalize all axis
    for ax in axis:
        norm_axis = normalize3D(ax) # However eigenvectors come normalized from scipy.linalg.eig
        one_d_points = numpy.dot(points, norm_axis)
        lengths.append(numpy.max(one_d_points) - numpy.min(one_d_points))

    return sorted(lengths, reverse=True)

def pca(points):
    """
    Simple implementation of PCA analysis.
    Sources:
        - http://stackoverflow.com/questions/13224362/pca-analysis-with-python, answer 2.
        - http://hewjunwei.wordpress.com/2013/01/26/obb-generation-via-principal-component-analysis/

    @param points: An array of points from which we want to extract the PCA.

    @result: the ordered eigenvectors (first is the axis of major variance, last is the axis of less variance).
    """
    mn = numpy.mean(points, axis=0)
    points -= mn
    C = numpy.cov(points.T)
    evals, evecs = scipy.linalg.eig(C, overwrite_a = True, check_finite= False)
    idx = numpy.argsort(evals)[::-1]
    evecs = evecs[:,idx].T #Traspose is correct as eigenvector[i] = evecs[:,i]
    return evecs


