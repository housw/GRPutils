#!/usr/bin/python


import sys, os, re, urllib
import subprocess
from Bio import SeqIO
import numpy as np
from collections import Counter
from itertools import izip
import collections
import copy




################### gff record and parser  ####################################

class GffRecord(object):
    """ This class used to represent gff3 record
    """

    __slots__ = ["seqid", "source", "feature", "start", "end", "score", "strand",
                 "phase", "attributes"]

    def __init__(self, seqid, source, feature, start, end, score, strand,
                 phase, attributes ):
        self.seqid = seqid
        self.source = source
        self.feature = feature
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attributes

    def __str__(self):
        """ this function will be used to give string representation of gff record
        """
        pretty_str = ""
        pretty_str = pretty_str + self.seqid +"\t"+ self.source +"\t"+\
                     self.feature +"\t"+ str(self.start) +"\t"+ str(self.end)\
                     +"\t"+ str(self.score) +"\t"+ self.strand +"\t" + \
                     str(self.phase)+ "\t"

        # add attributes info
        attributes_str_list = []
        for k, v in self.attributes.iteritems():
            attributes_str_list.append(k+"="+v)

        pretty_str += ";".join(attributes_str_list)

        return pretty_str





class Gene(GffRecord):

    __slots__ = ["seqid", "source", "feature", "start", "end", "score", "strand",
                 "phase", "attributes", "ID", "Name", "locus_tag", "children",
                 "product"]

    def __init__(self, seqid, source, feature, start, end, score, strand,
                 phase, attributes):
        GffRecord.__init__(self, seqid, source, feature, start, end, score, strand,
                 phase, attributes)
        self.ID = self.attributes['ID']
        self.locus_tag = self.attributes['locus_tag']
        self.Name = self.attributes["Name"] if "Name" in self.attributes else self.locus_tag
        self.children = []
        self.product = None


    def add_child(self, child):
        if isinstance(child, GffRecord):
            self.children.append(child)
            self._update_product()
        else:
            print "child should be GffRecord instance !"


    def _update_product(self):
        if self.product == None:
            product = self.children[0].attributes.get("product", None)
            self.product = product
        else:
            pass






class GffParser(object):
    """ This class used to parse gff3 file
    """

    def __init__(self, handle_or_fileStr):
        self.handle = handle_or_fileStr
        self.genome_info = None

    def _parseAttributes(self, attr_str):
        """Parse the GFF3 attribute column and return a dict
        """
        if attr_str == ".": return {}

        ret = {}

        for attribute in attr_str.split(";"):
            #print attribute
            key, value = attribute.split("=")
            ret[urllib.unquote(key)] = urllib.unquote(value)
        return ret

    def _line_list_parser(self, gff_line_list):
        seqid = gff_line_list[0]
        source = gff_line_list[1]
        feature = gff_line_list[2]
        start = int(gff_line_list[3])
        end = int(gff_line_list[4])
        score = gff_line_list[5]
        strand = gff_line_list[6]
        phase = gff_line_list[7]
        attributes = self._parseAttributes(gff_line_list[8])

        # initialize gff_record
        if not feature.lower() =="gene":
            gff_record = GffRecord(seqid, source, feature, start, end, score,
                                    strand, phase, attributes)
        # initialize gene record
        else:
            #print attributes
            gff_record = Gene(seqid, source, feature, start, end, score,
                                    strand, phase, attributes)
        return gff_record


    def _parse_block(self):
        # judge file opened or not
        if not hasattr(self.handle, "read"):
            handle = open(self.handle, "r")
        else:
            handle = self.handle

        record_lines = []

        while True:
            line = handle.readline()

            if not line:
                # close file handle
                try:
                    handle.close()
                except Exception as e:
                    print e
                break

            else:
                if line.startswith("#"):
                    continue
                else:
                    line_list = line.strip().split("\t")
                    #print line_list
                    # parse region
                    if line_list[2].lower() == "region":
                        self.genome_info = self._line_list_parser(line_list)
                        continue

                    # parse gene
                    elif line_list[2].lower() == "gene":
                        # record_lines already contains gene's info
                        if record_lines:
                            yield self._gene_parser(record_lines)
                            record_lines = [line_list]
                        # record_lines is empty, it's the first gene record in gff
                        else:
                            record_lines.append(line_list)

                    # append other info into gene
                    else:
                        record_lines.append(line_list)



    def _gene_parser(self, list_line_list):
        assert list_line_list[0][2] == "gene", "block line not starts with gene !"
        gene = self._line_list_parser(list_line_list[0])

        for line_list in list_line_list[1:]:
            gff_record = self._line_list_parser(line_list)
            gene.add_child(gff_record)
        return gene


    def __iter__(self):
        return self._parse_block()











########################  input data ##########################################

def read_grp(path2grp):
    """ This function used to read in grp file, input file format like follows

        # # BASE    fwd_coverage    fwd_tss    rev_coverage    rev_tss
        # # colour  255:0:0 128:0:0 0:0:255 0:0:128
        #
        0.00    0.00    0.00    -0.00   -0.00
        1.00    0.00    0.00    -0.00   -0.00
        2.00    0.00    0.00    -0.00   -0.00
        3.00    0.00    0.00    -0.00   -0.00

        and will return a fwd_cov, fwd_tss, rev_cov, rev_tss array object
    """
    fwd_cov = []
    fwd_tss = []
    rev_cov = []
    rev_tss = []

    with open(path2grp, "r") as ih:
        for line in ih:
            if line.startswith("#"):
                continue
            line = line.strip().split()
            fwd_cov.append(float(line[1]))
            fwd_tss.append(float(line[2]))
            rev_cov.append(abs(float(line[3])))
            rev_tss.append(abs(float(line[4])))
    fwd_cov = np.array(fwd_cov)
    fwd_tss = np.array(fwd_tss)
    rev_cov = np.array(rev_cov)
    rev_tss = np.array(rev_tss)

    return fwd_cov, fwd_tss, rev_cov, rev_tss






def get_gene_dict_from_gff(gff_file):
    """ this function used to generate fwd_gene_dict and rev_gene_dict from gff
    """
    fwd_gene_dict = {} # data format {pos: fwd_gene}
    rev_gene_dict = {}
    fwd_gene_starts_dict = {} # data format {pos: fwd_gene}
    rev_gene_starts_dict = {}
    # read gff file and store gene info
    gff_genes = GffParser(gff_file)
    for gene in gff_genes:
        if gene.strand == "+":
            fwd_gene_starts_dict.update({gene.start-1:gene})
            for x in xrange(gene.start-1, gene.end):
                if not fwd_gene_dict.has_key(x):
                    fwd_gene_dict.update({x:gene})
                else:
                    # if two fwd gene overlap, assign this pos to latter one
                    fwd_gene_dict[x] = gene
        else:
            assert gene.strand == "-", "gene strand should be +/- !"
            rev_gene_starts_dict.update({gene.end-1:gene})
            # gff is 1-based, python is 0-based
            for x in xrange(gene.start-1, gene.end):
                if not rev_gene_dict.has_key(x):
                    rev_gene_dict.update({x:gene})
                else:
                    # if two rev gene overlap, assign this pos to formmer one
                    continue

    return fwd_gene_dict, rev_gene_dict, fwd_gene_starts_dict, rev_gene_starts_dict









######################## group tss ############################################

def group_tss(tss_arr, cov_arr, strand, before=5, after=5):
    """ This function used to group raw given strand tss in tss_arr, in range
        regionLen.

        default regionLen = 10, 5 before and 5 after
    """


    # count how many tss we have before grouping
    oldNum = Counter(tss_arr > 0)[True]

    # use a dict to store tss
    tss_dict = {} # data format {index:TSS_instance}
    # index array for sorted tss_arr from small to large value
    sort_index = np.argsort(tss_arr)
    total_length = len(tss_arr)

    # inverse to make high TSS number's index in beginning
    sort_index = sort_index[::-1]

    for index in sort_index:
        # maybe this pos is already merged to a neighboring highly expressed tss
        if tss_arr[index] == 0:
            continue

        else:
            # get leftmost and rightmost pos of this region, rightmost is not inclusive
            fore_index = index-before if index-before >0 else 0
            back_index = index+after+1 if index+after+1 < total_length else total_length

            # get sum of this region, assign it to tss_arr[index]
            tss_arr[index] = tss_arr[index]+\
                             np.sum(tss_arr[fore_index:index])+\
                             np.sum(tss_arr[index+1:back_index])

            # set before and after region of tss_arr[index] to 0
            tss_arr[fore_index:index] = 0
            tss_arr[index+1:back_index] = 0

            # get max cov in this region
            cov = max(cov_arr[fore_index:back_index])

            assert not tss_dict.has_key(index), "index aready in tss_dict !"

            # initialize TSS instance, and store into tss_dict
            tss_dict.update({index:TSS(strand, index, tss_arr[index], cov)})


    # count how many tss we have after grouping
    newNum = Counter(tss_arr > 0)[True]
    print "\nSuccessfully grouped %d raw tss into %d raw tss in region %d before and %d after\n"%(oldNum, newNum, before, after)

    return tss_arr, tss_dict







######################### calculate tss/cov ratio  #############################

def calculate_ratio(input_tss, input_cov):
    """ this function used to calculate tss/cov ratio for input tss and cov
    """
    ratio = []
    for t, c in izip(input_tss, input_cov):
        if c == 0:
            ratio.append(0)
        else:
            ratio.append(t/float(c))
    return np.array(ratio)







######################## find nearest gene ####################################

def getLeft(pos, Dict):
    """Find gene in Dict, at left direction of pos, return gene instance or None"""
    left = False
    while not left:
        pos -= 1
        if Dict.has_key(pos):
            left=True
            return Dict[pos]
        if pos < 0:
            left=True
            return None


def getRight(pos, Dict, genomeLen):
    """Find gene in Dict, at right direction of pos, return gene instance or None"""
    right = False
    while not right:
        pos += 1
        if Dict.has_key(pos):
            right = True
            return Dict[pos]
        if pos > genomeLen:
            right=True
            return None



def find_nearest_gene(tss, fwd_gene_dict, rev_gene_dict, fwd_gene_starts_dict,
                                           rev_gene_starts_dict, genome_length):
    """ take tss object and fwd/rev gene dict as inputs, will look for current
        gene, left gene and right gene for this tss.
    """
    pos, strand = tss.pos, tss.strand


    # fwd tss
    if strand == "+":
        # update currentGene
        currentGene = fwd_gene_dict.get(pos, None)
        tss.currentGene = currentGene

        # check current pos if a start or not
        currentStart = fwd_gene_starts_dict.get(pos, None)

        # update downstreamGene
        if currentGene and not currentStart:
            downstreamGene = getRight(currentGene.end, fwd_gene_dict, genome_length)
        else:
            downstreamGene = getRight(pos, fwd_gene_dict, genome_length)
        tss.downstreamGene = downstreamGene

        # update antiNearestGene
        antiCurrent = rev_gene_dict.get(pos, None)
        if antiCurrent:
            antiNearestGene = antiCurrent
        else:
            antiLeft = getLeft(pos, rev_gene_dict)
            antiRight = getRight(pos, rev_gene_dict, genome_length)
            if antiLeft and antiRight:
                if pos - (antiLeft.end-1) >= (antiRight.start-1) - pos:
                    antiNearestGene = antiRight
                else:
                    antiNearestGene = antiLeft
            else:
                if antiLeft:
                    antiNearestGene = antiLeft
                elif antiRight:
                    antiNearestGene = antiRight
                else:
                    antiNearestGene = None
        tss.antiNearestGene = antiNearestGene

    else:
        assert strand == "-", "strand should be +/-, not %s !"%strand
        # update currentGene
        currentGene = rev_gene_dict.get(pos, None)
        tss.currentGene = currentGene

        # check current pos if a start or not
        currentStart = rev_gene_starts_dict.get(pos, None)

        # update downstreamGene
        if currentGene and not currentStart:
            downstreamGene = getLeft(currentGene.start-2, rev_gene_dict)
        else:
            downstreamGene = getLeft(pos, rev_gene_dict)
        tss.downstreamGene = downstreamGene

        # update antiNearestGene
        antiCurrent = fwd_gene_dict.get(pos, None)
        if antiCurrent:
            antiNearestGene = antiCurrent
        else:
            antiLeft = getLeft(pos, fwd_gene_dict)
            antiRight = getRight(pos, fwd_gene_dict, genome_length)
            if antiLeft and antiRight:
                if pos - (antiLeft.end-1) >= (antiRight.start-1)-pos:
                    antiNearestGene = antiRight
                else:
                    antiNearestGene = antiLeft
            else:
                if antiLeft:
                    antiNearestGene = antiLeft
                elif antiRight:
                    antiNearestGene = antiRight
                else:
                    antiNearestGene = None
        tss.antiNearestGene = antiNearestGene








#########################  classify tss into catagories  #######################
"""
TSS classification as follows:

----------TSS----Gene==============>------------------- gTSS,
                                             ATG-200 <= pos <= ATG
----------Gene==========TSS========>------------------- iTSS,
                                             ATG < pos < END
----------Gene==========SST========>------------------- aTSS,
                                             ATG-50 <= pos <= END+50
----------GeneA==>---------TSS-----------GeneB==>------ nTSS,
                                             any other condition except above.

"""


def gTSS_check(tss):
    """ this function used to check gTSS
    """
    downstreamGene = tss.downstreamGene
    if downstreamGene:
        # tss locates at fwd
        if tss.strand == "+":
            downstream_dist = (downstreamGene.start-1)-tss.pos
        # tss locates at rev
        else:
            assert tss.strand == "-", "tss strand should be +/- !"
            downstream_dist = tss.pos - (downstreamGene.end-1)

        # update tss type, description and product
        # gTSS, upstream of 200 nt
        if downstream_dist <= 200:
            tss.type = "gTSS"
            tss.description = "gTSS, %d bp upstream of gene %s."%\
                            (downstream_dist, downstreamGene.Name)
            tss.product = downstreamGene.product


        # gTSS, more than 200nt, the merged coverage region around TSS overlaps with downstream gene,
        # and the overlaps more than 1/3 of the gene
        elif tss.merged_cov_region != None:
            if tss.strand == "+" and abs(tss.merged_cov_region[1]) > downstreamGene.start:
            #if tss.strand == "+" and abs(tss.merged_cov_region[1] - downstreamGene.start) > 0.333*abs(downstreamGene.start - downstreamGene.end):
                #print "One gTSS long than 200 nt found on fwd strand"
                tss.type = "gTSS"
                tss.description = "gTSS, %d bp upstream of gene %s."%\
                                (downstream_dist, downstreamGene.Name)
                tss.product = downstreamGene.product
            elif tss.strand == "-" and tss.merged_cov_region[0] < downstreamGene.end:
            #elif tss.strand == "-" and abs(tss.merged_cov_region[0] - downstreamGene.end) > 0.333*abs(downstreamGene.start-downstreamGene.end):
                #print "One gTSS long than 200 nt found on rev strand"
                tss.type = "gTSS"
                tss.description = "gTSS, %d bp upstream of gene %s."%\
                                (downstream_dist, downstreamGene.Name)
                tss.product = downstreamGene.product


def iTSS_check(tss):
    """ this function used to check iTSS
    """
    currentGene = tss.currentGene
    if currentGene:
        # tss locates at fwd
        if tss.strand == "+":
            inside_dist = tss.pos - (currentGene.start-1)
        # tss locates at rev
        else:
            assert tss.strand == "-", "tss strand should be +/- !"
            inside_dist = currentGene.end-1 - tss.pos
        # update tss type, description and product
        tss.type = "iTSS"
        tss.description = "iTSS, %d bp inside of gene %s."%\
                        (inside_dist, currentGene.Name)
        tss.product = currentGene.product


def aTSS_check(tss):
    """ this function used to check aTSS
    """
    antiNearestGene = tss.antiNearestGene

    # antiRight
    if antiNearestGene.start-1 > tss.pos:
        anti_dist = (antiNearestGene.start-1) - tss.pos
        if anti_dist <= 50:
            tss.type = "aTSS"
            tss.product = antiNearestGene.product
            if tss.strand == "+":
                tss.description = "aTSS, %d bp downstream of reverse gene %s"%\
                                            (anti_dist, antiNearestGene.Name)
            else:
                tss.description = "aTSS, %d bp upstream of forward gene %s"%\
                                            (anti_dist, antiNearestGene.Name)

    # antiLeft
    elif antiNearestGene.end-1 < tss.pos:
        anti_dist = tss.pos-(antiNearestGene.end-1)
        if anti_dist <= 50:
            tss.type = "aTSS"
            tss.product = antiNearestGene.product
            if tss.strand == "+":
                tss.description = "aTSS, %d bp upstream of reverse gene %s"%\
                                            (anti_dist, antiNearestGene.Name)
            else:
                tss.description = "aTSS, %d bp downstream of forward gene %s"%\
                                            (anti_dist, antiNearestGene.Name)

    # antiCurrent
    else:
        tss.type = "aTSS"
        tss.product = antiNearestGene.product
        if tss.strand == "+":
            anti_dist = (antiNearestGene.end-1) - tss.pos
            tss.description = "aTSS, %d bp downstream of reverse gene %s"\
                            %(anti_dist, antiNearestGene.Name)
        else:
            anti_dist = tss.pos - (antiNearestGene.start-1)
            tss.description = "aTSS, %d bp downstream of forward gene %s"\
                            %(anti_dist, antiNearestGene.Name)


def nTSS_check(tss):
    """ this function used to check nTSS
    """
    assert tss.currentGene == None, "currentGene exists, should not be nTSS !"

    antiNearestGene, downstreamGene = tss.antiNearestGene, tss.downstreamGene
    # change to 1 based coordinate
    tss.description = "nTSS, at position %d on %s strand"%(tss.pos+1, tss.strand)

    if downstreamGene:
        if tss.strand == "+":
            assert (downstreamGene.start-1) - tss.pos > 200, "distance between tss\
                and downstreamGene is less than 200, should not be nTSS !"
        else:
            assert tss.pos - (downstreamGene.end -1) > 200, "distance between tss\
                and downstreamGene is less than 200, should not be nTSS !"
        tss.description += ", before gene %s"%(downstreamGene.Name)

    if antiNearestGene:
        assert tss.pos - (antiNearestGene.end-1) > 50 or \
                (antiNearestGene.start-1) - tss.pos > 50, "distance between tss \
                and antiNearestGene is less than 50, should not be nTSS !"
    tss.description += ", antistrand nearest gene is %s"%(antiNearestGene.Name)


    tss.type = "nTSS"
    tss.product = "nTSS"




def classify_tss_type_and_update_description(tss):
    """ this function will be used to classify tss type and update description field
    """
    if not tss.type:
        gTSS_check(tss)
    if not tss.type:
        iTSS_check(tss)
    if not tss.type:
        aTSS_check(tss)
    if not tss.type:
        nTSS_check(tss)

    assert tss.type != None, "One kind of TSS should be assigned !"







############################# local enrichment  ###############################

def do_local_TSS_enrichment(tss, tssDict, length):
    """This function used to calculate the local tss enrichment score for a given tss and enrichment length,
       the more a TSS obvious, the surrounding tss are more less obvious, in a certain region, like 300nt,
       if we assume only one TSS for each genes in this region.

    """
    currCount = tss.init
    if tss.strand == "+": # for the fwd tss, left is smaller, right is larger
        # calculate the sum of left part(not include tss.init)
        currPos = tss.pos
        sumleft = 0
        while True:
            currPos -= 1
            if tssDict.has_key(currPos):
                sumleft += tssDict[currPos].init
            if currPos <= tss.pos -length: # when currPos is smaller than the left boundary
                break
        # calculate the sum of right part(not include tss.init)
        currPos = tss.pos
        sumright = 0
        while True:
            currPos += 1
            if tssDict.has_key(currPos):
                sumright += tssDict[currPos].init
            if currPos >= tss.pos + length: # when currPos is larger than the right boundary
                break
        #le = (sumleft+currCount)/(sumleft+sumright+currCount)  # local enrichment
        le = currCount/(sumleft+sumright+currCount)
    else: # for the rev tss, left is larger, and the left is smaller
        # calculate the sum of left part(not include tss.init)
        currPos = tss.pos
        sumleft = 0
        while True:
            currPos += 1
            if tssDict.has_key(currPos):
                sumleft += tssDict[currPos].init
            if currPos >= tss.pos + length: # when currPos is larger than the left boundary
                break
        # calculate the sum of right part(not include tss.init)
        currPos = tss.pos
        sumright = 0
        while True:
            currPos -= 1
            if tssDict.has_key(currPos):
                sumright += tssDict[currPos].init
            if currPos <= tss.pos - length: # when currPos is larger than the right boundary
                break
        #le = (sumleft+currCount)/(sumleft+sumright+currCount) # local enrichment
        le = currCount/(sumleft+sumright+currCount)

    # update localEnrich and donwnstreamUpstreamRatio info
    tss.locTssEnrich = le



def do_local_cov_enrichment(tss, coverage, length=100):
    """This function used to calculate the local coverage enrichment score for a given tss and enrichment length,
       the more one TSS likely, the slope around this TSS should be more steep, so we can calculate the coverage
       after this TSS, devide by coverage around this TSS, to get an coverage enrichment score
    """
    if tss.strand == "+":
        sum_upstream = sum(coverage[tss.pos-length:tss.pos])
        sum_downstream = sum(coverage[tss.pos:tss.pos+length]) # downstream including the TSS position
    else:
        assert tss.strand == "-", "TSS strand should be -, not %s"%tss.strand
        sum_upstream = sum(coverage[tss.pos+1:tss.pos+1+length])
        sum_downstream = sum(coverage[tss.pos+1-length:tss.pos+1]) # downstream including the TSS position

    local_cov_enrichment = float(sum_downstream)/(sum_upstream+sum_downstream)
    tss.locCovEnrich = local_cov_enrichment



class TSS():

    def __init__(self, strand, pos, init, cov):
        self.strand = strand
        self.pos = pos
        self.init = init
        self.cov = cov
        # here if no cov, then we assume this should not be a tss, set ratio to 0
        self.ratio = 0 if self.cov == 0 else self.init/float(self.cov) # tss/cov
        self.ID = None
        self.locTssEnrich = None       # (upstream+TSS)/(upstream+TSS+downstream), the more sharp in up and down, the more like to be TSS
        self.locCovEnrich = None       # downstream/(upstream+downstream), the more adjacent to 1, the more steep after this TSS, the more likely
        self.type = None
        self.currentGene = None      # which gene tss locates in current strand
        self.downstreamGene = None   # gene downstream of tss in current strand
        self.antiNearestGene = None  # nearest gene in anti strand
        self.description = None
        self.product = None
        self.rev_ratio = None            # the ratio of TSS/rev_cov, the rev_cov is dRNA_RV, and take maximun of clustered region
        self.rdm_ratio = None            # the ratio of TSS/rdm_cov, the rdm_cov is merged rdm_FV and rdm_RV, and take maximun of clustered region
        self.merged_cov_region = None    # a tuple represent the start and end pos of merged coverage region, will be used in gTSS_check
        self.revSum_ratio = 0  # TSS/revSum, revSum is the sum of reads count in dRNA-Seq RV in 500 nt

    def update_ID(self):
        if self.type != None:
            self.ID = self.type + self.strand + str(self.pos+1)

    def _replaceNone(self, item):
        if item == None:
            return "None"
        else:
            return str(item)

    def _get_Name(self, item):
        try:
            name = item.Name
        except Exception as e:
            name = "None"
        return name

    def _get_Locus_tag(self, item):
        try:
            locus_tag = item.locus_tag
        except Exception as e:
            locus_tag = "None"
        return locus_tag


    def __str__(self):
        """ seqID, genome, feature, start, end, product, strand, other, attributes(add color=255 0 0 to make tss display in red)
        """
        # guide sentence for rdm_ratio
        rdm_ratio = self.rdm_ratio if self.rdm_ratio else 1

        return self.ID+"\t"+str(self.init)+"\t"+self.type+"\t"+str(self.pos+1)+"\t"+\
              str(self.pos+1)+"\t"+ self._replaceNone(self.description)+"\t"+\
              self.strand+"\t"+self._replaceNone(self.product)+"\t"+\
              "ID=%s;Name=%s;LocalTssEnrichmentScore=%.2f;LocalCoverageEnrichmentScore=%.2f;TssCovRatio=%.2f;TssRevSumRatio=%.2f;TssRevCovRatio=%.2f;TssRdmCovRatio=%.2f;CurrentGene=%s;DownstreamGene=%s;NearestAntiStrandGene=%s;DownstreamGeneLocusTag=%s;color=255 0 0"%\
              (self.ID, self.ID, self.locTssEnrich, self.locCovEnrich, self.ratio, self.revSum_ratio, self.rev_ratio, rdm_ratio,
              self._get_Name(self.currentGene), self._get_Name(self.downstreamGene),
              self._get_Name(self.antiNearestGene), self._get_Locus_tag(self.downstreamGene)
              )

def do_calculate_rdm_ratio(tss, rdm_cov):
    """ this function will be used to calculate the rdm_ratio, namely tss/rdm_cov ratio,
        the rdm_cov is the maximum covrage from the coresponding positions of where TSS
        clustered from, namely -5 to +5

        here the rdm_cov, is the average of rdm_FV_cov and rdm_RV_cov, it's a average of
        both libraries, to make the coverage more robust
    """
    # first, take care of the tss clustering region, calculate tss.rdm_ratio
    cluster_region = rdm_cov[tss.pos-6:tss.pos+5]
    cov = max(cluster_region)
    tss.rdm_ratio = 1 if cov == 0 else tss.init/cov



def do_calculate_rev_ratio(tss, rev_cov):
    """ this function will be used to calcuate the rev_ratio, namely tss/rev_cov ratio,
        the rdm_cov is the maximun coverage from the corresponding positions of where TSS
        clustered from, namely -5 to +5
    """
    # first, take care of the tss clustering region, calculate tss.rev_ratio
    cluster_region = rev_cov[tss.pos-6:tss.pos+5]
    cov = max(cluster_region)
    tss.rev_ratio = 1 if cov == 0 else tss.init/cov




def do_update_merged_cov_region(tss, merged_cov):
    """ this function will be used to record the coverage > 1 region downstream
        of tss, and will get a start and end tuple, recording the coordinate of
        coverage region, which will be used in gTSS prediction
    """
    # fwd tss, record coverage values after tss
    if tss.strand == "+":
        # then find all the pos that more than 1, to find the end region
        for i, v in enumerate(merged_cov[tss.pos:]):
            if v > 1:
                continue
            else:
                tss.merged_cov_region =tuple([tss.pos, tss.pos+i])
                break
    # rev tss, record coverage values downstream tss
    else:
        assert tss.strand == "-", "tss strand should be -, not %s"%tss.strand
        # then find all the pos that more than 1, append to the rdm_region list
        rev_region = merged_cov[0:tss.pos]
        for i, v in enumerate(rev_region[::-1]):
            if v > 1:
                continue
            else:
                tss.merged_cov_region = tuple([tss.pos-i, tss.pos])
                break


def do_update_revSum(tss, dRNA_RV_tss, dist=500):
    """ this function will be used to mark TSS has a dRNA-Seq REV region or not,
        it will searching 500nt after TSS, then add up all the REV reads count
        in this region, then get a TSS/revSum ratio
    """
    # fwd tss, check after 500 nt, else check before 500 nt, to see if has coverage in dRNA_RV
    RV_tss_region = None
    if tss.strand == "+":
        RV_tss_region = dRNA_RV_tss[tss.pos:tss.pos+dist+1]
    else:
        assert tss.strand == "-", "tss strand should be +, not %s"%tss.strand
        RV_tss_region = dRNA_RV_tss[tss.pos-dist:tss.pos+1]

    revSum = sum(RV_tss_region)
    # give an abitrary large number
    if revSum == 0:
        tss.revSum_ratio = 1e10
    else:
        tss.revSum_ratio = tss.init/revSum


def do_update_revMax(tss, dRNA_RV_cov, dist=500, minimum=2):
    """ this function will be used to mark TSS has a dRNA-Seq REV region or not,
        it will searching for a no-zero coverage region after TSS before dist,
        then get revMax in this region, and also TSS/revMax ratio
    """

    # fwd tss, check after 500 nt, else check before 500 nt, to see if has coverage in dRNA_RV
    RV_cov_region = None
    if tss.strand == "+":
        RV_cov_region = dRNA_RV_cov[tss.pos:tss.pos+dist+1]
    else:
        assert tss.strand == "-", "tss strand should be +, not %s"%tss.strand
        RV_cov_region = dRNA_RV_cov[tss.pos-dist:tss.pos+1]
        RV_cov_region = RV_cov_region[::-1]

    revMax = 0
    for i, v in enumerate(RV_cov_region):
        if v == 0:
            # before the first non-zero region
            if revMax == 0:
                continue
            # after the first non-zero region
            else:
                break
        # in the first non-zero region
        else:
            if abs(v) >= revMax:
                revMax = abs(v)

    # finally, if after looping through, give the tss.revMax and tss.revMax_ratio
    if revMax >=2:
        tss.revMax = revMax
        tss.revMax_ratio = tss.init/float(revMax)
    else:
        tss.revMax = 0
        tss.revMax_ratio = 1




def merge_grp_coverage(grp_lists, average=True):
    """ this function will be used to merge the coverage from grp_lists,
        if average, then the merged coverage will be divided by the length of grp_lists
    """
    # read in first one
    fwd_cov, fwd_tss, rev_cov, rev_tss = read_grp(grp_lists[0])

    for grp_file in grp_lists[1:]:
        tmp_fwd_cov, tmp_fwd_tss, tmp_rev_cov, tmp_rev_tss = read_grp(grp_file)
        fwd_cov += tmp_fwd_cov
        rev_cov += tmp_rev_cov

    if average:
        fwd_cov = fwd_cov/len(grp_lists)
        rev_cov = rev_cov/len(grp_lists)

    return fwd_cov, rev_cov





def main():

    if len(sys.argv) != 6:
        print "\nUsage: ./tss_prediction.py [input_dRNA_FV_grp_file] [input_dRNA_RV_grp_file] [input_rdm_FV_grp_file] [input_rdm_RV_grp_file] [input_gff_file]\n"
        print "Compared to v3:\n\n1) here we include dRNA_RV_grp, to make sure than every FV tss, will have a rev coverage in 500bp"
        print "   and take the sum of read number in this region, to give a TSS/revSum ratio."
        print "\n2) we also included a rev_ratio, that use TSS/rev_cov, rev_cov is the maximun of clustered TSS region on dRNA-Seq RV libraries,"
        print "   to test if a TSS is buried in the previous TSS' RV coverage, only when it exceed 0.5 of the coverage, we will take it a TSS."
        print "\n3) we also change the rdm_FV_grp, to use merged rdm_FV and rdm_RV grp, to get rdm_ratio\n"
        print "\n4) moreover, the calculation of previous local enrichment score was changed to TSS/(upstreamTSS+TSS+downstreamTSS), to get rid of neighbored small TSS\n"
        print "\n5) and coverage enrichment score were also calculated, downstreamCov/upstreamCov, which is more like the slope, to define the sharpness of 5'End\n"
        print "\n6) The gTSS_check was also modified to use merged_cov_region, then the coverage region will be more accurate in low reads condition\n"
        print "\n7) The overlap of merged coverage region and downstream gene should be more than 1/3 of the downstream gene, for long distance gTSS\n"
        sys.exit(0)


    #---------------------------------------
    # input parameters
    input_dRNA_FV_grp = sys.argv[1]
    input_dRNA_RV_grp = sys.argv[2]
    input_rdm_FV_grp = sys.argv[3]
    input_rdm_RV_grp = sys.argv[4]
    input_gff = sys.argv[5]




    #--------------------------------------------------------------------------\
    # get basename of input dRNA_FV grp file
    basename = os.path.splitext(os.path.split(input_dRNA_FV_grp)[-1])[0]

    # read in grp file, "BASE \t fwd_coverage \t fwd_tss \t rev_coverage \t rev_tss"
    fwd_cov, fwd_tss, rev_cov, rev_tss = read_grp(input_dRNA_FV_grp)

    # get genome length
    genome_length = len(fwd_cov)


    #--------------------------------------------------------------------------\
    # group tss, and initialize tss instances
    fwd_tss, fwd_tss_dict = group_tss(fwd_tss, fwd_cov, "+", before=5, after=5)
    rev_tss, rev_tss_dict = group_tss(rev_tss, rev_cov, "-", before=5, after=5)


    #--------------------------------------------------------------------------\
    # update the revSum ratio between TSS and the sum of RV reads in 500nt

    # read in dRNA_RV grp file, "# BASE fwd_coverage fwd_tss rev_coverage rev_tss"
    RV_fwd_cov, RV_fwd_tss, RV_rev_cov, RV_rev_tss = read_grp(input_dRNA_RV_grp)

    for tss in fwd_tss_dict.itervalues():
        # calculate revSum_ratio
        do_update_revSum(tss, RV_fwd_tss)
        # calculate rev_ratio
        do_calculate_rev_ratio(tss, RV_fwd_cov)

    for tss in rev_tss_dict.itervalues():
        # calculate revSum_ratio
        do_update_revSum(tss, RV_rev_tss)
        # calculate rev_ratio
        do_calculate_rev_ratio(tss, RV_rev_cov)


    #--------------------------------------------------------------------------\
    # update rdm_cov region, and merged dRNA and rdm coverage
    # first, merge the rdm_FV and rdm_RV coverage
    rdm_fwd_cov, rdm_rev_cov = merge_grp_coverage([input_rdm_FV_grp, input_rdm_RV_grp],average=True)
    merged_fwd_cov, merged_rev_cov = merge_grp_coverage([input_dRNA_FV_grp, input_dRNA_RV_grp, input_rdm_FV_grp, input_rdm_RV_grp])

    for tss in fwd_tss_dict.itervalues():
        # calculate rdm_ratio, tss/rdm_cov ratio
        do_calculate_rdm_ratio(tss, rdm_fwd_cov)
        # update merged_cov_region
        do_update_merged_cov_region(tss, merged_fwd_cov)
        #print tss.pos, tss.init, tss.ratio, tss.strand, tss.rdm_ratio, tss.rdm_region
    for tss in rev_tss_dict.itervalues():
        # calculate rdm_ratio, tss/rdm_cov ratio
        do_calculate_rdm_ratio(tss, rdm_rev_cov)
        # update merged_cov_region
        do_update_merged_cov_region(tss, merged_rev_cov)
        #print tss.pos, tss.init, tss.ratio, tss.strand, tss.rdm_ratio, tss.rdm_region


    #--------------------------------------------------------------------------\
    # update local enrichment score
    for tss in fwd_tss_dict.itervalues():
        do_local_TSS_enrichment(tss, fwd_tss_dict, 100)
        do_local_cov_enrichment(tss, fwd_cov, length=10)
    for tss in rev_tss_dict.itervalues():
        do_local_TSS_enrichment(tss, rev_tss_dict, 100)
        do_local_cov_enrichment(tss, rev_cov, length=10)



    #--------------------------------------------------------------------------\
    # get gene_dict from gff file
    # data format {pos: gene}
    fwd_gene_dict, rev_gene_dict, fwd_gene_starts_dict, rev_gene_starts_dict = \
                                             get_gene_dict_from_gff(input_gff)



    #--------------------------------------------------------------------------\
    # find nearest gene for each tss
    for pos, tss in fwd_tss_dict.iteritems():
        find_nearest_gene(tss, fwd_gene_dict, rev_gene_dict, fwd_gene_starts_dict,
                          rev_gene_starts_dict, genome_length)

    for pos, tss in rev_tss_dict.iteritems():
        find_nearest_gene(tss, fwd_gene_dict, rev_gene_dict, fwd_gene_starts_dict,
                          rev_gene_starts_dict, genome_length)



    #--------------------------------------------------------------------------\
    # classify tss type, and update tss description
    for tss in fwd_tss_dict.itervalues():
        classify_tss_type_and_update_description(tss)
        tss.update_ID()
    for tss in rev_tss_dict.itervalues():
        classify_tss_type_and_update_description(tss)
        tss.update_ID()



    #--------------------------------------------------------------------------\
    # set cutoff and writeout total tss
    total_tss = []
    bad_tss = []

    for tss in fwd_tss_dict.itervalues():
        # basic four parameters
        if tss.ratio > 0.5 and tss.locTssEnrich > 0.3 and tss.locCovEnrich > 0.5 and tss.revSum_ratio <= 30:# and tss.rev_ratio > 0.75 #and tss.revMax_ratio > 0.1 :
            # only tss that has lower than 100 dRNA_FV counts, are subjected to the fifth parameters
            if tss.init >= 100:
                total_tss.append(tss)
            else:
                # tss with less than 100 dRNA_FV counts, that pass rdm_ratio test to be good tss
                if tss.rdm_ratio > 0.85:
                    total_tss.append(tss)
                else:
                    bad_tss.append(tss)
        # tss that failed in the first five parameters
        else:
            bad_tss.append(tss)

    for tss in rev_tss_dict.itervalues():
        # basic four parameters
        if tss.ratio > 0.5 and tss.locTssEnrich > 0.3 and tss.locCovEnrich > 0.5 and tss.revSum_ratio <= 30:# and tss.rev_ratio > 0.75 #and tss.revMax_ratio > 0.1 :
            # tss with lower than 100 dRNA_FV counts, were subjected to the additional fifth parameters
            if tss.init >= 100:
                total_tss.append(tss)
            else:
                # tss with less than 100 dRNA_FV coutns
                if tss.rdm_ratio > 0.85:
                    total_tss.append(tss)
                else:
                    bad_tss.append(tss)
        # tss that failed in the first five parameters
        else:
            bad_tss.append(tss)

    #--------------------------------------------------------------------------\
    # get the cutoff from iTSS of bad tss

    # get all bad iTSS counts
    bad_iTSS_counts = [tss.init for tss in bad_tss if tss.type=="iTSS"]

    # use 95% percentile of bad iTSS, as cutoff for tss prediction
    #percentile = 80
    #cutoff = np.percentile(np.array(bad_iTSS_counts), percentile)
    cutoff = 4
    #print "reads count for %d %% percentile of bad iTSS is %.2f\n"%(percentile, cutoff)



    #--------------------------------------------------------------------------\
    # write out all good and bad tss
    good_tss_file = basename+"_cutoff_%.2f"%cutoff+"_good_v4.tss"
    bad_tss_file = basename+"_cutoff_%.2f"%cutoff+"_bad_v4.tss"
    with open(good_tss_file, "w") as good, open(bad_tss_file, "w") as bad:
        # write out good tss
        good.write("#ID\tExpression\tType\tStart\tEnd\tDescription\tStrand\tProduct\t\
                  Attributes\n")
        # sort by position
        total_tss = sorted(total_tss, key = lambda x: x.pos)
        for tss in total_tss:
            if tss.init >= cutoff:
                good.write(str(tss)+"\n")
            else:
                bad_tss.append(tss)
        # write out bad tss
        bad.write("#ID\tExpression\tType\tStart\tEnd\tDescription\tStrand\tProduct\t\
                  Attributes\n")
        # sort bad tss based on init counts
        bad_tss = sorted(bad_tss, key=lambda x: x.init, reverse=True)
        for tss in bad_tss:
            bad.write(str(tss)+"\n")



    tss_dict = {}
    with open(good_tss_file, "r") as ih:
        for line in ih:
            if line.startswith("#"):
                continue
            line = line.strip().split("\t")
            tssType = line[2]
            if tss_dict.has_key(tssType):
                tss_dict[tssType] += 1
            else:
                tss_dict.update({tssType:1})

    od = collections.OrderedDict(sorted(tss_dict.items()))
    print "\n\n\n---------------------------------\n"
    print " ".join(sys.argv)
    for k, v in od.iteritems():
        print k+"\t"+str(v)
    print "--------------------------------\n\n\n"



if __name__ == "__main__":
    main()
