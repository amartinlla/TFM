#!/usr/bin/env python
# -*- coding: utf-8 -*-

import collections
import re
import argparse
import sys
import codecs
import gzip
import itertools

import psycopg2, psycopg2.extras


try:
    import cparse
except ImportError:
    cparse = None

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict


from model import _Call, _Record, make_calldata_tuple
from model import _Substitution, _Breakend, _SingleBreakend, _SV

from config import config

# Metadata parsers/constants
RESERVED_INFO = {
    'AA': 'String', 'AC': 'Integer', 'AF': 'Float', 'AN': 'Integer',
    'BQ': 'Float', 'CIGAR': 'String', 'DB': 'Flag', 'DP': 'Integer',
    'END': 'Integer', 'H2': 'Flag', 'H3': 'Flag', 'MQ': 'Float',
    'MQ0': 'Integer', 'NS': 'Integer', 'SB': 'String', 'SOMATIC': 'Flag',
    'VALIDATED': 'Flag', '1000G': 'Flag',

    # Keys used for structural variants
    'IMPRECISE': 'Flag', 'NOVEL': 'Flag', 'SVTYPE': 'String',
    'SVLEN': 'Integer', 'CIPOS': 'Integer', 'CIEND': 'Integer',
    'HOMLEN': 'Integer', 'HOMSEQ': 'String', 'BKPTID': 'String',
    'MEINFO': 'String', 'METRANS': 'String', 'DGVID': 'String',
    'DBVARID': 'String', 'DBRIPID': 'String', 'MATEID': 'String',
    'PARID': 'String', 'EVENT': 'String', 'CILEN': 'Integer',
    'DPADJ': 'Integer', 'CN': 'Integer', 'CNADJ': 'Integer',
    'CICN': 'Integer', 'CICNADJ': 'Integer'
}

RESERVED_FORMAT = {
    'GT': 'String', 'DP': 'Integer', 'FT': 'String', 'GL': 'Float',
    'GLE': 'String', 'PL': 'Integer', 'GP': 'Float', 'GQ': 'Integer',
    'HQ': 'Integer', 'PS': 'Integer', 'PQ': 'Integer', 'EC': 'Integer',
    'MQ': 'Integer',

    # Keys used for structural variants
    'CN': 'Integer', 'CNQ': 'Float', 'CNL': 'Float', 'NQ': 'Integer',
    'HAP': 'Integer', 'AHAP': 'Integer'
}

# Conversion between value in file and Python value
field_counts = {
    '.': None,  # Unknown number of values
    'A': -1,  # Equal to the number of alternate alleles in a given record
    'G': -2,  # Equal to the number of genotypes in a given record
    'R': -3,  # Equal to the number of alleles including reference in a given record
}

# METADATA

# Types of tags in the header

Info_tag = collections.namedtuple('Info', ['id', 'num', 'type', 'desc', 'source', 'version'])
Filter_tag = collections.namedtuple('Filter', ['id', 'desc'])
Alt_tag = collections.namedtuple('Alt', ['id', 'desc'])
Format_tag = collections.namedtuple('Format', ['id', 'num', 'type', 'desc'])
SampleInfo_tag = collections.namedtuple('SampleInfo', ['samples', 'gt_bases', 'gt_types', 'gt_phases'])
Contig_tag = collections.namedtuple('Contig', ['id', 'length'])
GATK_tag = collections.namedtuple('GATKCommandLine', ['id','version','date'])

# Spec is a bit weak on which metadata lines are singular, like fileformat
# and which can have repeats, like contig
SINGULAR_METADATA = ['fileformat', 'fileDate', 'reference', 'source', 'phasing']

info_pattern = re.compile(r'''\#\#INFO=<
        ID=(?P<id>[^,]+),\s*
        Number=(?P<number>-?\d+|\.|[AGR])?,\s*
        Type=(?P<type>Integer|Float|Flag|Character|String),\s*
        Description="(?P<desc>[^"]*)"
        (?:,\s*Source="(?P<source>[^"]*)")?
        (?:,\s*Version="?(?P<version>[^"]*)"?)?
        >''', re.VERBOSE)
filter_pattern = re.compile(r'''\#\#FILTER=<
    ID=(?P<id>[^,]+),\s*
    Description="(?P<desc>[^"]*)"
    >''', re.VERBOSE)
alt_pattern = re.compile(r'''\#\#ALT=<
    ID=(?P<id>[^,]+),\s*
    Description="(?P<desc>[^"]*)"
    >''', re.VERBOSE)
format_pattern = re.compile(r'''\#\#FORMAT=<
    ID=(?P<id>.+),\s*
    Number=(?P<number>-?\d+|\.|[AGR]),\s*
    Type=(?P<type>.+),\s*
    Description="(?P<desc>.*)"
    >''', re.VERBOSE)
contig_pattern = re.compile(r'''\#\#contig=<
    ID=(?P<id>[^>,]+)
    (,.*length=(?P<length>-?\d+))?
    .*
    >''', re.VERBOSE)
meta_pattern = re.compile(r'''##(?P<key>.+?)=(?P<val>.+)''')

# -------------------------------------- FUNCTIONS -----------------------------------------------------


def vcf_field_count(num_str):
    """Cast vcf header numbers to integer or None"""
    if num_str is None:
        return None
    elif num_str not in field_counts:
        # Fixed, specified number
        return int(num_str)
    else:
        return field_counts[num_str]

def read_info(info_string):
    '''Read a meta-information INFO line.'''
    match = info_pattern.match(info_string)
    if not match:
        raise SyntaxError(
            "One of the INFO lines is malformed: %s" % info_string)

    num = vcf_field_count(match.group('number'))

    info = Info_tag(match.group('id'), num,
                    match.group('type'), match.group('desc'),
                    match.group('source'), match.group('version'))

    return (match.group('id'), info)

def read_filter(filter_string):
    '''Read a meta-information FILTER line.'''
    match = filter_pattern.match(filter_string)
    if not match:
        raise SyntaxError(
            "One of the FILTER lines is malformed: %s" % filter_string)

    filt = Filter_tag(match.group('id'), match.group('desc'))

    return (match.group('id'), filt)

def read_alt(alt_string):
    '''Read a meta-information ALTline.'''
    match = alt_pattern.match(alt_string)
    if not match:
        raise SyntaxError(
            "One of the ALT lines is malformed: %s" % alt_string)

    alt = Alt_tag(match.group('id'), match.group('desc'))

    return (match.group('id'), alt)

def read_format(format_string):
    '''Read a meta-information FORMAT line.'''
    match = format_pattern.match(format_string)
    if not match:
        raise SyntaxError(
            "One of the FORMAT lines is malformed: %s" % format_string)

    num = vcf_field_count(match.group('number'))

    form = Format_tag(match.group('id'), num,
                    match.group('type'), match.group('desc'))

    return (match.group('id'), form)

def read_contig(contig_string):
    '''Read a meta-contigrmation INFO line.'''
    match = contig_pattern.match(contig_string)
    if not match:
        raise SyntaxError(
            "One of the contig lines is malformed: %s" % contig_string)
    length = vcf_field_count(match.group('length'))
    contig = Contig_tag(match.group('id'), length)
    return (match.group('id'), contig)

def read_meta_hash(meta_string):
    # assert re.match("##.+=<", meta_string)
    items = meta_string.split('=', 1)
    # Removing initial hash marks
    key = items[0].lstrip('#')
    # N.B., items can have quoted values, so cannot just split on comma
    val = OrderedDict()
    state = 0
    k = ''
    v = ''
    for c in items[1].strip('[<>]'):

        if state == 0:  # reading item key
            if c == '=':
                state = 1  # end of key, start reading value
            else:
                k += c  # extend key
        elif state == 1:  # reading item value
            if v == '' and c == '"':
                v += c  # include quote mark in value
                state = 2  # start reading quoted value
            elif c == ',':
                val[k] = v  # store parsed item
                state = 0  # read next key
                k = ''
                v = ''
            else:
                v += c
        elif state == 2:  # reading quoted item value
            if c == '"':
                v += c  # include quote mark in value
                state = 1  # end quoting
            else:
                v += c
    if k != '':
        val[k] = v
    #print key, val
    return key, val

def read_meta(meta_string):
    if re.match("##.+=<", meta_string):
        return read_meta_hash(meta_string)
    match = meta_pattern.match(meta_string)
    if not match:
        # Spec only allows key=value, but we try to be liberal and
        # interpret anything else as key=none (and all values are parsed
        # as strings).
        return meta_string.lstrip('#'), 'none'
    return match.group('key'), match.group('val')

def read_vcf(filename,f_out):

    # initialization

    #filename=None
    compressed=None
    strict_whitespace=False
    encoding='ascii'
    separator = '\t'
    header_lines = []
    formats = {}
    csq = False

    if compressed is None:
            compressed = filename.endswith('.gz')
        
    vcf_file = open(filename, 'rb' if compressed else 'rt')
       
    if compressed:
        vcf_file = gzip.GzipFile(fileobj=vcf_file)
        if sys.version > '3':
            vcf_file = codecs.getreader(encoding)(vcf_file)

        if strict_whitespace:
            separator = '\t'
        else:
            separator = '\t| +'

    row_pattern = re.compile(separator)
    #alt_pattern = re.compile('[\[\]]')

    for attr in ('metadata', 'infos', 'filters', 'alts', 'contigs', 'formats'):      
        dic = dict(zip(attr,OrderedDict()))


    reader = (line.strip() for line in vcf_file if line.strip())
    #print 'LEYENDO...'

    line = next(reader)
    
    #print line
    
    while line.startswith('##'):
        header_lines.append(line)
   

        if line.startswith('##INFO'):
            key, val = read_info(line)         
            dic[key] = val

        elif line.startswith('##FILTER'):
            key, val = read_filter(line)
            dic[key] = val

        elif line.startswith('##ALT'):
            key, val = read_alt(line)
            dic[key] = val

        elif line.startswith('##FORMAT'):
            key, val = read_format(line)
            dic[key] = val
            formats[key] = val


        elif line.startswith('##contig'):
            key, val = read_contig(line)
            dic[key] = val

        else:
            key, val = read_meta(line)
            if key in SINGULAR_METADATA:
                dic[key] = val
            else:
                if key not in dic:
                    dic[key] = []
                dic[key].append(val)

        line = next(reader)


    # Here comes the fields of the data

    # Reading CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT SAMPLE1 SAMPLE2 ...
    fields = row_pattern.split(line[1:])
  
    # We keep only the fields until the samples
    #column_headers = fields[:9]
   
    # Names of the samples
    samples = fields[9:]
    sample_indexes = dict([(x,i) for (i,x) in enumerate(samples)])


    if ('CSQ' in dic.keys()):
        names = str(dic['CSQ'][3]).split('|') # Index 3 corresponds to 'desc'       
        csq = True
  
    # Now begins the data
    num_linea = 1
    for line in vcf_file:

        record = read_data(line,row_pattern, dic,samples,sample_indexes,formats)
        print 'linea', num_linea
        num_linea += 1
      
        # Insert data in BBDD
        #var_id = insert_variant(record['chrom'],record['pos'],record['ID'],record['ref'],record['alt'],record['qual'],record['filt']) 

        f_out.write('INSERT INTO VARIANT(var_id,chrom,pos,id,ref,alt,qual,filter) VALUES (')
        f_out.write(str(num_linea)) 
        f_out.write(',')
        f_out.write(str(record['chrom'])) 
        f_out.write(',')
        f_out.write(str(record['pos'])) 
        f_out.write(",'")
        f_out.write(str(record['ID'])) 
        f_out.write("','")
        f_out.write(record['ref']) 
        f_out.write("','")
        f_out.write(record['alt']) 
        f_out.write("',")
        f_out.write(str(record['qual'])) 
        f_out.write(",'")
        f_out.write(str(record['filt'])) 
        f_out.write("');\n")


        
        if (csq):    
        
            lon = len(str(record['info']['CSQ']).split(","))
            preds = {}
            mafs = {}
            #for value in values_raw: # for each group of values with CSQ tag
            for ind in range(lon): # for each group of values with CSQ tag

                f_out.write('INSERT INTO TRANSCRIPT(fk_var,allele,consequence,impact,symbol,gene,feature_type,feature,biotype) VALUES(')

                f_out.write(str(num_linea))
                f_out.write(",'")
                for name, val in itertools.izip(names, str(record['info']['CSQ'][ind]).split("|")):  # for each value  
              
                    #print 'insertando...' ,name , val 
                
                    # TRANSCRIPTS
                    
                    if name == "Consequence annotations from Ensembl VEP. Format: Allele":                    
                        f_out.write(val)
                        f_out.write("','")            
                    elif name == "Consequence":
                        f_out.write(val)
                        f_out.write("','")
                    elif name == "IMPACT":
                        f_out.write(val)
                        f_out.write("','")                               
                    elif name == "SYMBOL":
                        f_out.write(val)
                        f_out.write("','")                   
                    elif name == "Gene":
                        f_out.write(val)
                        f_out.write("','")                     
                    elif name == "Feature_type":
                        f_out.write(val)
                        f_out.write("','")           
                    elif name == "Feature":
                        f_out.write(val)
                        f_out.write("','")
                     

                    if name == "BIOTYPE": 
                        f_out.write(val)
                        f_out.write("');\n")
                      
                    # PREDICTORS          
                    

                    if name == "SIFT":   
                        preds['sift'] = "'" + val +  "'"                 
                        #f_out.write(val)
                        #f_out.write(',')
                    elif name == "PolyPhen":
                        preds['poly'] = "'" + val +  "'"  
                        #f_out.write(val)
                        #f_out.write(',')
                    elif name == "LoFtool":
                        preds['lof'] = "'" + val +  "'"  
                        #f_out.write(val)
                        #f_out.write(',')
                    elif name == "CADD_PHRED":
                        preds['phred'] = "'" + val +  "'"  
                        #f_out.write(val)
                        #f_out.write(',')
                    elif name == "CADD_RAW":
                        preds['raw'] = "'" + val +  "'"  
                        #f_out.write(val)
                        #f_out.write(',')

                    # POP_MAFS
                    
                    if name == "gnomAD_AF":   
                        mafs['af'] = "'" + val + "'"                
                        #f_out.write(val)
                        #f_out.write(',')
                    elif name == "gnomAD_AFR_AF":
                        mafs['afr'] =  "'" + val +  "'"
                        #f_out.write(val)
                        #f_out.write(',')
                    elif name == "gnomAD_AMR_AF":
                        mafs['amr'] = "'" + val +  "'"  
                        #f_out.write(val)
                        #f_out.write(',')
                    elif name == "gnomAD_ASJ_AF":
                        mafs['asj'] = "'" + val +  "'"  
                        #f_out.write(val)
                        #f_out.write(',')
                    elif name == "gnomAD_EAS_AF":
                        mafs['eas'] = "'" + val +  "'"
                    elif name == "gnomAD_FIN_AF":
                        mafs['fin'] = "'" + val +  "'"
                    elif name == "gnomAD_NFE_AF":
                        mafs['nfe'] = "'" + val +  "'"
                    elif name == "gnomAD_OTH_AF":
                        mafs['oth'] = "'" + val +  "'"
                    elif name == "gnomAD_SAS_AF":
                        mafs['sas'] = "'" + val +  "'"


                f_out.write('INSERT INTO PREDICTORS(fk_trans,sift,polyphen,loftool,cadd_phred,cadd_raw) VALUES (')
                f_out.write(str(ind+1))
                f_out.write(',')
                f_out.write(preds['sift'])
                f_out.write(',')
                f_out.write(preds['poly'])  
                f_out.write(',')
                f_out.write(preds['lof'])  
                f_out.write(',')
                f_out.write(preds['phred'])  
                f_out.write(',')
                f_out.write(preds['raw'])  
                f_out.write(');\n')
           

                f_out.write('INSERT INTO POP_MAFS(fk_trans,af,afr_af,amr_af,asj_af,eas_af,fin_af,nfe_af,oth_af,sas_af) VALUES(')
                f_out.write(str(ind+1))
                f_out.write(',')
                f_out.write(mafs['af'])
                f_out.write(',')
                f_out.write(mafs['afr'])  
                f_out.write(',')
                f_out.write(mafs['amr'])  
                f_out.write(',')
                f_out.write(mafs['asj'])  
                f_out.write(',')
                f_out.write(mafs['eas'])  
                f_out.write(',')
                f_out.write(mafs['fin'])  
                f_out.write(',')
                f_out.write(mafs['nfe'])  
                f_out.write(',')
                f_out.write(mafs['oth'])  
                f_out.write(',')
                f_out.write(mafs['sas'])  
                f_out.write(');\n')


            # VALUES(%s,%s,%s,%s,%s)                  
                    #f_out.write('INSERT INTO info(fk_var,name,value,groups) VALUES(',var_id,',',name,',',val,',''CSQ'');\n')
                    
        sample_id = 0   
        for i in sample_indexes: # for each sample

            if str(record['samples'][i]['GT']).find("|") >= 0: #if the genotype is phased for this sample            
                genotypes = str(record['samples'][i]['GT']).split("|") # gets the value of each genotype in a list
            else:
                genotypes = str(record['samples'][i]['GT']).split("/")

            #print 'INSERTANDO...'    
            #sample_id = insert_sample(var_id,i,genotypes[0],genotypes[1])
            sample_id = sample_id + 1
            f_out.write('INSERT INTO SAMPLE(fk_var,name,ale1,ale2) VALUES(')
            f_out.write(str(num_linea))
            f_out.write(",'")
            f_out.write(str(i))
            f_out.write("','")  
            f_out.write(genotypes[0])
            f_out.write("','")  
            f_out.write(genotypes[1])
            f_out.write("');\n")
            

            for field in record['samples'][i]: # for each sample
                if field != 'GT':   
                    #print 'insert GT'                                   
                    #insert_gt_fields(sample_id,field,record['samples'][i][field])
                    f_out.write('INSERT INTO GT_FIELDS (fk_sample,name,value) VALUES(')
                    f_out.write(str(sample_id))
                    f_out.write(",'")
                    f_out.write(field)
                    f_out.write("','")
                    value = record['samples'][i][field]
                    f_out.write(str(value))
                    f_out.write("');\n")
            

    return record



def parse_filter(filt_str):
        '''Parse the FILTER field of a VCF entry into a Python list
        NOTE: this method has a cython equivalent and care must be taken
        to keep the two methods equivalent
        '''
        if filt_str == '.':
            return None
        elif filt_str == 'PASS':
            return []
        else:
            return filt_str.split(';')

def parse_info(infos,info_str):
    '''Parse the INFO field of a VCF entry into a dictionary of Python
    types.
    '''

    if info_str == '.':
        return {}

    entries = info_str.split(';')
    retdict = {}

    for entry in entries:
        entry = entry.split('=', 1)
        ID = entry[0]
        try:
            entry_type = infos[ID].type
        except KeyError:
            try:
                entry_type = RESERVED_INFO[ID]
            except KeyError:
                if entry[1:]:
                    entry_type = 'String'
                else:
                    entry_type = 'Flag'

        if entry_type == 'Integer':
            vals = entry[1].split(',')
            try:
                val = map(int, vals)
            # Allow specified integers to be flexibly parsed as floats.
            # Handles cases with incorrectly specified header types.
            except ValueError:
                val = map(float, vals)
        elif entry_type == 'Float':
            vals = entry[1].split(',')
            val = map(float, vals)
        elif entry_type == 'Flag':
            val = True
        elif entry_type in ('String', 'Character'):
            try:
                vals = entry[1].split(',') # commas are reserved characters indicating multiple values
                val = map(str, vals)
            except IndexError:
                entry_type = 'Flag'
                val = True

        try:
            if infos[ID].num == 1 and entry_type not in ( 'Flag', ):
                val = val[0]
        except KeyError:
            pass

        retdict[ID] = val

    return retdict


def parse_sample_format(samp_fmt,formats):
    """ Parse the format of the calls in this _Record """
    samp_fmt = make_calldata_tuple(samp_fmt.split(':'))

    for fmt in samp_fmt._fields:
        try:
            entry_type = formats[fmt].type
            entry_num = formats[fmt].num
        except KeyError:
            entry_num = None
            try:
                entry_type = RESERVED_FORMAT[fmt]
            except KeyError:
                entry_type = 'String'
        samp_fmt._types.append(entry_type)
        samp_fmt._nums.append(entry_num)
    return samp_fmt


def parse_samples(samples,sample_indexes, samp_fmt, record,format_cache,formats):
    '''Parse a sample entry according to the format specified in the FORMAT
    column.
    NOTE: this method has a cython equivalent and care must be taken
    to keep the two methods equivalent
    '''


    # check whether we already know how to parse this format
    
    if samp_fmt not in format_cache:
        format_cache[samp_fmt] = parse_sample_format(samp_fmt,formats)
    samp_fmt = format_cache[samp_fmt]

    if cparse:
        return cparse.parse_samples(samples,sample_indexes, samp_fmt, samp_fmt._types, samp_fmt._nums, record)

    _map = map
    nfields = len(samp_fmt._fields)
    smpls = {}

    # for each sample
    for name, sample in itertools.izip(sample_indexes, samples):

     
        # parse the data for this sample
        sampdat = [None] * nfields
 

        for i, vals in enumerate(sample.split(':')):

            # short circuit the most common
            if samp_fmt._fields[i] == 'GT':
                sampdat[i] = vals
                continue
            # genotype filters are a special case
            elif samp_fmt._fields[i] == 'FT':
                sampdat[i] = parse_filter(vals)
                continue
            elif not vals or vals == ".":
                sampdat[i] = None
                continue

            entry_num = samp_fmt._nums[i]
            entry_type = samp_fmt._types[i]

            # we don't need to split single entries
            if entry_num == 1:
                if entry_type == 'Integer':
                    try:
                        sampdat[i] = int(vals)
                    except ValueError:
                        sampdat[i] = float(vals)
                elif entry_type == 'Float' or entry_type == 'Numeric':
                    sampdat[i] = float(vals)
                else:
                    sampdat[i] = vals
                continue

            vals = vals.split(',')
            if entry_type == 'Integer':
                try:
                    sampdat[i] = _map(int, vals)
                except ValueError:
                    sampdat[i] = _map(float, vals)
            elif entry_type == 'Float' or entry_type == 'Numeric':
                sampdat[i] = _map(float, vals)
            else:
                sampdat[i] = vals
     
        gt_field_names = record['fmt']
        gt_field_values = sampdat
        
        j = 0
        dic = {}

        for key in gt_field_names.split(":"): 
            dic[key] = gt_field_values[j]
            j += 1

        smpls[name] = dic
        
    return smpls


def read_data(line,row_pattern, dic, samples,sample_indexes,formats):
    '''Return the next record in the file.'''
   
    prepend_chr = False
    record = {}
    format_cache = {}
   
    row = row_pattern.split(line.rstrip())

    chrom = row[0]

    if prepend_chr:
        chrom = 'chr' + chrom
    pos = int(row[1])

    if row[2] != '.':
        ID = row[2]
    else:
        ID = None

    ref = row[3]
    
    alt = row[4]

    #alt = map(parse_alt, row[4].split(','))

    try:
        qual = int(row[5])
    except ValueError:
        try:
            qual = float(row[5])
        except ValueError:
            qual = None

    
    filt = parse_filter(row[6])
    
    info = parse_info(dic,row[7])

    try:
        fmt = row[8]
    except IndexError:
        fmt = None
    else:
        if fmt == '.':
            fmt = None

    
    #print (chrom, pos, ID, ref, alt, qual, filt, info, fmt, sample_indexes)

    record['chrom'] = chrom
    record['pos'] = pos
    record['ID'] = ID
    record['ref'] = ref
    record['alt'] = alt
    record['qual'] = qual
    record['filt'] = filt
    record['info'] = info
    record['fmt'] = fmt
    record['sample_indexes'] = sample_indexes


    #record = [chrom, pos, ID, ref, alt, qual, filt, info, fmt, sample_indexes]


    if fmt is not None:
        record['samples'] = parse_samples(row[9:], sample_indexes, fmt, record,format_cache,formats)
      
      
    return record




def create_output_file():

    f = open ('out.sql','w')
    return f

def write_output_file(f,text):
    f.write(text)

def close_output_file(f):
    f.close()

def drop_tables(f):
    f.write("DROP TABLE IF EXISTS VARIANT CASCADE;\n")  
    f.write("DROP TABLE IF EXISTS TRANSCRIPT CASCADE;\n")  
    f.write("DROP TABLE IF EXISTS PREDICTORS CASCADE;\n")  
    f.write("DROP TABLE IF EXISTS POP_MAFS CASCADE;\n")  
    f.write("DROP TABLE IF EXISTS GT_FIELDS CASCADE;\n") 
    f.write("DROP TABLE IF EXISTS SAMPLE CASCADE;\n")                
    

def create_tables(f):
    f.write("""CREATE TABLE VARIANT
                (
                var_id integer NOT NULL,
                chrom integer NOT NULL,
                pos integer NOT NULL,
                id varchar(255),
                ref varchar(100),
                alt varchar(100),
                qual decimal(20,2),
                filter varchar(255),
                PRIMARY KEY (var_id)
                );\n""")  
    f.write("""CREATE TABLE TRANSCRIPT
                (
                trans_id SERIAL,
                fk_var integer NOT NULL,
				allele varchar(100),
				consequence varchar(100),
				impact varchar(100),
				symbol varchar(100),
				gene varchar(100),
                feature_type varchar(100),
				feature varchar(100),
				biotype varchar(100),				
                PRIMARY KEY (trans_id),
                FOREIGN KEY (fk_var) REFERENCES VARIANT (var_id) ON UPDATE CASCADE ON DELETE CASCADE
                );\n""")  
    f.write("""CREATE TABLE PREDICTORS
                (
				pred_id SERIAL,
                fk_trans integer NOT NULL,
				sift varchar(100),
				polyphen varchar(100),
				loftool varchar(100),
				cadd_phred varchar(100),
				cadd_raw varchar(100),				
                PRIMARY KEY (pred_id),
                FOREIGN KEY (fk_trans) REFERENCES TRANSCRIPT (trans_id) ON UPDATE CASCADE ON DELETE CASCADE
                );\n""")  
    f.write("""CREATE TABLE POP_MAFS
                (
				mafs_id SERIAL,
                fk_trans integer NOT NULL,
				af varchar(100),
				afr_af varchar(100),
				amr_af varchar(100),
				asj_af varchar(100),
				eas_af varchar(100),
				fin_af varchar(100),
				nfe_af varchar(100),
				oth_af varchar(100),
				sas_af varchar(100),
                PRIMARY KEY (mafs_id),
                FOREIGN KEY (fk_trans) REFERENCES TRANSCRIPT (trans_id) ON UPDATE CASCADE ON DELETE CASCADE
                );\n""")  
    f.write("""CREATE TABLE SAMPLE
                (
                sample_id SERIAL,
                fk_var integer NOT NULL,
                name varchar(255) NOT NULL,
                ale1 varchar(255),
                ale2 varchar(20),
                PRIMARY KEY (sample_id),
                FOREIGN KEY (fk_var) REFERENCES VARIANT (var_id) ON UPDATE CASCADE ON DELETE CASCADE
                );\n""")                
    f.write("""CREATE TABLE GT_FIELDS
                (
                gt_fields_id SERIAL,
                fk_sample integer NOT NULL,
                name varchar(20) NOT NULL,
                value varchar(255),
                PRIMARY KEY (gt_fields_id),
                FOREIGN KEY (fk_sample) REFERENCES SAMPLE (sample_id) ON UPDATE CASCADE ON DELETE CASCADE
                );\n""") 

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="Show debugger information", action="store_true")
    parser.add_argument("-f", "--file", help="Name of the file to process")
    args = parser.parse_args()
    
    # Aquí procesamos lo que se tiene que hacer con cada argumento
    if args.verbose:
        print "debbug mode on!"
    
    if args.file:
        print "The name of the file to process is: ", args.file
    
    #connect()

    f_out = create_output_file()
    
    drop_tables(f_out)
    create_tables(f_out)

    read_vcf(args.file,f_out)

    close_output_file(f_out)
   
    return 0


main()

