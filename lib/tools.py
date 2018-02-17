# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 16:27:06 2013

@author: Gael P. Alamancos
@email: gael.perez[at]upf.edu
"""

import sys
import logging
from abc import ABCMeta, abstractmethod
from lib.event import *


#Setting logging preferences
logging.basicConfig(level=logging.INFO)    
logger = logging.getLogger(__name__)


def nextel(my_object):
    """
    Returns next element of generator. Python2 and Python3 compatible.
    """
    if sys.version_info >= (3, 0):
        return my_object.__next__()
    else:
        return my_object.next()


def setToolsLoggerLevel(mode):
    try:
        logger.setLevel(eval(mode))
    except:
        logger.warning("Could not set up tools logger into %s level. Using INFO." %
            mode)


class FormatError(Exception):
    def __init__(self, line = None, *args):
        '''Exception raised for errors in the format of the input file.
        
        Keyword arguments:
        
        line            -- input line number where the error occured
        message (*args) -- accumulated errors explanation, i.e., tuple with the 
        explanations of the errors
        '''
        self.line = line        
        self.message = args  


class Parser(object):
    '''Abstract Class containing the common behaviour of all the parsers'''
    __metaclass__=ABCMeta
    
    def parseLine(self, line, lineNumber):
        '''General line parsing.
        
        Keyword arguments:
        
        line -- the input line
        lineNumber -- the line number of the input file
        
        Observations:
        
        The following attributes should be define in a subclass since it 
        depends on the filetype.
        
        self.MIN_FIELDS, self.SEQ_NAME, self.FEATURE_INDEX, self.START_INDEX, 
        self.END_INDEX, self.STRAND_INDEX, self.ATTR_INDEX 
        
        other attributes may be define if appropriate. 
        
        '''
        try:
            msg = []    #It will containg the error message acummulated.
            fields = line.rstrip('\n').split('\t')
            #Parsing the number of fields
            if len(fields) < self.MIN_FIELDS:
                msg.append("Unexpeced number of fields.")
            else:
                #Checking for empty or blank(\s) fields
                for empty in (i for i in range(self.MIN_FIELDS) if not \
                    fields[i] or fields[i].isspace()):                    
                    msg.append("Field %i empty or blank" % (empty+1))
                #Checking format of some fields
                if not msg:
                    try:
                        if int(fields[self.START_INDEX]) > int(
                            fields[self.END_INDEX]):
                            msg.append("Start coordinate greater than end "+ \
                                "coordinate.")
                    except ValueError:
                        msg.append("Start, end or both coordinates are not "+ \
                            "integers.")
                    if fields[self.STRAND_INDEX] not in ('+', '-', '.'):
                        msg.append("Unknown strand format in field: %i" % (
                            self.STRAND_INDEX+1))
            if msg: 
                raise FormatError(lineNumber, msg)
            return True
        except FormatError:
            logger.error("%s, in line %i. Skipping line..." % (sys.exc_info()[1].args[0], lineNumber + 1))
            return False
        except BaseException:
            logger.error("Unknown error: %s" % sys.exc_info()[1].args[0])
            sys.exit(1)


class IoeParser(Parser):
    def __init__(self):    
        '''It contains the correspondace of field in the IOE file starting 
            with 0.
        '''
        self.MIN_FIELDS = 5 #Number of fields
        self.SEQNAME_INDEX = 0 
        self.GENE_ID_INDEX = 1
        self.EVENT_ID_INDEX = 2
        self.ALT_ISO_INDEX = 3
        self.TOTAL_ISO_INDEX = 4
        
    def parseLine(self, line, lineNumber):
        try:
            msg = []    #It will containg the error message acummulated.
            fields = line.rstrip('\n').split('\t')
            #Parsing the number of fields
            if len(fields) < self.MIN_FIELDS:
                msg.append("Unexpeced number of fields.")
            else:
                #Checking for empty or blank(\s) fields
                for empty in (i for i in range(self.MIN_FIELDS) if not \
                    fields[i] or fields[i].isspace()):                    
                    msg.append("Field %i empty or blank" % (empty+1))
                #Checking format of some fields
                if not msg:
                    alt1 = set(fields[self.ALT_ISO_INDEX].split(','))
                    total = set(fields[self.TOTAL_ISO_INDEX].split(','))
                    if not total.issuperset(alt1):
                        msg.append ("Not all alternative_trancripts (field %i) present in the total_trancripts (field %i)" % (
                            self.ALT_ISO_INDEX+1, self.TOTAL_ISO_INDEX+1))
            if msg: 
                raise FormatError(lineNumber, msg)
            return True
        except FormatError:
            logger.error("%s, in line %i. Skipping line..." % (sys.exc_info()[1].args[0], lineNumber + 1))
            return False
        except BaseException:
            logger.error("Unknown error: %s" % sys.exc_info()[1].args[0])
            sys.exit(1) 
            
            
class ExpressionParser(Parser):
    def __init__(self):
        self.MIN_FIELDS = 2 #Number of fields
        self.TRANSCRIPT_ID_INDEX = 0
        self.EXPRESSION_INDEX = 1
        
#    def parseLine(self, line, lineNumber, tableFormat = False):
    def parseLine(self, line, lineNumber):
        '''Parse a expression line 
        
        keyword arguments:

        line -- line of an expression file.
        lineNumber -- line number to display in case of format errors.
        tableFormat -- true when there's more than one expression field.        
        '''
        try:
            msg = []    #It will containg the error message acummulated.
            fields = line.rstrip('\n').split('\t')
            #Parsing the number of fields
            if len(fields) < self.MIN_FIELDS:
                msg.append("Unexpeced number of fields. %i expected, %i given." % \
                    (self.MIN_FIELDS, len(fields)))
            else:
                #Checking for empty or blank(\s) fields
                for empty in (i for i in range(self.MIN_FIELDS) if not \
                    fields[i] or fields[i].isspace()):                    
                    msg.append("Field %i empty or blank" % (empty+1))
                
                #Checking format of all expression fields
                for i in range(1, self.MIN_FIELDS):
                    try:
                        float(fields[i])
                    except ValueError:
                        msg.append("Field %i is not a float." % (i+1))
                        break
            if msg: 
                raise FormatError(lineNumber, msg)
            return True
        except FormatError:
            logger.error("%s, in line %i. Skipping line..." % (sys.exc_info()[1].args[0], lineNumber + 1))
            return False
        except BaseException:
            logger.error("Unknown error: %s" % sys.exc_info()[1].args[0])
            sys.exit(1) 
            

class PsiParser(Parser):
    def __init__(self):
        self.MIN_FIELDS = 2 #Number of fields
        self.EVENT_ID_INDEX = 0
        self.PSI_INDEX = 1
        
#    def parseLine(self, line, lineNumber, tableFormat = False):
    def parseLine(self, line, lineNumber):
        try:
            msg = []    #It will containg the error message acummulated.
            fields = line.rstrip('\n').split('\t')
            #Parsing the number of fields
            if len(fields) < self.MIN_FIELDS:
                msg.append("Unexpeced number of fields.")
            else:
                #Checking for empty or blank(\s) fields
                for empty in (i for i in range(self.MIN_FIELDS) if not \
                    fields[i] or fields[i].isspace()):                    
                    msg.append("Field %i empty or blank" % (empty+1))
                
                #Checking format of all psi fields
                for i in range(1, self.MIN_FIELDS):
                    try:
                        float(fields[i])
                    except ValueError:
                        if fields[i] != "NA":
                            msg.append("Field %i is not a float or \"NA\"." % (i+1))
                            break                
                if msg: 
                    raise FormatError(lineNumber, msg)
            return True
        except FormatError:
            logger.error("%s, in line %i. Skipping line..." % (sys.exc_info()[1].args[0], lineNumber + 1))
            return False
        except BaseException:
            logger.error("Unknown error: %s" % sys.exc_info()[1].args[0])
            sys.exit(1)


class TpmParser(Parser):
    def __init__(self):
        self.MIN_FIELDS = 2  # Number of fields
        self.EVENT_ID_INDEX = 0
        self.PSI_INDEX = 1

    #    def parseLine(self, line, lineNumber, tableFormat = False):
    def parseLine(self, line, lineNumber):
        try:
            msg = []  # It will containg the error message acummulated.
            fields = line.rstrip('\n').split('\t')
            # Parsing the number of fields
            if len(fields) < self.MIN_FIELDS:
                msg.append("Unexpeced number of fields.")
            else:
                # Checking for empty or blank(\s) fields
                for empty in (i for i in range(self.MIN_FIELDS) if not \
                        fields[i] or fields[i].isspace()):
                    msg.append("Field %i empty or blank" % (empty + 1))

                # Checking format of all psi fields
                for i in range(1, self.MIN_FIELDS):
                    try:
                        float(fields[i])
                    except ValueError:
                        if fields[i] != "NA":
                            msg.append("Field %i is not a float or \"NA\"." % (i + 1))
                            break
                if msg:
                    raise FormatError(lineNumber, msg)
            return True
        except FormatError:
            logger.error("%s, in line %i. Skipping line..." % (sys.exc_info()[1].args[0], lineNumber + 1))
            return False
        except BaseException:
            logger.error("Unknown error: %s" % sys.exc_info()[1].args[0])
            sys.exit(1)


class Reader(object):
    '''Abstract Class with the common behaviour of all the readers.'''
    __metaclass__ = ABCMeta
    
    def openFile(self, fileName):
        '''It creates a pipe for reading the input file'''
        try:
            self._pipe = open(fileName, 'r')
            logger.info("File %s opened in reading mode." % fileName)
        except IOError:
            logger.error("I/O Error: %s" % (sys.exc_info()[1].args[1]))
            sys.exit(1)
        except BaseException:
            logger.error("Unknown error: %s" % sys.exc_info()[1].args[0])
            sys.exit(1)

    def readLine(self, header = False):
        '''This method requires to have a Reader sublass that inherits also from a 
        Parser class'''
        try:        
            for lineNumber, line in enumerate(self._pipe):
                #Parsing the line. Line skipped if there's any error. 
                #fields[0:8]
                if lineNumber == 0 and header: #Skip the header line
                    continue
                if line.startswith('#'):
                    logger.debug("Line %i starts with #. Skipping line..." % (
                    lineNumber+1))
                    continue
                if not self.parseLine(line, lineNumber):
                    continue
                else:
                    arguments = {}
                    fields = line.rstrip("\n").split("\t")
                    arguments["seqname"] = fields[self.SEQNAME_INDEX]
                    arguments["feature"] = fields[self.FEATURE_INDEX]
                    arguments["start"] = fields[self.START_INDEX]
                    arguments["end"] = fields[self.END_INDEX]
                    arguments["strand"] = fields[self.STRAND_INDEX]
                    attributes = self.parseAttributes(
                        fields[self.ATTR_INDEX], lineNumber,
                        arguments["feature"])
                    #Parsing attributes. Line skipped if there's any error. 
                    #fields[8]                  
                    if attributes == False:
                        continue
                    elif isinstance(attributes, dict):
                        arguments.update(attributes)
                    yield arguments
            self._pipe.close()
            logger.info("File %s closed." % self._pipe.name)
        except BaseException:
            logger.error("Unknown error: %s" % sys.exc_info()[1].args[0])
            sys.exit(1)

        
class IoeReader(Reader, IoeParser):
    def __init__(self):
        IoeParser.__init__(self)
        
    def readLine(self, header = True):
        '''This method reads and parses an ioe line.'''
        try:        
            for lineNumber, line in enumerate(self._pipe):
                if lineNumber == 0 and header: #Skip the header line
                    continue
                if line.startswith('#'):
                    logger.debug("Line %i starts with #. Skipping line..." % (
                    lineNumber+1))
                    continue
                if not self.parseLine(line, lineNumber):
                    continue
                else:
                    arguments = {}
                    fields = line.rstrip("\n").split("\t")
                    arguments["seqname"] = fields[self.SEQNAME_INDEX]
                    arguments["gene_id"] = fields[self.GENE_ID_INDEX]
                    arguments["event_id"] = fields[self.EVENT_ID_INDEX]
                    arguments["alt_iso"] = fields[self.ALT_ISO_INDEX]
                    arguments["total_iso"] = fields[self.TOTAL_ISO_INDEX]
                    yield arguments
            self._pipe.close()
            logger.info("File %s closed." % self._pipe.name)
        except BaseException:
            logger.error("Unknown error: %s" % sys.exc_info()[1].args[0])
            sys.exit(1)


class ExpressionReader(Reader, ExpressionParser):
    def __init__(self):
        ExpressionParser.__init__(self)
        
#    def readLine(self, header = True, tableFormat = False):
    def readLine(self, header = True):
        '''This method reads and parses an "expression" line'''
        try:
            min_fields = 0
            colIds = list() 
            for lineNumber, line in enumerate(self._pipe):
                fields = line.rstrip("\n").split("\t")
                if lineNumber == 0 and header: #Skip the header line
                    #Calculating the number of fields required
                    min_fields = (len(line.rstrip("\n").split("\t")) + 1)
                    #Storing column_id for the expression fields
                    colIds = fields
                    continue
                if line.startswith('#'):
                    logger.debug("Line %i starts with #. Skipping line..." % (
                        lineNumber+1))
                    continue
#                if tableFormat and min_fields > 2:
                if min_fields > 2:
                    #update the number of fields to request in each line
                    self.MIN_FIELDS = min_fields    
                #Parse the line
#                if not self.parseLine(line, lineNumber, tableFormat):
                if not self.parseLine(line, lineNumber):
                     continue
                else:
                    arguments = {}
                    arguments["transcript_id"] = fields[self.TRANSCRIPT_ID_INDEX]
                    arguments["expression"] = {}
                    arguments["colIds"] = colIds
                    for x in range(1, self.MIN_FIELDS):
                        arguments["expression"][colIds[(x-1)]] = fields[x]                        
                    yield arguments
            self._pipe.close()
            logger.info("File %s closed." % self._pipe.name)    
        except BaseException:
            logger.error("Unknown error: %s" % sys.exc_info()[1].args[0])
            sys.exit(1)


class FactoryReader(object):
    def getReader(self, fileFormat):
        if fileFormat == "ioe":
            return IoeReader()
        elif fileFormat == "expression":
            return ExpressionReader()
        else:
            logger.error("Unknown reader:" + fileFormat)
            sys.exit(1)


class Writer(object):
    '''Abstract class that encapsulates all the common behaviour of the 
    writers
    '''
    __metaclass__ = ABCMeta
    
    def openFile(self, outputFile):
        if (not os.path.exists(os.path.dirname(outputFile))) and (os.path.dirname(outputFile) != '') :
            logger.error('Can not save to file "{}". Path does not exist.'.format(outputFile))
            sys.exit(1)
        try:
            self._pipe = open(outputFile, 'w')
            self.lineNumber = 0    
        except BaseException:
            logger.error("Unknown error: %s" % sys.exc_info()[1].args[0])
            sys.exit(1)

    def writeLine(self, line, parse = True):
        self.write(line+"\n", parse)
        
    def write(self, line, parse = True):
        '''This method requires to have a Writer sublass that inherits also 
        from a Parser class'''
        try:
            if not line.startswith('#') and parse:
                if self.parseLine(line, self.lineNumber):
                    self._pipe.write(line)    
                    self.lineNumber += 1
            else:
                self._pipe.write(line)
                self.lineNumber += 1
        except BaseException:
            logger.error("Unknown error: %s" % sys.exc_info()[1].args[0])
            sys.exit(1)

    @staticmethod
    def getWriter(fileFormat):
        try:
            if fileFormat == "GTF" or fileFormat == "IOE":
                logger.error("Error: These GTF and IOE writers no longer exist.")
                sys.exit(1)
            elif fileFormat == "PSI":
                return PsiWriter()
            elif fileFormat == "TPM":
                return TpmWriter()
            else:
                logger.error("Unknown writer:" + fileFormat)
                sys.exit(1)
        except BaseException:
            logger.error("Unknown error: %s" % sys.exc_info()[1].args[0])
            sys.exit(1)

    def closeFile(self):
        try:
            self._pipe.close()
        except BaseException:
            logger.error("Unknown error: %s" % sys.exc_info()[1].args[0])
            sys.exit(1)


class PsiWriter(Writer, PsiParser):
    '''It conatains all the utilities required to write an PSI file.'''
    def openFile(self, outputFile):
        outputFile += ".psi"
        Writer.openFile(self, outputFile)
        
    @staticmethod
    def lineGenerator(key, value, colIds):
        '''Generates a PSI line from a single event
        
        key -- id of the event.
        value -- the dictionary containning all column names and psi values calculated.
        colIds -- list containing the column names
        '''
        line = []
        line.append(key)
        for x in colIds:
            line.append(str(value[x]))
        return line


def ioi_writer(genome_obj, output_file):

    with open(output_file+".ioi", "w+") as fh:
        fh.write("seqname\tgene_id\tisoform_id\tinclusion_transcripts\ttotal_transcripts\n")
        for obj in genome_obj:
            for trans_id in obj[0].transcripts.keys():
                chrm = obj[1]
                gene_id = obj[0].name
                all_trans_ids = ",".join(obj[0].sortedTranscripts)
                line = "%s\t%s\t%s\t%s\t%s\n" % \
                       (chrm, gene_id, gene_id+";"+trans_id, trans_id, all_trans_ids)
                fh.write(line)

                
class TpmWriter(Writer, TpmParser):
    '''It conatains all the utilities required to write an PSI file.'''

    def openFile(self, outputFile):
        outputFile += ".tpm"
        Writer.openFile(self, outputFile)

    @staticmethod
    def lineGenerator(key, value, colIds):
        '''Generates a PSI line from a single event

        key -- id of the event.
        value -- the dictionary containning all column names and psi values calculated.
        colIds -- list containing the column names
        '''
        line = []
        line.append(key)
        for x in colIds:
            line.append(str(value[x]))
        return line
