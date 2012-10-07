#Parser for Sequence files

import re
import time

class SequenceError(Exception):
    '''Syntax Error in Sequence. Prohibits code execution.'''
    def __init__(self, seq, position, info):
        self.position = position
        self.info = info
        self.seq = seq
        
    def __str__(self):
        #Print error message e
        message = 'Sequence Error: {0}, Position: {1}\n'.format(self.info, self.position)
        #Print annotated sequence
        before = self.seq[:self.position]
        after = self.seq[self.position+1:]
        max_len = 15
        if len(before) > max_len:
            before = str(len(before)-max_len) + ' ... ' + before[-max_len:]
        if len(after) > max_len:
            after = after[:max_len] + ' ... ' + str(len(after)-max_len)
            
        message +=  before + '<' + self.seq[self.position] + '>' + after + '\n'
        
        return message
    
class SequenceWarning(Exception):
    '''Syntax Error in Sequence. Allows sequence execution.'''
    def __init__(self, sequence, warnings):
        self.sequence = sequence
        self.warnings = warnings
        
    def __str__(self):
        self.message = 'Sequence Warnings:\n'
        for warning in self.warnings:
            self.message += '{0}, Position: {1}\n'.format(warning[0], warning[1])
            before = self.sequence[:warning[1]]
            after = self.sequence[warning[1]+1:]
            max_len = 15
            if len(before) > max_len:
                before = str(len(before)-max_len) + ' ... ' + before[-max_len:]
            if len(after) > max_len:
                after = after[:max_len] + ' ... ' + str(len(after)-max_len)
                
            self.message += before + '<' + self.sequence[warning[1]] + '>' + after + '\n'
            
        return self.message

class Parser:        
    def load_file(self, filename):
        #Load the file
        fil = open(filename)
        content = ''
        for line in fil:
            content += line
        fil.close()
        self.sequence = content
        
    def save_file(self, filepath, names, verbose):
        
        if verbose == True:
            
            to_save = '# Sequence Processing File\n\n# Generated {0}\n\n# This sequence can be modeled using SeqGUI.py\n'.format(time.asctime())
            
            if len(names) > 0:
                to_save += '\n# Original names for sequences: '
                not_first = False
                for name in names:
                    if not_first:
                        to_save += ' & '
                    to_save += name
                    not_first = True
            
            to_save += '\n'
            
            descriptions = {}
            unit_types = {'N': 'Translation','R': 'Restriction','L': 'Ligation','F': 'Transcription Factor',
                          'M': 'Mutation', 'E': 'Sequence Export','I': 'Sequence Import','O': 'Copy',
                          'Q': 'Sequence Protection','D': 'Sequence Degradation','Z': 'Cell Disintegration','U': 'Tag Transport',
                          'V': 'Substance Conversion','Y': 'Cell Splitting'}
            
            #Create regular expression pattern
            #matching for unit description
            pattern = r'([NRLFMEIOQDZUVY])'               #unit type [0]
            pattern += r'([\-]*\d+[\.]?\d*)'            #max action parameter (can be floating point) [1]
            pattern += r'/?(\d*[\.]?\d*)'               #optional energy parameter after / (can be floating point) [2]
            pattern += r'\(?([a-z_,\*]*)\)?'              #optional tags in () [3]
            pattern += r'\[[GTACXP0-9]*\]'              #contained sequence in [] 
            pattern += r'([GTACX\!]*)'                  #optional target sequence [4]
            pattern += r'\(?([a-z_,\*]*)\)?'              #optional tags targets [5] (or starting substances for conversion unit V)
            pattern += r'\(?([a-z_,\*]*)\)?'              #optional tags resulting substances (only used by conversion unit) [6]
            pattern += r'([0-9]+)'                      #life parameter [7]
            
            regex = re.compile(pattern)
            new_units = regex.finditer(self.sequence)
            
            for unit in new_units:
                params = unit.groups()
                desc = '\n\n#{0} Unit, a: {1}, l: {2}\n'.format(unit_types[params[0]], params[1], params[7])
                end_pos = unit.start() + 2
                for i in range(len(params)):
                    end_pos += len(params[i])
                    
                descriptions[unit.start()] = desc
                descriptions[end_pos] = '\n'
                
            for pos in range(len(self.sequence)):
                if pos in descriptions.keys():
                    to_save += descriptions[pos]
                if self.sequence[pos] == ';':
                    to_save += '\n'
                to_save += self.sequence[pos]
        else:   
            to_save = self.sequence
                                                                                    
        try:
            # This will create a new file or **overwrite an existing file**.
            f = open(filepath, "w")
            try:
                f.write(to_save) # Write a string to a file
            finally:
                f.close()
        except IOError, e:
            print e
                
    def clean_sequence(self):
        '''Splits lines, removes whitespace, newlines and comments from sequence.'''  
        #Remove whitespaces and newlines and comments from program
        seq = self.sequence.splitlines()
        sequence = ''
        for sequ in seq:
            #remove comments
            if len(sequ) > 0:
                if sequ[0] != '#':
                    sequence += sequ
        
        sequence = re.sub(r'\s', '', sequence)
        
        self.sequence = sequence
        
    def check_syntax(self):
        '''Crude syntax check'''
        
        #Check if semicolon syntax is correct
        #if any of these patterns matches raises Sequence error_type 'semicolon'
        if self.sequence[0] != ';':
            raise SequenceError(self.sequence, 0, 'Missing semicolon at start.')
        
        patterns = ['[^;0-9][0-9]*(;)[0-9]+[^;0-9]']
        result = re.compile(patterns[0]).search(self.sequence, 1)
        if result != None:
            raise SequenceError(self.sequence, result.start(), 'Missing semicolon at overhanging end.')
        
        result = re.compile('[^;0-9]([0-9]+)$').search(self.sequence, 1)
        if result != None:
            position = result.start() + len(result.groups()[0])
            raise SequenceError(self.sequence, position, 'Detected illegal end of sequence.')
        
        result = re.compile('([^;0-9])$').search(self.sequence, 1)
        if result != None:
            raise SequenceError(self.sequence, result.start(), 'Missing semicolon at end of sequence.')
        
        #detect illegal multiples of semicolon
        patterns = [';;;', ';;$']
        for pattern in patterns:
            result = re.compile(pattern).search(self.sequence, 1)
            if result != None:
                raise SequenceError(self.sequence, result.start(), 'Too many semicolons.')
                
        if self.sequence[:2] == ';;':
            raise SequenceError(self.sequence, 0, 'Too many semicolons.')
        
        #detect illegal floating point number
        patterns = [
                    ';\d+\.\d+;', ';;\d+\.\d+', ';\d+\.\d+', 
                    '[^NMRLQFIEODPZUV\d\-/]\d+\.\d+',
                    '\.\.', '\.[^\d]'
                    ]
        for pattern in patterns:
            result = re.compile(pattern).search(self.sequence, 1)
            if result != None:
                raise SequenceError(self.sequence, result.start(), 'Floating point number only allowed as max action- or promoter-parameter.')
        
        #Check if bracket [ syntax is correct
        #raises SequenceError_type 'bracket'
        result = re.compile('[^0-9\)](\[)').search(self.sequence, 1)
        if result != None:
            raise SequenceError(self.sequence, result.start(), 'Bracket preceded by other than number or ).')
        
        #Check if number of '[' matches number of ']'
        count = 0
        for pos in range(len(self.sequence)):
            if self.sequence[pos] == '[':
                count += 1
            elif self.sequence[pos] == ']':
                count -= 1
                
        if count > 0:
            raise SequenceError(self.sequence, len(self.sequence)-1, 'Too many open brackets. Remove \'[\'.') 
        elif count < 0:
            raise SequenceError(self.sequence, len(self.sequence)-1, 'Too many closed brackets. Remove \']\'.')
        
        #Check if illegal characters are in sequence using whitelist of allowed characters
        result = re.compile(r'([^0-9a-zPSGATCXNDFMQRLEIOZUVY,_\*;.!\)\(\[\]\-/])').search(self.sequence, 1)
        if result != None:
            message = 'Found illegal character: ' + result.groups(0)[0] + '.'
            raise SequenceError(self.sequence, result.start(), message)
        
        #Check tag syntax
        patterns = ['[^0-9GTCAX\!\]\)]\(', '[^a-z_\*\(]\)', '\)[^\[\(\d]']
        for pattern in patterns:
            result = re.compile(pattern).search(self.sequence, 1)
            if result != None:
                raise SequenceError(self.sequence, result.start(), 'Invalid tag syntax.')
            
        complete_tags = [(m.groups()[0], m.start()) for m in re.compile('\(([^\)]*)\)').finditer(self.sequence)]
        for tag in complete_tags:
            result = re.compile('[^a-z_,\*]').search(tag[0], 1)
            if result != None:
                raise SequenceError(self.sequence, (tag[1] + result.start()+1), 'Character not allowed in tag.')
        
        #Check if no parameter exceeds the three digit limit
        result = re.compile('[^.](\d{4})').search(self.sequence, 1)
        if result != None:
            raise SequenceError(self.sequence, result.start() + 1, 'Parameter exceeds three digit limit.')
        
        #Check if forward slash is surrounded by digits
        result = re.compile('[^\d]/').search(self.sequence, 1)
        res = re.compile('/[^\d]').search(self.sequence, 1)
        if result != None or res != None:
            raise SequenceError(self.sequence, result.start() + 1, 'Forward slash has to be surrounded by digits.')
        
        #Check if no parameter stands after stop codon or [
        result = re.compile('[S\[]([0-9])').search(self.sequence, 1)
        if result != None:
            raise SequenceError(self.sequence, result.start(), 'Parameter not allowed at this position.')
        
        #Negative parameters are only allowed for F,U,I or E action parameter
        patterns = ['[^FUIE\!](\-)', '(\-)[^0-9]', '[^FUIE]\!(\-)']
        for pattern in patterns:
            res0 = re.compile(pattern).search(self.sequence, 1)
            if res0 != None:
                raise SequenceError(self.sequence, res0.start(), 'Negative parameters only allowed for transcription factor (F) or transport unit (U) maximal action.')
        
        #Check if ! syntax is correct
        result = re.compile('([SNMQDLNIOEZRUVFPY\[])\!').search(self.sequence, 1)
        if result != None:
            raise SequenceError(self.sequence, result.start(), '! not allowed after {0}.'.format(result.groups(0)[0]))
        
        result = re.compile('\!([\[\]\(])').search(self.sequence, 1)
        if result != None:
            raise SequenceError(self.sequence, result.start(), '! not allowed before {0}.'.format(result.groups(0)[0]))
        
        result = re.compile('\d(\!)\d').search(self.sequence, 1)
        if result != None:
            raise SequenceError(self.sequence, result.start(), '! between parameters not allowed.')
        
        result = re.compile('(\!\!)').search(self.sequence, 1)
        if result != None:
            raise SequenceError(self.sequence, result.start(), 'Pattern \'!!\' not allowed.')
        
        
        #find all separate sequences
        pattern = ';([0-9]*)([^;]+);([0-9]*)'
        regex = re.compile(pattern)
        sequences = [(m.groups(), m.start()) for m in regex.finditer(self.sequence)]
        
        #Check for illegal overhanging ends
        for sequence in sequences:
            tot_overlap = 0
            if sequence[0][0]:
                tot_overlap += int(sequence[0][0])
            if sequence[0][2]:
                tot_overlap += int(sequence[0][2])
            if len(sequence[0][1]) < tot_overlap:
                raise SequenceError(self.sequence, sequence[1] + 1, 'Sequence has too much overlap.')
        
        # Check for sequence warnings (allow simulation of sequence)
        # Therefore only one SequenceWarning is raised listing all the
        # occurring warnings.
        
        warnings = []
        
        #Create regular expression pattern
        #matching for unit description
        pattern = r'([NRLFMEIOQDZUVY])'               #unit type [0]
        pattern += r'([\-]*\d+[\.]?\d*)'            #max action parameter (can be floating point) [1]
        pattern += r'/?(\d*[\.]?\d*)'               #optional energy parameter after / (can be floating point) [2]
        pattern += r'\(?([a-z_,\*]*)\)?'              #optional tags in () [3]
        pattern += r'\[[GTACXP0-9]*\]'              #contained sequence in [] 
        pattern += r'([GTACX\!]*)'                  #optional target sequence [4]
        pattern += r'\(?([a-z_,\*]*)\)?'              #optional tags targets [5] (or starting substances for conversion unit V)
        pattern += r'\(?([a-z_,\*]*)\)?'              #optional tags resulting substances (only used by conversion unit) [6]
        pattern += r'([0-9]+)'                      #life parameter [7]
        
        regex = re.compile(pattern)
        new_units = regex.finditer(self.sequence)
        
        for unit in new_units:
            g = unit.groups()
            if g[0] in ('Z', 'N', 'M', 'Y'):
                if len(g[4]) != 0:
                    msg = g[0] + ' unit with target'
                    warnings.append((msg, unit.start()))
                    break
                    
            if g[0] in ('O','E','F','R','L','I','Q','D'):
                if len(g[4]) == 0:
                    msg = g[0] + ' unit without target'
                    warnings.append((msg, unit.start()))
                    break
                    
                if g[0] in ('I','Q','D'):
                    try:
                        g[4].index('!')
                        msg = g[0] + ' unit with ! in target'
                        warnings.append((msg, unit.start()))
                        break
                        
                    except:
                        pass
                
                if g[0] in ('O','E','F','R','L'):
                    try:
                        g[4].index('!')
                    except:
                        msg = g[0] + ' unit without ! in target'
                        warnings.append((msg, unit.start()))
                        break
                    
            if g[0] == 'U':
                if not g[5]:
                    msg = 'U unit without tag target'
                    warnings.append((msg, unit.start()))
                    break
        
        #raise sequence warning
        if warnings:
            raise SequenceWarning(self.sequence, warnings)
        
    def sequence_info(self):
        '''!!!FIX: This is very broken...
        Returns a formated string containing information
        about the content of current sequence.'''
        
        #Find all units in sequence
        pattern = r'([NRLFQMEIDOZUV])([\-]*[0-9]+)\[([GTACXP0-9]*)\]([GTACX\!]*)([0-9]+)'
        units = re.compile(pattern).finditer(self.sequence)

        unit_dict = {}
        for unit in units:
            #Save parameters
            pos = unit.start()
            unit_target = []
            unit_type = unit.groups()[0]
            target = unit.groups()[3]
            if unit_type in ('R', 'O', 'E', 'F', 'Q'):
                target_pos = target.index('!')
                if unit_type == 'F':
                    target = '(P[0-9]*)' + target
                #Find target positions
                target = '[^/!]' + target.replace('!', '')
                target = target.replace('X', '[GATCX]') + '[^\!]'
                
                results = re.compile(target).finditer(self.sequence)
                for result in results:
                    if unit_type == 'F':
                        unit_target.append(result.start() + 1 + target_pos + len(result.groups()[0]))
                    else:
                        unit_target.append(result.start() + 1 + target_pos)
                
            elif unit_type == 'L':
                pass
                        
            unit_dict[pos] = [unit_type, int(unit.groups()[1]), unit.groups()[2], unit.groups()[3], int(unit.groups()[4]), unit_target]
                
        #Find all promoters and 
        prom_dict = {}
        pattern = 'P([0-9]*)([^S;]*)'
            
        promoters = re.compile(pattern).finditer(self.sequence)
            
        for promoter in promoters:
            start = promoter.start()
            end = start + 1 + len(promoter.groups()[0]) + len(promoter.groups()[1])
            try:
                num = int(promoter.groups()[0])
            except:
                num = 0
                
            prom_dict[start] = [end, num]
        
        resp = ''
        for pos in range(len(self.sequence)):
            resp += '\n' + str(pos) + ': ' + self.sequence[pos]
            target_pos = []
            for pos1 in unit_dict.keys():
                if pos in unit_dict[pos1][5]:
                    target_pos.append(pos1)
                    
            if pos in prom_dict or pos in unit_dict or len(target_pos) > 0:
                resp += ' ==> '
            if pos in prom_dict:
                found_prom = True
                if prom_dict[pos][1] != 0: 
                    resp += 'Promoter ({0} activations/round) until position {1}'.format(prom_dict[pos][1], prom_dict[pos][0])
                else:
                    resp += 'Promoter until position ' + str(prom_dict[pos][0])
            else:
                found_prom = False
                
            if pos in unit_dict:
                found_unit = True
                if found_prom:
                    resp += ' & '
                if unit_dict[pos][0] == 'R': 
                    resp += 'Restriction Unit'
                elif unit_dict[pos][0] == 'L':
                    resp += 'Ligation Unit'
                elif unit_dict[pos][0] == 'O':
                    resp += 'Copy Unit'
                elif unit_dict[pos][0] == 'I':
                    resp += 'Import Unit'
                elif unit_dict[pos][0] == 'E':
                    resp += 'Export Unit'
                elif unit_dict[pos][0] == 'M':
                    resp += 'Mutation Unit'
                elif unit_dict[pos][0] == 'F':
                    resp += 'Transcription Factor'
                elif unit_dict[pos][0] == 'N':
                    resp += 'Translation Unit'
                elif unit_dict[pos][0] == 'Q':
                    resp += 'Sequence Protection Unit'
                    
                resp += ', Actions/Round: {0}, Life: {1}, Target: {2}, Target Positions: {3}'.format(unit_dict[pos][1], unit_dict[pos][4], unit_dict[pos][3], unit_dict[pos][5])
            
            else:
                found_unit = False
            
            
            found_target = False 
            
            if len(target_pos) > 0:
                if found_prom or found_unit:
                    resp += ' & '
                for target in target_pos:
                    if found_target:
                        resp += ' & '
                    resp += 'Target of {1} Unit @ Position {0}'.format(target, unit_dict[target][0])
                    found_target = True
                
        return resp