import re
import random
import csv
from collections import Counter

class Cell(object):
    def __init__(self, sequence = '', name = '', ind = 0, 
                 sample_mode='standard', normal_factor = 0.2, 
                 units = [], substances = [], cells = [],
                 is_top_cell = True):
        
        self.is_top_cell = is_top_cell
        self.id, self.sequence, self.name = ind, sequence, name
        self.rounds = 0
        self.pr_pos = []
        self.target_sequences, self.partner_target_sequences = [], []
        
        self.partner_sequence = ''
        self.partner_units = {}
        self.partner_substances = {}
        
        self.log = {}
        
        self.interact_with_subcells = True
        
        #substances which will be added every round
        self.subs_const = {}
        
        #max action and promoter parameter sample mode
        self.sample_mode = sample_mode
        self.normal_factor = normal_factor
        
        self.cells = {}
        self.units = {}
        self.substances = {}
        
        self.unit_count, self.cell_count, self.subs_count = 0, 0, 0
        
        #add previous units, substances and subcells
        if units:
            for unit in units:
                self.add_unit(unit_dict = unit)
        if substances:
            self.add_subs(tuple(substances))
        if cells:
            for c in cells:
                self.new_sub_cell(c.sequence, c.name, units = c.units.values(), substances = c.substances.values(), cells = c.cells.values())
        
        self.reset_target()
        self.write_log()
        
    def new_sub_cell(self, sequence = '', name = '', units = [], substances = [], cells = []):
        '''Creates a new sub cell in the current cell from a sequence.
        Returns the ID of the new sub_cell'''
        
        #Create a new instance of the cell class
        self.cells[self.cell_count] = Cell(sequence = sequence, name = name, ind = self.cell_count, 
                                           units = units, substances = substances, cells = cells, 
                                           is_top_cell = False)
        
        #update cell counter
        self.cell_count += 1
        
        return self.cell_count - 1
    
    def add_unit(self, unit_dict):
        '''add a unit by dictionary'''
        self.units[self.unit_count] = unit_dict
            
        self.unit_count += 1
    
    def get_subcell_by_adress(self, adress):
        'Returns False if not successful'
        try:
            adress = list(adress)
            if len(adress) > 1:
                adress.pop(0)
                sub_cell = adress[0]
                #redirect call to subcell
                answer = self.cells[sub_cell].get_subcell_by_adress(tuple(adress))
                return answer
            else:
                return self
            
        except Exception, e:
            return False
        
    def get_cell_sequence(self, adress):
        cell = self.get_subcell_by_adress(adress)
        return cell.sequence
    
    def set_cell_sequence(self, adress, sequence):
        if adress == (self.id,):
            self.sequence = sequence
        else:
            adress = list(adress)
            adress.pop(0)
            self.cells[adress[0]].change_parameters(tuple(adress), sequence)
    
    def delete_sub_cell(self, adress):
        '''Delete a sub cell by id, Return True if successful and False if not.'''
        if adress[0] == self.id and len(adress) == 2:
            #find cell_id in current cells
            try:
                del self.cells[adress[1]]
                return True
            except:
                return False
        
        elif len(adress) < 2:
            return False
        
        elif len(adress) > 2:
            adress = list(adress)
            adress.pop(0)
            answer = self.cells[adress[0]].delete_sub_cell(adress = tuple(adress))
            return answer
        
        #i think this never triggers unless you use invalid adresses...
        return False
    
    def change_energy(self, amount):
        '''changes energy by amount, returns False if impossible'''
        en_id = [key for key, sub in self.substances.items() if sub == 'energy']
        if amount < 0:
            amount = abs(amount)
            if len(en_id) >= amount:
                for ind in en_id[:amount]:
                    del self.substances[ind]
                return True
            else:
                return False
        else:
            for i in range(amount):
                self.add_subs(u'energy')
            return True
    
    def has_energy(self, amount = 1, total = False):
        en_id = [key for key, sub in self.substances.items() if sub == 'energy']
        if len(en_id) >= amount:
            if total:
                return len(en_id)
            return True
        else:
            if total:
                return len(en_id)
            return False
            
    def change_parameters(self, adress=False, sample_mode= False, sigma_factor=False, subs_const=False, cont = False, interact_with_subcells = False):
        if adress == (self.id,) or adress == False:
            #change all requested parameters
            if sigma_factor:
                self.normal_factor = sigma_factor
            if sample_mode:
                self.sample_mode = sample_mode
            if subs_const:
                self.subs_const = subs_const
            if interact_with_subcells:
                if interact_with_subcells == 'False':
                    self.interact_with_subcells = False
                elif interact_with_subcells == 'True':
                    self.interact_with_subcells = True
            
            if cont:
                for cell_id in self.cells.keys():
                    self.cells[cell_id].change_parameters(False, sample_mode, sigma_factor, subs_const, True, interact_with_subcells)
        else:
            adress = list(adress)
            adress.pop(0)
            self.cells[adress[0]].change_parameters(tuple(adress), sample_mode, sigma_factor, subs_const, cont, interact_with_subcells)
            
    def add_subs(self, substances):
        if type(substances) == tuple:
            for substance in substances:
                self.substances[self.subs_count] = substance
                self.subs_count += 1
        elif (type(substances) == str) or (type(substances) == unicode):
            self.substances[self.subs_count] = substances
            self.subs_count += 1
            
    def find_tags(self, tags, substances = True, units = False, remove = False, use_partner = False):
        '''First checks if all the substances are present,
        then removes the substances from the list. 
        Returns True if successful and False if not.
        also returns True when empty tuple is passed...'''
        
        subs_found = []
        units_found = []
        objects_to_return = []
        
        if use_partner:
            mixed_subs = self.partner_substances.items()
            mixed_units = self.partner_units.items()
        else:
            mixed_subs = self.substances.items()
            mixed_units = self.units.items()
            
        random.shuffle(mixed_subs)
        random.shuffle(mixed_units)
        
        for tag in tags:
            found = False
            tag = tag.replace('*', '[a-z_\*]')
            tag = '^' + tag + '$'
            reg = re.compile(tag)
            
            #Look for substances
            if substances:
                for ind in range(len(mixed_subs)):
                    key, sub = mixed_subs[ind]
                    res = reg.search(sub)
                    if res:
                        subs_found.append(key)
                        objects_to_return.append((key, sub))
                        del mixed_subs[ind]
                        found = True
                        break
                     
            #if substance is not found search for unit
            if found == False and units == True:
                for ind in range(len(mixed_units)):
                    key, unit = mixed_units[ind]
                    for sub_tag in unit['tags']:
                        res = reg.search(sub_tag)
                        if res:
                            units_found.append(key)
                            objects_to_return.append((key, unit))
                            del mixed_units[ind]
                            found = True
                            break
                    if found:
                        break
                
        if len(tags) == len(units_found) + len(subs_found):
            if remove and not use_partner:
                for sub in subs_found:
                    del self.substances[sub]
                for unit in units_found:
                    del self.units[unit]   
            return tuple(objects_to_return)
        else:
            return False
        
    def run_all(self):
        '''Run actions from the cell and all sub cells in random order'''
        
        #Get action list from all sub cells & own actions
        action_list = self.start_round()      
        #Shuffle all actions:
        random.shuffle(action_list)
            
        for action in action_list:
            #Top cell may only interact with a mirror of itself
            self.perform_action(action, '', {}, {})
                
        self.degrade_units()
        self.write_log(cont = True)
        
    def get_requested_reads(self):
        '''Find all positions that are requested for translation this round.'''
        requests, neg_requests = [], []
        #find number of possible translations by adding up all actions from all translation units
        tot_tr = 0
        tr_list = [id for id in self.units.keys() if self.units[id]['type'] == 'N']
        tf_list = [id for id in self.units.keys() if self.units[id]['type'] == 'F']
        
        for tr in tr_list:
            tot_tr += self.sample_probability(self.units[tr]['action'])
        
        if tot_tr > 0:
            #Translations
            #Find all requestet translations
            for unit in tf_list:
                max_action = self.sample_probability(self.units[unit]['action'])
                
                #The target sequence without !
                target = self.units[unit]['target'].replace('!', '')
                
                occurences = self.find_occurences('prom', 0, len(self.sequence), target, abs(max_action))
                
                for occ in occurences:
                    if max_action > 0:
                        requests.append(occ)
                    else:
                        neg_requests.append(occ)
            
            #Find constitutive promoters
            pattern = 'P([0-9]+\.?[0-9]*)'
            
            regex = re.compile(pattern)
            occurences = [(m.start(), m.groups()[0]) for m in regex.finditer(self.sequence)]
            
            for occ in occurences:
                tr_start = occ[0] + len(occ[1]) + 1
                #add found const prom to request list
                max_tr = self.sample_probability(float(occ[1]))
                for i in range(max_tr):
                    requests.append(tr_start)
                    
            #remove negative regulated positions from the positive list if they are present
            for req in neg_requests:
                if req in requests:
                    ind = requests.index(req)
                    requests.pop(ind)
                    
            # pop starts that are protected
            todel = []
            for i in range(len(requests)):
                if self.pos_protected(requests[i]):
                    todel.append(requests[i])
            
            for tod in todel:
                requests.pop(requests.index(tod))
                   
            random.shuffle(requests)
            if len(requests) > tot_tr:
                requests = requests[:tot_tr]
                
        return requests
        
    def start_round(self):
        '''Initializes the round, withdrawing life from units.
        returns a list of actions that are requested this round'''
        
        self.rounds += 1
        
        #reset target sites
        self.reset_target()
        
        #reset protected sites
        self.reset_protected()
        
        #add constant substances to cell (also energy)
        new_tags = []
        for tag,amount in self.subs_const.items():
            for i in range(amount):
                new_tags.append(tag)
        self.add_subs(tuple(new_tags))
        
        #Withdraw life from units
        for unit in self.units.keys():
            self.units[unit]['life'] -= 1
        
        action_list = []
        
        #Collect all actions for most units and new_translations in random order
        #Do NOT use transcription factors, protection units, translation units
        #as they are handled separately.
        for unit_id in self.units.keys():
            unit_type = self.units[unit_id]['type']
            if unit_type not in ('F', 'N', 'Q'):
                max_action = self.sample_probability(self.units[unit_id]['action'])
                for i in range(abs(max_action)):
                    action_list.append([[self.id], unit_type, unit_id])
        
        #Get all requested translations
        read_ind = self.get_requested_reads()
        if read_ind:
            for ind in read_ind:
                action_list.append([[self.id], 'read', ind])
        
        #update dict of all translation units and their actions
        tr_list = [id for id in self.units.keys() if self.units[id]['type'] == 'N']
        self.tr_count = {}
        for tr in tr_list:
            self.tr_count[tr] = self.units[tr]['action']
        
        #start round in ALL sub cells
        for cell in self.cells.keys():
            sub_action_list = self.cells[cell].start_round()
            for adress, act_type, param in sub_action_list:
                adress.insert(0, self.id)
                action_list.append([adress, act_type, param])
        
        return action_list
    
    def perform_action(self, action, top_sequence, top_units, top_subs):
        '''Performs an action or redirects it to sub cell'''
        self.transported = {'unit_export': None,
                            'unit_import': None,
                            'subs_export': None,
                            'subs_import': None,
                            'seq_export': None,
                            'seq_import': None,
                            'partner_cell': -1}
        
        adress= action[0]
        a_type = action[1]
        param = action[2]
        
        if adress != [self.id] and a_type not in ('Y', 'Z'):
            #remove own cell id from adress
            adress.pop(0)
            action = [adress, a_type, param]
            
            seq = ''
            subs = {}
            units = {}
            #parameters are only required if action is happening in next cell
            #and a transport can happen
            if len(adress) == 1:
                if a_type in ('E', 'I'):
                    seq = self.sequence
                if a_type == 'U':
                    units = self.units
                    subs = self.substances
                    
            #redirect action to lower layer
            try:
                self.cells[adress[0]].perform_action(action, seq, units, subs)
                tr = self.cells[adress[0]].transported
            except (KeyError, AttributeError), e:
                return
            
            #Perform all transports in into itself
            if tr['partner_cell'] == -1:
                if tr['seq_export'] != None:
                    self.sequence += tr['seq_export']
                if tr['seq_import'] != None:
                    pos = tr['seq_import']
                    self.sequence = self.sequence[:pos[0]] + self.sequence[pos[1]:]
                if tr['unit_export'] != None:
                    self.add_unit(tr['unit_export'])
                if tr['unit_import'] != None:
                    del self.units[tr['unit_import']]
                if tr['subs_export'] != None:
                    self.add_subs(tr['subs_export'])
                if tr['subs_import'] != None:
                    del self.substances[tr['subs_import']]
            
        #these actions are conducted by mother cell (one layer up...)
        elif action[0][0] == self.id and action[1] in ('Y', 'Z') and len(action[0]) == 2:
            try:
                energy_used = self.cells[action[0][1]].units[action[2]]['energy']
                energy_used = self.sample_probability(energy_used)
            except:
                #if the unit was moved during this round 
                #it will no longer perform any actions.
                return False
            
            if action[1] == 'Y':
                self.split_sub_cell(tuple(action[0]), energy_used)
                return True
            elif action[1] == 'Z':
                self.disintegrate_sub_cell(tuple(action[0]), energy_used)
                return True
                
        #action happens in this cell
        else:
            if self.interact_with_subcells:
                #randomly choose among top cell and all subcells
                choices = self.cells.keys()
                choices.append(-1)
                cell_id = random.choice(choices)
            else:
                cell_id = -1
            
            self.transported['partner_cell'] = cell_id
            
            #change the partner sequence, units & substances if needed
            if cell_id == -1:
                if a_type in ('E', 'I'):
                    self.partner_sequence = top_sequence
                if a_type == 'U':
                    self.partner_units = top_units
                    self.partner_substances = top_subs
                    
            else:
                if a_type in ('E', 'I'):
                    self.partner_sequence = self.cells[cell_id].sequence
                if a_type == 'U':
                    self.partner_units = self.cells[cell_id].units
                    self.partner_substances = self.cells[cell_id].substances
            
            try:
                needs = self.units[param]['needs']
                if not self.find_tags(needs, units=True):
                    return
            except KeyError:
                pass
            
            #reset target sites
            self.reset_target()
            
            num_act, energy_used = 0, 0
            
            if a_type not in ('read','Y','Z'):
                try:
                    energy_used = self.sample_probability(self.units[param]['energy'])
                except KeyError:
                    #the unit is not in the cell anymore therefore
                    #no action can be performed
                    return
                
                if  self.has_energy(total = True) >= energy_used:
                    #Perform action, num_act = numbers of time energy was used. 
                    #This is only needed because translation unit can perform multiple actions
                    if a_type == 'R':
                        num_act = self.cut_seq(param)
                    elif a_type == 'L':
                        num_act = self.ligate_seq(param)
                    elif a_type == 'E':
                        num_act = self.copy_transport_seq(param)
                    elif a_type == 'I':
                        num_act = self.transport_seq(param)
                    elif a_type == 'M':
                        num_act = self.mutate_seq(param)
                    elif a_type == 'O':
                        num_act = self.copy_seq(param)
                    elif a_type == 'D':
                        num_act = self.delete_seq(param)
                    elif a_type == 'V':
                        num_act = self.convert_subs(param)
                    elif a_type == 'U':
                        num_act = self.transport_by_tag(param)
                
            #If it is a read action select a 
            elif a_type == 'read':
                tries = 0
                #randomly select a translation unit
                while tries < len(self.tr_count):
                    tries += 1
                    unit = random.choice(self.tr_count.keys())
                    if self.tr_count[unit] == 0:
                        del self.tr_count[unit]
                    else:
                        energy_used = self.sample_probability(self.units[unit]['energy'])
                        self.tr_count[unit] -= 1
                        num_act = self.read_seq(param, energy_per_unit = energy_used)
                        break
                
            #deduct energy for every action performed if not unlimited energy
            self.change_energy(-energy_used * num_act)
            
            #perform changes in sub cell if partner cell is not top cell
            if self.transported['partner_cell'] != -1:
                tr = self.transported
                c_id = tr['partner_cell']
                if tr['seq_export'] != None:
                    self.cells[c_id].sequence += tr['seq_export']
                if tr['seq_import'] != None:
                    pos = tr['seq_import']
                    self.cells[c_id].sequence = self.cells[c_id].sequence[:pos[0]] + self.cells[c_id].sequence[pos[1]:]
                if tr['unit_export'] != None:
                    self.cells[c_id].add_unit(tr['unit_export'])
                if tr['unit_import'] != None:
                    del self.cells[c_id].units[tr['unit_import']]
                if tr['subs_export'] != None:
                    self.cells[c_id].add_subs(tr['subs_export'])
                if tr['subs_import'] != None:
                    del self.cells[c_id].substances[tr['subs_import']]
        
    def find_occurences(self, template, start, stop, target, maxi):
        '''finds sequence (only sequence data no enzyme target sites)
        finds maximal number of occurences, fills up to maximum'''
        
        target = target.replace('X', '[GATCX]')
        
        if target[len(target)-1] == ';':
            target += '([^0-9]|$)'
        
        target_seq = lambda x: self.pos_target(x, 'own')
        seq = self.sequence
        
        if template == 'prom':
            target = '(P[0-9]*\.?[0-9]*)' + target
        
        elif template == 'ext':
            seq = self.partner_sequence
            target_seq = lambda x: self.pos_target(x, 'ext')
        
        results = re.compile(target).finditer(seq[start:stop])
        
        occ = []
        for res in results:
            pos = res.start() + start
            if template == 'prom':
                pos += len(res.groups()[0])
                
            if not target_seq(pos):
                occ.append(pos)
        
        #Limit/extend the number of result until max_action is reached  
        choices = []
        if len(occ)>0:
            random.shuffle(occ)
            if maxi < len(occ):
                choices = occ[:maxi]
            else:
                #Add random elements until max is reached
                choices = occ
                missing = maxi - len(choices)
                while missing > 0:
                    if len(occ) < missing:
                        choices += occ
                    else:
                        choices += occ[:missing]
                    missing = maxi - len(choices)
            
        return choices
    
    def read_seq(self, start = 0, energy_per_unit = 0, adress = False):
        '''Parse a sequence into unit(s) and add to pool'''
        
        if adress != (self.id,) and adress != False:
            adress = list(adress)
            adress.pop(0)
            try:
                self.cells[adress[0]].read_seq(start = start, energy_per_unit = energy_per_unit, adress = tuple(adress))
                return True
            except KeyError:
                #cell isn't here anymore
                return False
        else:
            try:
                l = len(self.sequence)
                if l == 0:
                    return False
            except IndexError:
                return False
            
            try:   
                while self.sequence[start] == ';' and start != len(self.sequence) -1:
                    start += 1
            except IndexError:
                return False
            
            #Find stop sign 'S' or end of sequence ';'
            sequence = self.sequence[start:]
            
            try:
                stop = sequence.index('S') + start
                seq_end = sequence.index(';') + start
                if seq_end < stop:
                    stop = seq_end
            except ValueError:
                try:
                    stop = sequence.index(';') + start
                except ValueError:
                    stop = len(self.sequence)
                    
            # Check if any position between start and stop is protected,
            #if so this moves the stop up to that position
            
            for test in range(start, stop):
                if self.pos_protected(test):
                    stop = test
                    break
                    
            sequence = self.sequence[start:stop]
            #Create regular expression pattern
            #matching for unit description
            pattern = r'([NRLFMEIOQDZUVY])'               #unit type [0]
            pattern += r'([\-]*\d+[\.]?\d*)'            #max action parameter (can be floating point) [1]
            pattern += r'/?(\d*[\.]?\d*)'               #optional energy parameter after / (can be floating point) [2]
            pattern += r'\(?([a-z_,\*]*)\)?'              #optional tags in () [3]
            pattern += r'\[[GTACXP0-9]*\]'              #contained sequence in [] 
            pattern += r'([GTACX\!]*)'                  #optional target sequence [4]
            pattern += r'\(?([a-z_,\*]*)\)?'              #optional tags targets [5] (or starting substances for conversion unit V)
                                                        #this also serves as conditional substances for all except V or U
            pattern += r'\(?([a-z_,\*]*)\)?'              #optional tags resulting substances (V) [6]
                                                        #or conditional substances for U
            pattern += r'([0-9]+)'                      #life parameter [7]
        
            regex = re.compile(pattern)
            new_units = regex.finditer(sequence)
            
            #little helper function to parse the target content into a 
            #tuple of strings
            def parse_tags(tag):
                pattern = r'([a-z_\*]+)'
                found_tags = re.compile(pattern).finditer(tag)
                all_tags = []
                for t in found_tags:
                    all_tags.append(t.groups()[0])
                    
                return tuple(all_tags)
            
            units_translated = 0
            #Add new units to pool and respective lists
            for unit in new_units:
                #Check if enough energy is present to translate next unit
                if self.has_energy(total = True) - (energy_per_unit * (units_translated+1)) >= 0:
                    params = unit.groups()
                    
                    action = float(params[1])
                    life = int(params[7])
                    tags = parse_tags(params[3])
                    target = params[4]
                    u_type = params[0]
                    needs = parse_tags(params[5])
                    
                    #negative parameters only allowed for transport (U) and Transcription Factor (F)
                    if u_type not in ('F', 'U', 'I', 'E'):
                        action = abs(action)
                    
                    if len(params[2]) == 0:
                        energy = 0
                    else:
                        energy = float(params[2])
                    
                    if u_type in ('I', 'Q', 'D'):
                        try:
                            target = target.replace('!', '')
                        except:
                            pass
                    
                    elif u_type == 'U':
                        target = parse_tags(params[5])
                        needs = parse_tags(params[6])
                        
                    elif u_type == 'V':
                        target = parse_tags(params[6])
                        needs = parse_tags(params[5])
                        
                    unit = {'type': u_type, 'action': action, 'energy': energy, 'tags': tags, 'target': target, 'life': life}
                    
                    if len(needs) > 0:
                        unit['needs'] = needs 
                    
                    self.add_unit(unit)
                    
                    units_translated += 1
                    
                #if not enough energy is present return    
                else:
                    return units_translated
            return units_translated
        return False
    
    def sample_probability(self, value):
        '''This function converts a probability (max action parameter) into 
        the actual number of actions performed.
        E.g.     0.5 ==> 50% 0, 50% 1
                0.2 ==> 80% 0, 20% 1
                4.5 ==> 50% 4, 50% 5
                4.9 ==> 10% 4, 90% 5
                You can also use: 7.332, 0.11, 13.601 or even -1.224 for transcription factors
        '''
        if self.sample_mode == 'standard' or (-1 < value < 1):
            
            neg= False
            if value < 0:
                value = abs(value)
                neg = True
                
            new_value = int(value)
            
            #Check if value is floating point
            if new_value != value:
                rest = float(value - new_value)
                #select a random number
                resolution = 100000
                num = random.randrange(resolution)
                if num < rest * resolution:
                    new_value += 1
                    
            if neg:
                new_value = -new_value
            
        elif self.sample_mode == 'normal':
            
            neg = False
            if value < 0:
                neg = True
                value = abs(value)
                
            sigma = value * self.normal_factor
            new_value = round(random.normalvariate(value, sigma))
            
            if neg:
                new_value = -new_value
                
            if (new_value < 0 and not neg) or (new_value > 0 and neg):
                new_value = 0
                
        return int(new_value)
    
    def degrade_units(self):
        '''Delete all units with 0 live'''
        for unit in self.units.keys():
            if self.units[unit]['life'] <= 0:
                #delete unit entry...
                del self.units[unit]
                
        #degrade units in all sub cells
        for ind in self.cells.keys():
            self.cells[ind].degrade_units()
                
    def reset_target(self):
        '''Make a list of all target sequences (start, end). these sites can not be acted on.
        This also includes a list of targets for the current external sequence.'''
        
        pattern = r'\]([GTACX\!]*\(?[a-z_]*\)?)([0-9]+)'
        target_sequences = []
        targets = re.compile(pattern).finditer(self.sequence)
        for target in targets:
            start = target.start() + 1
            end = target.start() + len(target.groups()[0])
            target_sequences.append((start, end))
        
        self.target_sequences = target_sequences
        
        targets = re.compile(pattern).finditer(self.partner_sequence)
        target_sequences = []
        for target in targets:
            length = len(target.groups()[0])
            if length > 0:
                start = target.start() + 1
                end = target.start() + length
                target_sequences.append((start, end))
        
        self.partner_target_sequences = target_sequences
        
    def pos_target(self, position, template):
        flag = False
        
        if template == 'own':
            for start, end in self.target_sequences:
                if start <= position <= end:
                    flag = True
                    break
                
        elif template == 'ext':
            for start, end in self.partner_target_sequences:
                if start <= position <= end:
                    flag = True
                    break
            
        return flag
    
    def reset_protected(self):
        '''Reset protected sites (at beginning of every round)'''
        self.pr_pos = []
        for unit_id in self.units.keys():
            if self.units[unit_id]['type'] == 'Q':
                target = self.units[unit_id]['target']
                
                max_action = self.sample_probability(self.units[unit_id]['action'])
                end = len(target)
                
                found_pos = self.find_occurences('own', 0, len(self.sequence), target, max_action)
                found_pos = list(set(found_pos))
                for pos in found_pos:
                    new_pos = (pos, pos + end -1)
                    self.pr_pos.append(new_pos)
    
    def change_sequence_length(self, start, length_diff, end=-1):
        '''Moves translation start positions and protected sites according to length_diff'''
        if end == -1:
            end = len(self.sequence)
        
        #update protected
        for i in range(len(self.pr_pos)):
            pr = self.pr_pos[i]
            if start <= pr[0] <= end and start <= pr[1] <= end:
                self.pr_pos[i] = (pr[0] + length_diff, pr[1] + length_diff)
        #update read positions
       
    def pos_protected(self, position):
        '''Checks if an available protection unit inhibits action at position'''
        protected = False
        for start, end in self.pr_pos:
            if start <= position <= end:
                protected = True
                break    
        return protected
    
    def disintegrate_sub_cell(self, adress, energy_used = 0):
        if adress[0] == self.id and len(adress) == 2:
            old_id = adress[1]
            #Check if enough energy is in the cell and it still exists
            try:
                success = self.cells[old_id].change_energy(-energy_used)
            except KeyError:
                return
            if success:
                self.sequence += self.cells[old_id].sequence
                #for every substance & unit and subcell in old cell 
                #transfer to new with a probability decided by asymmetry factor
                for unit in self.cells[old_id].units.values():
                    #add unit to new cell
                    self.add_unit(unit_dict = unit)
                
                #add all substances   
                subs = tuple(self.cells[old_id].substances.values())
                self.add_subs(subs)
                        
                for sub_id in self.cells[old_id].cells.keys():
                    #create subcell in mother_cell
                    seq = self.cells[old_id].cells[sub_id].sequence
                    name = self.cells[old_id].cells[sub_id].name
                    units = self.cells[old_id].cells[sub_id].units.values()
                    subs = self.cells[old_id].cells[sub_id].substances.values()
                    cells = self.cells[old_id].cells[sub_id].cells.values()
                    
                    self.new_sub_cell(sequence = seq, name = name, units = units, substances = subs, cells = cells)
                    
                #delete old cell
                del self.cells[old_id]
                
                return True
                
        else:
            adress = list(adress)
            adress.pop(0)
            try:
                success = self.cells[adress[0]].disintegrate_sub_cell(tuple(adress), energy_used)
                return success
            except KeyError:
                #cell isn't here anymore
                return False
        
        return False
    
    def split_sub_cell(self, adress, energy_used = 0, asymetry_factor = 0.5):
        if adress[0] == self.id and len(adress) == 2:
            old_id = adress[1]
            #Check if enough energy is in cell and if present
            try:
                success = self.cells[old_id].change_energy(-energy_used)
            except KeyError:
                #cell could not be found (got moved) ==> no action
                return False
            
            if success:
                #create new sub cell with identical sequence
                seq = self.cells[old_id].sequence
                name = 'Sub: ' + self.cells[old_id].name
                units = []
                subs = []
                cells = []
                
                #for every substance & unit and subcell in old cell 
                #transfer to new with a probability decided by asymmetry factor
                for unit_id, unit in self.cells[old_id].units.items():
                    if random.random() > asymetry_factor:
                        #add unit to new cell
                        units.append(unit)
                        #Remove unit from old cell
                        del self.cells[old_id].units[unit_id]
                
                for subs_id, subst in self.cells[old_id].substances.items():
                    if random.random() > asymetry_factor:
                        #add substance to new cell
                        subs.append(subst)
                        #Remove substance from old cell
                        del self.cells[old_id].substances[subs_id]
                
                for sub_id, cell in self.cells[old_id].cells.items():
                    if random.random() > asymetry_factor:
                        cells.append(cell)
                        #Remove subcell
                        del self.cells[old_id].cells[sub_id]
                
                new_id = self.new_sub_cell(sequence = seq, name = name, substances = subs, units = units, cells = cells)
                return new_id
            
        else:
            adress = list(adress)
            adress.pop(0)
            try:
                id = self.cells[adress[0]].split_sub_cell(tuple(adress), energy_used)
                return id
            except KeyError:
                #cell isn't here anymore
                return False
            
    def cut_seq(self, unit):
        '''Cut sequence once.'''
        
        target = self.units[unit]['target']
        #Find first '!' and remove
        pos1 = target.index('!')
        try:
            pos2 = pos1 + target[pos1+1:].index('!')
        except ValueError:
            pos2 = -1
        
        target = target.replace('!', '')
        tot = 0
        found = False
        #Again very inelegant...
        black_list = []
        while tot < 25 and not found:
            tot += 1
            occ = self.find_occurences('own', 0, len(self.sequence), target, 1)
            if len(occ) > 0:
                oc = occ[0]
                if oc not in black_list:
                    if not self.pos_protected(oc):
                        #Blunt cut site
                        if pos2 == -1:
                            self.sequence = self.sequence[:oc+pos1] + ';;' + self.sequence[oc+pos1:]
                            self.change_sequence_length(oc+pos1, 2)
                            found = True
                        #overlapping end cut site
                        else:
                            overlap = pos2-pos1
                            
                            forward = self.sequence[oc+pos1:]
                            backward = self.sequence[:oc+pos2]
                            
                            #check if end does not create too much overlap
                            target_forward = '(;)([0-9]*)'
                            #revert backward sequence to check for next match
                            target_backward = '([0-9]*)(;)'
                            back_rev = backward[::-1]
                            
                            #Perform regex search
                            res_forw = re.compile(target_forward).search(forward ,1)
                            res_backw = re.compile(target_backward).search(back_rev ,1)
                            
                            len_seq1 = res_forw.start()
                            len_seq2 = res_backw.start()
                            
                            try:
                                overlap_tot1 = overlap + int(res_forw.groups()[1])
                            except ValueError:
                                overlap_tot1 = overlap
                            try:
                                overlap_tot2 = overlap + int(res_backw.groups()[0])
                            except ValueError:
                                overlap_tot2 = overlap
                                
                            #Not too much overlapp
                            if len_seq1 > overlap_tot1 and len_seq2 > overlap_tot2:
                                overlap = str(overlap)
                                #Create overlapping ends 
                                self.sequence = backward + ';' + overlap + ';' + overlap + forward
                                #length of introduced sequence:
                                len_diff = 2 + 2 * len(overlap)
                                position = len(backward)
                                self.change_sequence_length(position, len_diff)
                                found = True
                    else:
                        black_list.append(oc)
        if found:                
            return 1
        else:
            return 0
                     
    def mutate_seq(self, unit):
        '''Introduce mutations; select random position
        choose one of following scenarios randomly:
        1. Delete position
        2. Insert random letter:
            Same algorithm used as change position
        3. Change position:
            50% chance: Normal sequence (G, A, T, C or X)
            40% chance: Units, P, S or digit
            10% chance: ;;, [random_letter], (random_tag_letter)
        '''
        letters = ['G','A','T','C','X']
        units = ['N','R','F','L','M','Q','E','I','O','D','Z','U','V','P','S','Y']
        for i in range(10):
            units.append(str(i))
        
        tags = []
        tag_string = 'abcdefghijklmnopqrstuvwxyz_'
        for letter in tag_string:
            tags.append(letter)
            
        #Get a position which is not protected
        i = 0
        found = False
        while i < 10 and not found:
            i += 1
            try:
                pos = random.randrange(len(self.sequence))
            except ValueError:
                return 0
            if not self.pos_protected(pos) and pos != 0 and pos != (len(self.sequence) -1):
                found = True
                break
        
        if found:
            sit = random.random()
            #delete position
            if sit < 0.3333:
                self.sequence = self.sequence[:pos] + self.sequence[pos+1:]
                self.change_sequence_length(pos, -1)
                return 1
            
            #choose random letter
            else:
                letter_type = random.random()
                #its a normal sequence letter (randomly chosen)
                if letter_type < 0.5:
                    new_letter = random.choice(letters)
                #its a unit, P, S or a digit
                elif 0.5 < letter_type < 0.9:
                    new_letter = random.choice(units)
                #its a special character
                else:
                    #a strand break
                    if letter_type < 0.93333:
                        new_letter = ';;'
                    #a contained sequence
                    elif 0.93333 < letter_type < 0.96666:
                        new_letter = '['+ random.choice(letters) + ']'
                    #a tag
                    else:
                        new_letter = '(' + random.choice(tags) + ')'
                
                #insert a random letter
                if sit < 0.6666:
                    self.sequence = self.sequence[:pos] + new_letter + self.sequence[pos:]
                #change a letter
                elif 0.6666 < sit:
                    self.sequence = self.sequence[:pos] + new_letter + self.sequence[pos+1:]
                    
                self.change_sequence_length(pos, len(new_letter))
                
                return 1
            
        return 0
        
    def ligate_seq(self, unit):
        target = self.units[unit]['target']
        
        try:
            ex_1 = target.index('!')
        except IndexError:
            return 0
        
        try:
            #its a overlapping target
            ex_2 = target[ex_1+1:].index('!') + ex_1
            overlap = ex_2 - ex_1
            target = target.replace('!', '')
            target2 = '(' + target[:ex_2] + ';' + str(overlap) + ')'
            target1 = '(;' + str(overlap) + ')' + target[ex_1:]
            len_before = ex_1
            
        except ValueError:
            #its a blunt end target
            overlap = 0
            target2 = '(' + target[:ex_1] + ';)'
            target1 = '(;' + target[ex_1+1:] + ')'
        
        reg1 = re.compile(target1.replace('X', '[GATCX]')).finditer(self.sequence)
        reg2 = re.compile(target2.replace('X', '[GATCX]')).finditer(self.sequence)
        
        occ1 = [(m.start(), m.groups()[0]) for m in reg1]
        occ2 = [(m.start(), m.groups()[0]) for m in reg2]
        
        if len(occ1) > 0 and len(occ2) > 0:
            count = 0
            while count < 10:
                count += 1
                sel1 = random.choice(occ1)
                sel2 = random.choice(occ2)
                
                if overlap == 0:
                    start1 = sel1[0]
                    res = re.compile('(;[0-9]*)').search(self.sequence[start1:], 1)
                    try:
                        end1 = start1 + res.start() + len(res.groups()[0]) - 1
                    except AttributeError:
                        return 0
                    
                    end2 = sel2[0] + len(sel2[1]) - 1
                    seq_rev = self.sequence[0:end2]
                    seq_rev = seq_rev[::-1]
                    start2 = end2 - seq_rev.index(';') - 1
                else:
                    start1 = sel1[0] + len(sel1[1]) - 1
                    res = re.compile('(;[0-9]*)').search(self.sequence[start1+3:], 1)
                    try:
                        end1 = start1 + res.start() + len(res.groups()[0]) + 3
                    except AttributeError:
                        return 0
                    
                    end2 = sel2[0]
                    seq_rev = self.sequence[0:end2-3]
                    seq_rev = seq_rev[::-1]
                    start2 = end2 - seq_rev.index(';') - 4
                    
                #Check if ligation does not happen on same strand (no circularization supported)
                if abs(start1 - start2) > 2:
                    if start1 > start2:
                        b = self.sequence[start1+1:]
                        
                        if overlap == 0:
                            a = self.sequence[:end2]
                            c = self.sequence[end2+1:start1]
                        else:
                            a = self.sequence[:end2+len_before]
                            c = self.sequence[end2+len(sel2[1]):start1-1]
                            
                        self.sequence = a + b + c
                        
                    else:
                        if overlap == 0:
                            a = self.sequence[:start1]
                            b = self.sequence[end1+1:end2]
                            c = self.sequence[start1+1:end1+1]
                            d = self.sequence[end2+1:]
                        else:
                            a = self.sequence[:start1-1]
                            b = self.sequence[end1:end2+len_before]
                            c = self.sequence[start1+1:end1]
                            d = self.sequence[end2+len(sel2[1]):]
                        
                        self.sequence = a + b + c + d 
                        
                    return 1
        return 0
         
    def convert_subs(self, unit):
        #try to find and delete the starting substances
        #if successful add resulting sequences
        try:
            start_subs = self.units[unit]['needs']
        except KeyError, e:
            start_subs = tuple()
            
        res_subs = self.units[unit]['target']
        if start_subs != tuple():
            success = self.find_tags(start_subs, remove = True)
            if success:
                self.add_subs(res_subs)
                return 1
            else:
                return 0
        else:
            self.add_subs(res_subs)
            return 1
    
    def transport_by_tag(self, unit):
        '''if there are any substances to export they always get preference over units'''
        all_targets = list(self.units[unit]['target'])
        random.shuffle(all_targets)
        
        #Exporting a substance or unit by tag
        if self.units[unit]['action'] > 0:
            for target in all_targets:
                success = self.find_tags((target,), units = True, remove = True)
                if success:
                    success = success[0]
                    if type(success[1]) == str:
                        self.transported['subs_export'] = success[1]
                    elif type(success[1]) == dict:
                        self.transported['unit_export'] = success[1]
                    return 1
        
        #Importing a substance or unit        
        elif self.units[unit]['action'] < 0:
            for target in all_targets:
                success = self.find_tags((target,), units = True, use_partner = True)
                if success:
                    success = success[0]
                    if type(success[1]) == str:
                        self.add_subs(success[1])
                        self.transported['subs_import'] = success[0]
                    elif type(success[1]) == dict:
                        self.transported['unit_import'] = success[0]
                        self.add_unit(success[1])
                    return 1
        return 0
    
    def copy_transport_seq(self, unit):
        '''Does all the actions of one export enzyme. Copy sequence and label for export.
        Caution: The export unit creates only double-stranded sequences even from 
        a overhanging template (overhang will be copied).
        Exported sequences will therefore always be blunt ended.
        '''
        #The target sequence without !
        target = self.units[unit]['target']
        
        if self.units[unit]['action'] > 0:
            sequence = self.sequence
            mode = 'own'
            length = len(self.sequence)
            
        elif self.units[unit]['action'] < 0:
            sequence = self.partner_sequence
            mode = 'ext'
            length = len(self.partner_sequence)
        
        pos = target.index('!')
        target = target.replace('!', '')
        i = 0
        while i < 15:
            i += 1
            occ = self.find_occurences(mode, 0, length, target, 1)
            if len(occ) > 0:
                occ = occ[0]
                occ += pos
                if mode == 'own':
                    if self.pos_protected(occ):
                        return 0
                
                #search the next end
                pattern = '(;[0-9]*)'
                seq = sequence[occ:]
                result = re.compile(pattern).search(seq, 1)
                    
                if result:
                    end_pos = result.start() + len(result.groups(0))
                    if mode == 'own':            
                        self.transported['seq_export'] = ';' + seq[:end_pos]
                        return 1
                    if mode == 'ext':
                        self.sequence += ';' + seq[:end_pos]
                        return 1
        return 0
    
    def copy_seq(self, unit):
        '''Copy sequence and attach to end of sequence.
        Caution: The copy unit creates only double-stranded sequences even from 
        a overhanging template (overhang will be copied). 
        Copied sequences will therefore always be blunt ended.
        !!! Possible FIX: Nearly identical function to export unit.
        '''
        #The target sequence without !
        target = self.units[unit]['target']
        pos = target.index('!')
        target = target.replace('!', '')
        count = 0
        black_list = []
        while count < 25:
            count += 1
            occ = self.find_occurences('own', 0, len(self.sequence), target, 1)
            if len(occ) > 0:
                start = occ[0] + pos
                if start not in black_list:
                    if not self.pos_protected(start):
                        #search the next end
                        pattern = '(;[0-9]*)'
                        seq = self.sequence[start:]
                        result = re.compile(pattern).search(seq, 1)
                        end_pos = result.start() + len(result.groups(0))
                        
                        #Add copied string to sequence
                        self.sequence += ';' + seq[:end_pos]
                        return 1
                    else:
                        black_list.append(start)
        return 0 
            
    def transport_seq(self, unit):
        '''Makes exactly one transport'''
        #The target sequence without !
        target = self.units[unit]['target'].replace('!', '')
        count = 0
        if self.units[unit]['action'] < 0:
            sequence = self.partner_sequence
            mode = 'ext'
            length = len(self.partner_sequence)
            
        elif self.units[unit]['action'] > 0:
            sequence = self.sequence
            mode = 'own'
            length = len(self.sequence)
        
        else:
            return 0
            
        while count < 25:
            count += 1
            occ = self.find_occurences(mode, 0, length, target, 1)
            if len(occ) > 0:
                occ = occ[0]
                #search the next forward end
                target_forward = '(;[0-9]*)'
                #revert backward sequence to check for next match
                back_rev = sequence[:occ]
                back_rev = back_rev[::-1]
                #search results in forward direction
                res_forw = re.compile(target_forward).search(sequence[occ:], 1)
                #search results in backward direction with inverted positions
                start_pos = occ - back_rev.index(';') - 1
                if not res_forw:
                    return 0
                
                end_pos = occ + res_forw.start() + len(res_forw.groups()[0])
                new_sequence = sequence[start_pos : end_pos]
                old_seq = sequence[:start_pos] + sequence[end_pos:]
                if old_seq:
                    if old_seq[len(old_seq)-1] != ';':
                        end_pos += 1
                
                if mode == 'ext':
                    #add new sequence to own sequence
                    self.sequence += new_sequence
                            
                    #remove from top cell
                    self.transported['seq_import'] = (start_pos, end_pos,)
                    return 1
                
                elif mode == 'own':
                    #remove sequence from own sequence
                    self.sequence = old_seq
                
                    #export in partner cell
                    self.transported['seq_export'] = new_sequence
        return 0
                
    def delete_seq(self, unit):
        '''Does all the actions of one sequence degradation unit.
        Very similar too import unit'''
        
        #The target sequence
        target = self.units[unit]['target']
        count, found = 0, False
        
        while count < 3 and not found:
            count += 1
            occ = self.find_occurences('own', 0, len(self.sequence), target, 1)
            
            if len(occ) > 0:
                occ = occ[0]
                #Check if any letter in the target is protected
                protected = False
                for test in range(occ, occ+len(target)):
                    if self.pos_protected(test):
                        protected = True
                        break
                    
                if not protected:
                    found = True
                    #search the next forward end
                    target_forward = '(;[0-9]*)'
                    
                    #revert backward sequence to check for the beginning of the strand
                    back = self.sequence[:occ]
                    start_pos = occ - back[::-1].index(';') - 1
                    
                    #search results in forward direction
                    res_forw = re.compile(target_forward).search(self.sequence[occ:], 1)
                    end_pos = occ + res_forw.start() + len(res_forw.groups()[0])
                    #remove sequence
                    before = self.sequence[:start_pos]
                    after = self.sequence[end_pos:]
                    if len(before) > 0 and before[len(before)-1] != ';':
                        before += ';'
                    if len(after) > 0 and after[0] != ';':
                        after = ';' + after
                    self.sequence = before + after
                    
                    #remove all the reads and protected sites that were on the strand
                    to_del = []
                    for i in range(len(self.pr_pos)):
                        start, end = self.pr_pos[i]
                        if start_pos < start < end_pos or start_pos < end < end_pos:
                            to_del.append(i)
                    
                    for i in to_del:
                        del self.pr_pos[i]
                     
                    #!!! update positions   
                    #for action_id in self.round_actions.keys():
                    #    action = self.round_actions[action_id]
                    #    if action[0] == 'read' and start_pos < action[1] < end_pos:
                    #        del self.round_actions[action_id]
                            
                    #update length
                    self.change_sequence_length(start_pos, end_pos - start_pos)
                    
                    return 1
                
        return 0
              
    def parse_status(self):
        '''Parse units into string'''
        
        if self.name != '':
            resp_str = 'Name: ' + self.name + '\n'
        else:
            resp_str = 'Unnamed Cell\n'
            
        resp_str += 'ID: {0}, Round: {1}\n\n'.format(self.id, self.rounds)
            
        resp_str += 'Sample mode: ' + str(self.sample_mode)
        if self.sample_mode == 'normal':
            resp_str += ', Normal factor: ' + str(self.normal_factor)
        
        resp_str += '\nAllow transport with subcells: {0}\n'.format(self.interact_with_subcells)
        
        if self.subs_const != {}:
            resp_str += 'Constantly added substances:\n'
            for sub, amount in self.subs_const.items():
                resp_str += '    ' + sub + ': ' + str(amount) + ' per round\n'
        
        if len(self.substances) > 0:   
            resp_str += '\nSubstances:\n'
            
            subs = [value for key, value in self.substances.items()]
            subs = Counter(subs)
            for sub in subs.keys():
                resp_str += sub + ': ' + str(subs[sub]) + 'x\n'
        
        if len(self.units) > 0:
            resp_str += '\nUnits:\n'
            
            for ind in self.units.keys():
                resp_str += '- ID: ' + str(ind) + ', '
                resp_str += 'Type: ' + self.units[ind]['type']
                if self.units[ind]['tags'] != '':
                    resp_str += ', Tags: '
                    for tag in self.units[ind]['tags']:
                        resp_str += tag + ', '
                resp_str += 'Life: ' + str(self.units[ind]['life'])
                resp_str += ', Max: ' + str(self.units[ind]['action'])
                if self.units[ind]['energy'] != 0:
                    resp_str += ', Energy/Action: ' + str(self.units[ind]['energy'])
                    
                if self.units[ind]['type'] != 'V':
                    if self.units[ind]['target'] != '':
                        if type(self.units[ind]['target']) == tuple:
                            if  self.units[ind]['target'][0] != '':
                                resp_str += ', Target: '
                                for targ in self.units[ind]['target']:
                                    resp_str += ', ' + targ
                            else:
                                pass
                        else:
                            resp_str += ', Target: ' + self.units[ind]['target']
                else:
                    resp_str += ', From: '
                    for target in self.units[ind]['target'][0]:
                        resp_str += target + ', '   
                    resp_str += 'To: '
                    for target in self.units[ind]['target'][1]:
                        resp_str += target + ', '
                    
                resp_str += '\n'
        
        if len(self.sequence) > 0:
            resp_str += '\nSequence:\n'
            resp_str += self.sequence + '\n'
        
        return resp_str
    
    def write_log(self, cont = False):
        #write log in cell and all sub cells
        #cast lists of unit and substance values (otherwise it will be dynamic
        #and the values would always reflect the same round...
        
        units = list(self.units.values())
        substances = list(self.substances.values())
        self.log[self.rounds] = {'units': units, 'energy': self.has_energy(total = True), 'sequence': self.sequence, 'substances': substances}
        
        if cont:
            for cell in self.cells.keys():
                self.cells[cell].write_log(cont = True)
                
    def save_cell_log_csv(self, filename, adress = False):
        '''Adress has to contain the current cell id as first parameter'''
        
        if adress != False:
            #redirect call to subcell 
            cell = self.get_subcell_by_adress(adress)
            cell.save_cell_log_csv(filename, adress = False)
        else:
            f = open(filename, 'wt')
            try:
                writer = csv.writer(f, delimiter = '\t')
                #set titles of columns
                writer.writerow( ('Round', 'Length Sequence', 'Energy', '#Substances','#F', '#N','#R','#L','#M','#E','#I','#O','#Q','#D','#Z','#U','#V','#Y','Sequence'))
                
                for rnd, entry in self.log.items():
                    #Count number of ocurences of each unit type
                    u = Counter([unit['type'] for unit in entry['units']])
                    
                    #write data in column
                    writer.writerow((str(rnd), len(entry['sequence']), entry['energy'], len(entry['substances']), 
                                         u['F'], u['N'], u['R'], u['L'], u['M'], u['E'], 
                                         u['I'], u['O'], u['Q'], u['D'], u['Z'], u['U'], u['V'], u['Y'], entry['sequence']))
                    
            finally:
                f.close()