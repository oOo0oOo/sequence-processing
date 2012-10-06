import unittest, operator
import time
import cell

class Benchmarks(unittest.TestCase):
    
    def test_find_occurence(self):
        '''Test how long a basic function takes.'''
        
        tries = 1000
        print 'Time per operation: find_occurences (average over {0} tries):\n'.format(tries)
        
        sequences = [';' + 'AG' * 10 + ';',
                     ';' + 'AG' * 500 + ';']
        
        for sequence in sequences:
            c = cell.Cell(sequence)
            before = time.time()
            for i in range(tries):
                    c.find_occurences('own', 0, len(sequence), 'AGA', 1000)
            time_passed = 1000 * (time.time() - before)/tries
            print 'Normal: Length sequence: {0}, Time/Action: {1} ms'.format(len(sequence), time_passed)
            
        prom_sequences = [';' + 'PAGAC' * 5 + ';',
                          ';' + 'PAGAC' * 190 + ';']
        
        for sequence in prom_sequences:
            c = cell.Cell(sequence)
            before = time.time()
            for i in range(tries):
                    c.find_occurences('prom', 0, len(sequence), 'AGA', 1000)
            time_passed = 1000 * (time.time() - before)/tries
            print 'Promoter: Length sequence: {0}, Time/Action: {1} ms'.format(len(sequence), time_passed)
            
        time.sleep(0.1)
        
    def test_one_operation(self):
        '''Test how long one operation takes. Subtract baseline of cell creation (empty sequence)...'''
        sequences = [
                    ('Baseline',            ';XXXXXXXXXXXXXXXXXXXXXX;'), #baseline: cell and unit creation, one round. no action
                    ('Protection',          ';Q1[]TT1SXXXXXXXXTTXXXX;'), #one protection
                    ('Promoter',            ';N1[]1SP1N1[]1XXXXXXXXX;'), #one promoter
                    ('Restriction',         ';R1[]T!T1SXXXXXXXXXTTXX;'), #one cut
                    ('Ligation',            ';L1[]T!T1SXXXXXT;;TXXXX;'), #one ligation
                    ('Mutation',            ';M1[]1SXXXXXXXXXXXXXXXX;'), #one mutation
                    ('Copy',                ';O1[]T!T1SXXXXXXXTTXXXX;'), #one copy
                    ('Degradation',         ';D1[]TT1XXXXXXXXXTTXXXX;'), #one degradation
                    ('Factor',              ';N1[]1F1[]T!T1SPTTXXXXX;'), #one transcription factor
                    ('Export',              ';E1[]T!T1SXXXTTXXXXXXXX;'), #one export
                    ('Four Units',          ';N0[]1N0[]1N0[]1N0[]1XX;'), #a unit exporting itself into the pool
                    ('Unit Export',         ';U1(self)[](self)2XXXXX;'),
                    ('Substance Creation',  ';V1[]()(a)1XXXXXXXXXXXX;'),
                    ('Multi Creation',      ';V1[]()(a,b,c,d,e,f,g)1;'),
                    ('3x Promoter',         ';N3[]1SP3N1[]1XXXXXXXXX;')
                    ]
        
        num_cells = 2500
        
        print '\nTime per operation of a specific unit (average over {0} cells):\n'.format(num_cells)
        
        res_dict = {}
        for label, seq in sequences:
            before = time.time()
            for i in range(num_cells):
                c = cell.Cell(seq)
                c.read_seq()
                c.run_all()
                
            res_dict[label] = ((time.time() - before), seq)
            
        baseline = res_dict['Baseline'][0]
        
        for test in res_dict.keys():
            per_cell = 1000 * (res_dict[test][0] - baseline) / num_cells
            res_dict[test] = ( per_cell , res_dict[test][1] )
        
        res = {}
        for k in res_dict.keys():
            res[k] = res_dict[k][0]
        
        sorted_res = sorted(res.iteritems(), key=operator.itemgetter(1))
        
        print 'Baseline (subtracted from all other): {1} ms / cell\nSequence: {0}'.format(res_dict['Baseline'][1], baseline)
        print 'All tests ordered by time / action:\n'
        
        for i in range(1, len(sorted_res)):
            current = sorted_res[i]
            percentage = round(current[1]/(baseline / 100),1)
            position = '{0}. {1}: {2}\n'.format(i, current[0], res_dict[current[0]][1])
            position += '\t{0}% of baseline time, {1} ms / action'.format(percentage, round(current[1], 4))
            print position
        
        #wait, so that end message of test does not conflict with printing of last results
        time.sleep(0.1)
        
    def test_sub_cells(self):
        '''Run a set of sequences in different levels of subcells'''
        
        sequences = [
                    (';F-3[]!T2P1N10[]3SP3AF-3[]!T2SP3TF-3[]!C2SP3CF-3[]!A2S;', 'Triple Oscillator'),
                    (';' + 'G' * 1000000 + ';', 'Very long sequence (1000000 letters)'),
                    (';N10[]1P1N10[]2O1[]A!A2AAAAA;', 'Endless copy; sequence grows every round'),
                    (';N10[]1P1N10[]2R2[]A!A2L1[]A!A2S' + 'A' * 100 + ';', 'Endless restriction/ligation')
                     ]
        
        print '\nTest different layers of sub cells:\n'
        rounds = 500
        
        for seq, name in sequences:
            print name + ': '
            
            seq_res = {}
            
            #top cell
            c = cell.Cell(seq)
            c.read_seq()
            before = time.time()
            for i in range(rounds):
                c.run_all()
            seq_res[1] = time.time() - before
            
            #first layer
            c = cell.Cell()
            ind = c.new_sub_cell(seq)
            c.cells[ind].read_seq()
            before = time.time()
            for i in range(rounds):
                c.run_all()
            seq_res[2] = time.time() - before
            
            #Second layer
            c = cell.Cell()
            ind1 = c.new_sub_cell()
            ind2 = c.cells[ind1].new_sub_cell(seq)
            c.cells[ind1].cells[ind2].read_seq()
            before = time.time()
            for i in range(rounds):
                c.run_all()
            seq_res[3] = time.time() - before
            
            #Third layer
            c = cell.Cell()
            ind1 = c.new_sub_cell()
            ind2 = c.cells[ind1].new_sub_cell()
            ind3 = c.cells[ind1].cells[ind2].new_sub_cell(seq)
            c.cells[ind1].cells[ind2].cells[ind3].read_seq()
            before = time.time()
            for i in range(rounds):
                c.run_all()
            seq_res[4] = time.time() - before
            
            #Fourth layer
            c = cell.Cell()
            ind1 = c.new_sub_cell()
            ind2 = c.cells[ind1].new_sub_cell()
            ind3 = c.cells[ind1].cells[ind2].new_sub_cell()
            ind4 = c.cells[ind1].cells[ind2].cells[ind3].new_sub_cell(seq)
            c.cells[ind1].cells[ind2].cells[ind3].cells[ind4].read_seq()
            before = time.time()
            for i in range(rounds):
                c.run_all()
            seq_res[5] = time.time() - before
            
            #Fifth layer
            c = cell.Cell()
            ind1 = c.new_sub_cell()
            ind2 = c.cells[ind1].new_sub_cell()
            ind3 = c.cells[ind1].cells[ind2].new_sub_cell()
            ind4 = c.cells[ind1].cells[ind2].cells[ind3].new_sub_cell()
            ind5 = c.cells[ind1].cells[ind2].cells[ind3].cells[ind4].new_sub_cell(seq)
            c.cells[ind1].cells[ind2].cells[ind3].cells[ind4].cells[ind5].read_seq()
            before = time.time()
            for i in range(rounds):
                c.run_all()
            seq_res[6] = time.time() - before
            
            base_time = (float(seq_res[1])*1000/rounds)
            res_str = '\tIn top cell: {0} ms / round\n'.format(round(base_time, 4))
            
            for level, time_needed in seq_res.items():
                if level != 1:
                    this_time = float(time_needed)*1000/rounds
                    diff = abs(this_time - base_time)
                    perc = 100 * (diff/base_time)
                    if diff > 0:
                        state = 'slower'
                    else:
                        state = 'faster'
                    
                    res_str += '\tIn sub cell {0}: {1} ms ({3}%) {2}\n'.format(level-1, round(diff, 3), state, round(perc, 3))
            print res_str
        
        time.sleep(0.1)
if __name__ == '__main__':
    unittest.main()