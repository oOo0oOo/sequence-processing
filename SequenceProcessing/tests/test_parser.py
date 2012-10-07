import unittest
import parse
from parse import SequenceError, SequenceWarning

class TestSequenceWarnings(unittest.TestCase):
    
    def check_sequences(self, sequences):
        '''Checks sequence for SequenceWarnings during parsing and returns the
        total number of sequences triggering SequenceWarnings.'''
        
        test = 0
        for sequence in sequences:
            #print sequence
            try:
                parser = parse.Parser()
                parser.sequence = sequence
                parser.check_syntax()
            except SequenceWarning, e:
                #print e
                test += 1
        return test
    
    def test_sequence_warning_unit_targets(self):
        
        sequences = [';I3(a)[]G!C10;', ';I3(a)[]10;', ';D3(a)[]4;', ';Q3(a)[]G!C10;', ';Q3(a)[]10;', #only target without !
                     ';M1(a)[]CCCT2;',';M1(a)[]X!X2;', ';N3(a)[]XCXC54;',';N3(a)[]XC!XC54;', #no target
                     ';Z1(a)[]CCCT2;',';Z1(a)[]X!X2;', ';Y1(a)[]CCCT2;',';Y1(a)[]X!X2;',
                     ';F3(a)[]4;', ';R3(a)[]4;', ';L3(a)[]4;', ';U3(a)[]4;', ';E3(a)[]4;' #normal target
                     ';F3(a)[]XX4;', ';R3(a)[]XX4;', ';L3(a)[]XX4;', ';E-3(a)[]TTC4;']
        
        res = self.check_sequences(sequences)         
        self.assertEqual(res, len(sequences), 'Sequence failed to trigger SequenceWarning.')

class TestParserSequenceInfo(unittest.TestCase):
    
    def test_sequence_info(self):
        #;N10[]1F1[]G!AT1SXXXP3GATXXXR2[]CCC!TT!CCC3;
        sequences = [
                     ';N10[]1;', ';P3AAAGGGM10[TTC]3AASTTR1[XXC]X!C3TT;',
                     ';N2[GTAXXXTAAAXXXXXXXXA]10PATTACAF5[GATACA]GA!XTACA5SPGATTACAPGAXTACAR5[GAATACA]GA!XT!ACA4SGATTACAP3GATTACAGAGGATTACAGP45ATTACAGATTACA;'
                     ]
        for sequence in sequences:
            parser = parse.Parser()
            parser.sequence = sequence
            ret = parser.sequence_info()
            #print '\nSequence info:\n', ret
        
class TestParserCheckSequence(unittest.TestCase):
    
    def check_sequences(self, sequences):
        '''Checks sequence for SequenceErrors during parsing and returns the
        total number of SequenceErrors produced by all sequences.'''
        
        test = 0
        for sequence in sequences:
            #print 'test_sequence: ', sequence
            try:
                parser = parse.Parser()
                parser.sequence = sequence
                parser.check_syntax()
                #print sequence
            except SequenceError, e:
                #print e, sequence
                test += 1
        return test   
    
    def test_valid_sequence(self):
        sequences = [
                     ';AGTCT;',';3GTACP3PC;2;CCPTGA;1',
                     ';F-3[]ATT!G5;', ';N1[]3SF-1[]ATT!G5;',
                     ';R0.733[]A!TTG3;',';L0.42[]AT!T!G3;',
                     ';P0.2R0.733[]A!TTG3;', ';N2.11[]3;',
                     ';Z0.733[]3;', ';N2.11[]3Z1.313[]3S;',
                     ';Q1.234[]ATTG3;', ';M5.432[]3;',
                     ';N2.11[](a)3Z1.313[](a)3S;',
                     ';F-3[]ATT!G()5;', ';N1[](a)3SF-1[]ATT!G()5;',
                     ';R0.733[]A!TTG(a)3;',';L0.42[]AT!T!G(a)3;',
                     ';Q1.234[]ATTG(a)3;', ';M5.432[](a)3;',
                     ';P5N3[]7SR15[]GATTA!CA5GGX;5;XXTAC;',
                     ';D3[]ATTG5;', ';D1[]XXCAC1;;CCCA;;GTCCCA;',
                     ';P1N5[]15R15[]GATTA!CA5F5[]GATT!ACA15L15[]GAT!TACA5;',
                     ';L10[]!CC!A3SAAAG;2;2CCAG;;3TTTTT;GGCC;2;TTGA;1',
                     ';R10[]A!GG!3AAAAGG;2;2GGGG;', ';AGGTC;3;TTTCCCCCCCCCCCCCCCCCC;12',
                     ';Q10[]AAATTTG3;', ';Q1[]X1;',
                     ';P21.33N10.3[]1F-76.44[]AT!2;',
                     ';2AAR10[]A!GG!3AAAAGG;2;2GGGG;1',
                     ';Q3[]XACC3O3[]!XACC2;;AACC;2',
                     ';Q3/5[]XACC3O3/2.3[]!XACC2;;AM0.955/0.775[]1ACC;2',
                     ';50GGGGGGGTTTTTAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGAGGGGGGGGGGTTTTTTT;',
                     ';N1/1(translation_unit_one)[]1;',
                     ';N1/1(t)[]1M1.22/1.223(t)[]1M1/1()[]1;',
                     ';N1(translation_*nit)[]1R1(restriction_unit)[]A!T1M3.33/3.333(mutation_unit)[]3;',
                     ';U1[](translation)1U-1/6.29[](mutation_unit)2;',
                     ';V1[](glucose,glucose*)(mannose)3;',
                     ';Z1[](energy,*,**)2;',
                     ';Y1/5[]1;', ';V1[]()()1;', ';V1[]()1;'
                     ';I-5[]AAA10;', ';I-3[]TTC10;',
                     ';I5/3.5[]ACCT1;',
                     ';E-5[]AA!A10;', ';E-3[]T!TC10;',
                     ';E5/3.5[]A!CCT1;'
                     ]
        
        test = self.check_sequences(sequences)         
        self.assertEqual(test, 0, 'Valid sequence triggered SequenceError or SequenceWarning.')
        
       
    def test_semicolon_syntax(self):
        sequences = ['ATTGG;', ';ATTGG', ';ATTGGP3', ';ATGGC;3GTCA;' , ';AAA;;;AAA;',
                     ';;CCAAAA;', ';GTCCCA;;', ';AGCXX;;;3CCCC;']
        test = self.check_sequences(sequences)    
        self.assertEqual(test, len(sequences), 'Semicolon failed to trigger SequenceError...')
        
    def test_brackets_syntax(self):
        sequences = [';AT[[TGG;', ';AT[[TG]G;', ';N3AX[]5;']
        test = self.check_sequences(sequences)         
        self.assertEqual(test, len(sequences), 'Bracket failed to trigger SequenceError...')
    
    def test_forward_slash_syntax(self):
        sequences = [';AT//TGG;', ';AT/TGG;', ';N3/[]5;', ';N/33[]5;']
        test = self.check_sequences(sequences)         
        self.assertEqual(test, len(sequences), 'Forward slash failed to trigger SequenceError...')
    
    def test_point_syntax(self):
        sequences = [';AA.AA;', ';AATTTC;3.2', ';AAAAA;;1.33TTTCC;', ';AAAAA;1.33;TTTCC;',
                     ';AAP3..1AA;', ';R1[]AA!AA2.3;', ';M1[]0.25;', ';N1[]3.92;', ';I1[]TTCC0.12;',
                     ';TS0.1AA;']
        test = self.check_sequences(sequences)    
        self.assertEqual(test, len(sequences), 'Invalid point (.) syntax failed to trigger SequenceError...')
        
    def test_illegal_characters(self):
        
        chars = '%+&=`:<>#?${}\\@^BHJKW'
        count = 0
        for char in chars:
            sequence = ';XAAA' + char + 'GGG;'
            try:
                parser = parse.Parser()
                parser.sequence = sequence
                parser.check_syntax()
            except SequenceError:
                count += 1
        
        self.assertEqual(count, len(chars), 'Failed to detect illegal character.')
        
    def test_tag_syntax(self):
        sequences = [';(sdf)AAT;', ';F(sdf)1/1[]A!AT2;', ';(sdf)F1/1[]A!AT2;', ';F1/1[]A(sdf)!AT2;',
                     ';N1(sdf[]1AAT;',';N1)[]1AAT;', ';N1(a[]1N1b)[]1;', ';N1([]1N1)[]1;',
                     ';M1/1(a16a)[]1;',';M1/1(aAAa)[]1;',';M1/1(aj1_A_TTC_a)[]1;']
        test = self.check_sequences(sequences)
        self.assertEqual(test, len(sequences), 'Wrong tag syntax failed to trigger SequenceError...')
        
    def test_exclamation_syntax(self):
        sequences = [';AS!TGG;', ';AXTTCG;!', '!;AXTTCG;', ';CP3!4TC;',
                     ';R1[]GTC!!!CCT3;', ';N!10[]3;', ';N10[!AA]3;', ';N10[AA!]3;',
                     ';R!1[]GTC!CCT3;']
        test = self.check_sequences(sequences)
        self.assertEqual(test, len(sequences), 'Wrong ! failed to trigger SequenceError...')
        
    def test_illegal_overlap(self):
        sequences = [';3AAATG;3', ';3AAPTG;3', ';ATGA;5', ';6ATGAC;', ';16ATGACTTTTTTTTTT;']
        test = self.check_sequences(sequences)
        self.assertEqual(test, len(sequences), 'Illegal overlap failed to trigger SequenceError...')
        
    def test_parameter_syntax(self):
        sequences = [';AAPTS2G;', ';M10[3AA]1;', ';ATGP1234AC;', ';R-10.12[]AT!3;', ';M-3[]AT!3;', ';M3[]-3;',
                     #';F-1[]A!-2;', 
                     ';D-1.3[]A!G2;', ';D1[]A!G-22;']
        test = self.check_sequences(sequences)
        self.assertEqual(test, len(sequences), 'Illegal parameter failed to trigger SequenceError...')

if __name__ == '__main__':
    unittest.main()   