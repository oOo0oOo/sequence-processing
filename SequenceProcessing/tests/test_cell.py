import unittest
import cell
from collections import Counter

class TestCellGeneration(unittest.TestCase):
    
    def test_generation(self):
        c = cell.Cell()
        self.assertIsInstance(c, cell.Cell)
        self.assertEqual(c.sequence, '')
        self.assertEqual(c.substances, {})
        self.assertEqual(c.units, {})
        c = cell.Cell('')
        self.assertIsInstance(c, cell.Cell)
        self.assertEqual(c.units, {})
        c = cell.Cell(';N50[]2;')
        self.assertIsInstance(c, cell.Cell)
        self.assertEqual(c.units, {})
        
        for i in range(3):
            c.read_seq()
            c.run_all()
        
        c1 = cell.Cell(units = c.units.values())
        self.assertIsInstance(c1, cell.Cell)
        self.assertEqual(c1.units.values().sort(), c.units.values().sort())
        
        c1.new_sub_cell(units = c1.units.values())
        c1.new_sub_cell(units = c1.units.values(), cells = c1.cells.values())
        c1.new_sub_cell(units = c1.units.values(), cells = c1.cells.values())
        c1.new_sub_cell(units = c1.units.values(), cells = c1.cells.values())
        c1.new_sub_cell(units = c1.units.values(), cells = c1.cells.values())
        
        for ce in c1.cells.values():
            self.assertEqual(ce.units, c1.units)
            self.assertEqual(ce.unit_count, c1.unit_count)
            for subcell in ce.cells.values():
                self.assertEqual(subcell.units, c1.units)
                self.assertEqual(subcell.unit_count, c1.unit_count)
                
class TestGetSubcell(unittest.TestCase):
    def test_get_subcell(self):
        '''test 4th layers of sub cell'''
        
        c = cell.Cell(';AA;')
        c.new_sub_cell(';GG;')
        c.cells[0].new_sub_cell(';TT;')
        c.cells[0].cells[0].new_sub_cell(';CC;')
            
        c.cells[0].cells[0].cells[0].read_seq(0)
        
        adresses = [(0,0,0,0),(0,0,0),(0,0),(0,)]
        
        for adress in adresses:
            sub_cell = c.get_subcell_by_adress(adress)
            self.assertIsInstance(sub_cell, cell.Cell)

class TestSubCell(unittest.TestCase):
    def test_subcell_generation(self):
        '''test 4th layers of sub cell'''
        
        test_cases = [(';R10[]A!C1AAAACCCCCCC;', ';R10[]A!C1AAAA;;CCCCCCC;'),
                      (';E1[]G!A1I1[]AA2SGGGAAAA;', ';E1[]G!A1I1[]AA2SGGGAAAA;;AAAA;'),
                      (';E1[]G!A1I1[]TT2SGGGACCCCCCTT;', ';E1[]G!A1I1[]TT2SGGGACCCCCCTT;;ACCCCCCTT;'),
                      (';E1[]CC!1I1[]AA3GGCCTTAATT;', ';E1[]CC!1I1[]AA3GGCCTTAATT;;TTAATT;'),
                      (';E1[]G!A1I1[]CC1SGGGAAAA;',';E1[]G!A1I1[]CC1SGGGAAAA;'), #test bad roundtrips (expecting no change in sequence)
                      (';E1[]G!A1I1[]XG1SGGGAAAA;',';E1[]G!A1I1[]XG1SGGGAAAA;'),
                      (';E1[]C!A1I1[]AA1SGGGAAAA;',';E1[]C!A1I1[]AA1SGGGAAAA;')]
        
        for test in test_cases:
            
            c = cell.Cell()
            c.new_sub_cell()
            c.cells[0].new_sub_cell()
            c.cells[0].cells[0].new_sub_cell(test[0])
            
            c.cells[0].cells[0].cells[0].read_seq(0)
            
            c.run_all()
            c.run_all()
            self.assertEqual(c.cells[0].cells[0].cells[0].sequence, test[1])
            
    def test_cell_creation_with_parameters(self):
        tests = [';M10[]10;', ';R10[]A!A11L12[]A!A13Q14[]AA15;']
        
        for seq in tests:
            c1 = cell.Cell(seq)
            c1.read_seq()
            c2 = c1.new_sub_cell(seq, units = c1.units.values())
            self.assertDictEqual(c1.cells[c2].units, c1.units)
            self.assertEqual(c1.cells[c2].unit_count, c1.unit_count)
            
            c1.cells[c2].new_sub_cell(seq, units = c1.cells[c2].units.values())
            c2_1 = c1.new_sub_cell(seq, units = c1.units.values(), cells = c1.cells.values())
            
    def test_complex_actions(self):
        '''This test takes more than 95% of time...
            Its fine as long as it does not trigger an error.'''
        
        tests = [';AAAAI3[]AA3E3[]A!A3I5[]AA3Q3[]CC3;',
                 ';R10(re)[]A!C10U10[](re)10AAAACCCCCCC;', 
                 ';V3[]()(crazy_stuff)3;', ';N3[]4F-5.0[]!AAA4SP6AAAN1.0[]2S;',
                 ';Y1[]1R10[]A!A10L10[]A!A10Q10[]AA10P50SF-100[]20;', '', ';;']
        
        adresses = [(0,),(0,0,),(0,0,0,),(0,0,0,0,),(0,0,0,0,0,)]
        
        for seq in tests:
            #create cells
            c = cell.Cell(seq)
            c.new_sub_cell(seq)
            c.cells[0].new_sub_cell(seq)
            c.cells[0].cells[0].new_sub_cell(seq)
            c.cells[0].cells[0].cells[0].new_sub_cell(seq)
            
            #make transport into subcells possible
            c.change_parameters(cont = True, interact_with_subcells = 'True')
            
            for i in range(2):
                for adress in adresses:
                    c.read_seq(adress = adress)
                    c.read_seq(adress = adress)
                c.run_all()
            
            #split all cells
            for adr in adresses[1:]:
                c.split_sub_cell(adr)
                
            #make transport into subcells possible
            c.change_parameters(cont = True, interact_with_subcells = 'True')
                
            c.run_all()
            c.run_all()
            for adress in adresses[1:]:
                c.read_seq(adress = adress)
                c.disintegrate_sub_cell(adress)
            
            c.run_all()
            
            #make transport into subcells possible
            c.change_parameters(cont = True, interact_with_subcells = 'True')
            
            #and run 5 more times
            for i in range(3):
                c.run_all()
        
class TestUnitCreation(unittest.TestCase):
    
    def test_unit_creation(self): 
        test_cases = [#Units with energy and probabilistic parameters
                      #as well as conditional tags
                      (';M3.5/3[]1;',         [{'type': 'M', 'action': 3.5, 'energy': 3, 
                                                'target': '', 'life': 1, 'tags': tuple()}]),
                      (';Z1[](number_one)1;',         [{'type': 'Z', 'action': 1, 'energy': 0, 
                                            'target': '', 'life': 1, 'tags': tuple(), 'needs': ('number_one',)}]),
                      (';Z3.5/3[]1;',         [{'type': 'Z', 'action': 3.5, 'energy': 3, 
                                                'target': '', 'life': 1, 'tags': tuple()}]),
                      (';Q1/1.5[]TT(signal)1;',     [{'type': 'Q','action': 1,'energy': 1.5,
                                              'target': 'TT','life': 1,'tags': tuple(), 'needs': ('signal',)}]),
                      (';N10/0.33[]4;',     [{'type': 'N','action': 10,'energy': 0.33,
                                              'target': '','life': 4,'tags': tuple()}]),
                      
                      #Ligation, Restriction
                      (';L11/15.09[]X!X3;', [{'type': 'L','action': 11,'energy': 15.09,
                                              'target': 'X!X','life': 3,'tags': tuple()}]),
                      (';R11/0.0[]X!X(glucose)3;', [{'type': 'R','action': 11,'energy': 0,
                                            'target': 'X!X','life': 3,'tags': tuple(), 'needs': ('glucose',)}]),
                      
                      #Transcription Factor
                      (';F-0.733[]A!TTG3;', [{'type': 'F','action': -0.733,'energy': 0,
                                              'target': 'A!TTG','life': 3,'tags': tuple()}]),
                      (';F-10.23[]!ATT3;', [{'type': 'F','action': -10.23,'energy': 0,
                                             'target': '!ATT','life': 3,'tags': tuple()}]),
                      
                      #Sequence Degradation
                      (';D3[]AAAACCCTTT1;', [{'type': 'D','action': 3,'energy': 0,
                                              'target': 'AAAACCCTTT','life': 1,'tags': tuple()}]),
                      (';D10/3.9[]XXXCCC3;', [{'type': 'D','action': 10,'energy': 3.9,
                                               'target': 'XXXCCC','life': 3,'tags': tuple()}]),
                      (';D111[]XXX!CCC3;', [{'type': 'D','action': 111,'energy': 0,
                                             'target': 'XXXCCC','life': 3,'tags': tuple()}]),
                      
                      #Cell splitting
                      (';Y3.5/3[]1;',                [{'type': 'Y', 'action': 3.5, 'energy': 3, 
                                                       'target': '', 'life': 1, 'tags': tuple()}]),
                      (';Y3.5/3(cell_splitter)[]1;', [{'type': 'Y', 'action': 3.5, 'energy': 3, 
                                                       'target': '', 'life': 1, 'tags': ('cell_splitter',)}]),
                      
                      #Multiple Units
                      (';N10[]3F1[]A!A(stuff)2F-3[]A!A2SM0[]2;', [{'type': 'N','action': 10,'energy': 0,
                                                            'target': '','life': 3,'tags': tuple()}, 
                                                           {'type': 'F','action': 1,'energy': 0,
                                                            'target': 'A!A','life': 2,'tags': tuple(), 'needs': ('stuff',)}, 
                                                           {'type': 'F','action': -3,'energy': 0,
                                                            'target': 'A!A','life': 2,'tags': tuple()}]),
                      
                      (';N10[]4F-5.0[]!AAA4SP5AAAN1.0[]2S;', [{'type': 'N','action': 10,'energy': 0,
                                                               'target': '','life': 4,'tags': tuple()}, 
                                                              {'type': 'F','action': -5,'energy': 0,
                                                               'target': '!AAA','life': 4,'tags': tuple()}]),
                      
                      #tagged units
                      (';N1(a*,b,c,*d,e,f)[]1;', 
                       [{'type': 'N','action': 1,'energy': 0,'target': '',
                         'life': 1,'tags': ('a*', 'b', 'c', '*d', 'e', 'f')}]),
                      (';N1(transl*tion)[]1M1(mutation)[]1L1(ligation)[]A!A1;', 
                       [{'type': 'N','action': 1,'energy': 0,'target': '',
                         'life': 1,'tags': ('transl*tion',)}, 
                        {'type': 'M','action': 1,'energy': 0,'target': '',
                         'life': 1,'tags': ('mutation',)}, 
                        {'type': 'L','action': 1,'energy': 0,'target': 'A!A',
                         'life': 1,'tags': ('ligation',)}]),
                      
                      #units with tag as target
                      (';U1(bb)[](a*a)1U-2/1.1(cc)[](bb)2U3/2.2(a*a)[](cc)3;', 
                       [{'type': 'U','action': 1,'energy': 0,'target': ('a*a',),
                         'life': 1,'tags': ('bb',)}, 
                        {'type': 'U','action': -2,'energy': 1.1,'target': ('bb',),
                         'life': 2,'tags': ('cc',)}, 
                        {'type': 'U','action': 3,'energy': 2.2,'target': ('cc',),
                         'life': 3,'tags': ('a*a',)}]),
                      (';U-0.733[](**,b_t,a)3;', [{'type': 'U','action': -0.733,'energy': 0,
                                              'target': ('**', 'b_t', 'a'),'life': 3,'tags': tuple()}]),
                      
                      #conversion unit with multiple tags as target
                      (';V1/1(conversion_unit)[](substance_a)(substance_b)3;', 
                       [{'type': 'V','action': 1,'energy': 1,'target': ('substance_b',), 'needs': ('substance_a',),
                         'life': 3,'tags': ('conversion_unit',)}]),
                      (';V1/1(conversion_unit)[](substance_a)(substance_b,substance_c)3;', 
                       [{'type': 'V','action': 1,'energy': 1,'target': ('substance_b', 'substance_c'), 'needs': ('substance_a',),
                         'life': 3,'tags': ('conversion_unit',)}]),
                      (';V1/1(conversion_unit)[](substance_d,sub_c,b*,**)(aa,bb)3;', 
                       [{'type': 'V','action': 1,'energy': 1,'target': ('aa', 'bb'), 'needs': ('substance_d', 'sub_c' ,'b*', '**'),
                         'life': 3,'tags': ('conversion_unit',)}]),
                      (';V1[]()(a,b)3;', [{'type': 'V', 'action': 1,'energy': 0,'target': ('a', 'b'),
                                           'life': 3,'tags': tuple()}]),
                      (';V1[](a_b_c)()3;', [{'type': 'V','action': 1,'energy': 0,'target': tuple(), 'needs': ('a_b_c', ),
                                             'life': 3,'tags': tuple()}]),
                      (';V1/1.5[](a_b_c)5;', [{'type': 'V','action': 1, 'energy': 1.5, 'target': tuple(), 'needs': ('a_b_c', ), 
                                               'life': 5,'tags': tuple()}])
                      ]
        
        for test_case in test_cases:
            c = cell.Cell()
            c.new_sub_cell(test_case[0])
            c.cells[0].read_seq(0)
            units = c.cells[0].units.values()
            self.assertEqual(units, test_case[1])
            
    def test_unit_not_created(self):
        '''all of these sequences are generate using a strong mutation unit... they should not trigger an
        error but simply not translate any units'''
        
        test_cases = [#these sequences are generate using a strong mutation unit...
                      ';P05A0C[]DR1XXXXC(q)XQXX(5v)GATXAARXXXCXAXX([G]R)T[X)]Q93MGXXAXCXX;;SXT;;82TXUTXT;',
                      ';MZ0[CTL]AAPGX]1XGXXTXXX6CTXUXUQGAX7G75XZXGSAMA([T]c)PXMXX08XXXXV(b)8TX;',
                      ';E0[X1X(h)XXTTXXXXXXMXPSSQ(z)QQQ;',';TX[NCX]I0;;OX23T1XeTXAGAXAXX;5XXQXXXFXLC)TTTX5XXXL89C]XXXP(eTECXVXD3GD;',
                      ';ZMI1A[AAXXV6XXX5X(GtXNX2PFAC)0XCXCA8XGT811XXOGC(y)XXXGFFPGT9X]UGGGTRAXXXQXGXQX;',
                      ';SCGA0((l)i)0[]XGXI(uOXZXXXD3;7XTZG;;IXAVX;XXXGXX7AXXXXCD;',
                      ';X37[]0Z[P1X1XGXQATXXX8AG;XXXXX9n)1T[A];GXVCF(rMXG5XTXXGT;X8TX;;M(w)TX;',
                      ';US06X7G[C7XT4XXX(h)XXXDXX0FXANT[A[C3TPNPTd)RXXXGXAOXX657X[A];X8A]1XXMXSXXX;',
                      ';QTANCC;;F;;T0;VXXTD(i)XXXXZ;8A9;X5CAXXIAA8XGXXCX(c)CXXGAXCX(f)ECXSXXGXU8XTUNXXX;',
                      ';X]XGGXAINZXXGGXXj)XU[XCCSXXXXC[EXXNNU;XTXXX4X4XIIGT(;;XEXTRX;',
                      ';P(oC7CXXXCGRRCCA0;;;XXCUMXGG9LXXX6XXX;]XTOT[A]7XTXXCXXT)N;',
                      ';10TX1ATA(w;[2TXs;;XTX9XXX(p)XGXXCG(Oh8A)V;;AX(gTTGTX4AA]GGZXXXX2XASXAXX71(ZAAT);',
                      ';M(GhVL0A]T41Z8GX(u)T;AX;9XAX;V(DAA828A3GXLCGGGXATXZA;QA;ARXQXTXXCXAIXPXX;',
                      ';GFS9V]I;XVQV5XXAR8AUf)XXXXA6CXCXXRG9XXXDA(j3XXXX[G]TX7AX8A0XXX;',
        
                      #This is a selection of test cases from test_parser.py
                      #Not selected test cases did not raise Exceptions but
                      #would have created units. Not the following scenarios:
                      'ATTGG;', ';ATTGG', ';ATTGGP3', ';ATGGC;3GTCA;' , ';AAA;;;AAA;', 
                      ';;CCAAAA;', ';GTCCCA;;', ';AGCXX;;;3CCCC;', ';AT[[TGG;', ';AT[[TG]G;', 
                      ';N3AX[]5;',';AT//TGG;', ';AT/TGG;', ';N/33[]5;', ';AA.AA;', 
                      ';AATTTC;3.2', ';AAAAA;;1.33TTTCC;', ';AAAAA;1.33;TTTCC;', ';AAP3..1AA;',
                      ';TS0.1AA;', ';(sdf)AAT;', ';F(sdf)1/1[]A!AT2;', ';F1/1[]A(sdf)!AT2;',
                      ';M13/1(a16a)[]1;',';M14/1(aAAa)[]1;',';M15/1(aj1_A_TTC_a)[]1;',
                      ';AS!TGG;', ';AXTTCG;!', '!;AXTTCG;', ';CP3!4TC;', ';N!10[]3;', ';N11[!AA]3;', 
                      ';N12[AA!]3;', ';R!12[]GTC!CCT3;', ';3AAATG;3', ';3AAPTG;3', ';ATGA;5', 
                      ';6ATGAC;', ';16ATGACTTTTTTTTTT;', ';AAPTS2G;', ';ATGP1234AC;', ';M3[]-3;',
                      ';F-1[]A!-2;', '']
        
        #Invalid Characters
        chars = '%+&=`*:<>#?${}\\@^BHJKW'
        for char in chars:
            test_cases.append(';XAAA' + char + 'GGG;')
        
        for test in test_cases:
            c = cell.Cell(test)
            c.read_seq(0)
            self.assertEqual(c.units, {}, '{0}, {1}'.format(test, c.units))
            
class TestCellSplittingUnit(unittest.TestCase):
    
    def test_split_no_units(self):
        '''Check if Splitting unit creates the right number of cells with an equal sequence'''
        tests = [';Y1[]1;', ';Y1[]1;;AAAAAA;', ';Y3[]1;', ';Y50[]1;']
        
        for seq in tests:
            c = cell.Cell()
            ind = c.new_sub_cell(seq)
            c.cells[ind].read_seq()
            c.run_all()
            self.assertEqual(len(c.cells) > 1, True)
            for ind in c.cells.keys():
                self.assertEqual(c.cells[ind].sequence, seq)
                
    def test_split_units_subs(self):
        '''also testing: calling cell by adress'''
        tests = [(';Z1[]1R10[]A!A10L10[]A!A10Q10[]AA10P50SF-100[]20;', ('a','b','c','d','d','e','f','g'))]
        
        for seq, subs in tests:
            c = cell.Cell()
            i = c.new_sub_cell(seq)
            c.cells[i].read_seq()
            c.cells[i].add_subs(subs)
            b_units = c.cells[i].units.values()
            b_subs = c.cells[i].substances.values()
            j = c.split_sub_cell((0,0))
            
            #Check if units after and before splitting are the same
            #They did not age as we performed this via direct call and not via
            #perform action...
            a_units = c.cells[j].units.values() + c.cells[i].units.values()
            a_subs = c.cells[j].substances.values() + c.cells[i].substances.values()
            
            a_units = [tuple([u['type'], u['life'], u['target'], u['tags'], u['action'], u['energy']]) for u in a_units]
            b_units = [tuple([u['type'], u['life'], u['target'], u['tags'], u['action'], u['energy']]) for u in b_units]
            
            self.assertEqual(Counter(a_units), Counter(b_units))
            self.assertEqual(Counter(a_subs), Counter(b_subs))

class TestMutationUnit(unittest.TestCase):
    def test_sequence_change(self):
        test_cases = [';M1[]1XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX;', ';M10[]1AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA;',
                      ';M10[]1XXXXXXXXXXXXXPSSQQQQQQ;',
                      ';M100[]1XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX;',]
        
        for test in test_cases:
            c = cell.Cell(test)
            c.read_seq(0)
            c.run_all()
            self.assertNotEqual(c.sequence, test, 'Don\'t worry if this test fails from time to time...\n\tIt occurs when the mutation unit replaces a letter with that same letter...')
            
class TestSubstances(unittest.TestCase):
    
    def test_substance_addition(self):
        
        test_cases = [(('a', 'b', 'a', 'c'), {0: 'a', 1: 'b', 2: 'a', 3: 'c'}),
                      (('a', 'b', 'heres_another_one', 'and_more'), {0: 'a', 1: 'b', 2: 'heres_another_one', 3: 'and_more'}),
                      ((), {}),
                      (('glucos*', '**', '*_*'), {0: 'glucos*', 1: '**', 2: '*_*'}),
                      (('*',), {0: '*'}), (('**',), {0: '**'}), (('***',), {0: '***'}),
                      (('*****',), {0: '*****'})]
    
        for test in test_cases:
            c = cell.Cell()
            c.add_subs(test[0])
            self.assertEqual(c.substances, test[1])
            
    def test_substance_find_remove(self):
        test_cases = [(('a', 'b', 'a', 'c'), ('a', 'a', 'c'), (('substance', 'a'), ('substance', 'a'), ('substance', 'c')), {1: 'b'}),
                      (('a', 'b', 'a', 'c'), ('a', 'a', 'a'), False, {0: 'a', 1: 'b', 2: 'a', 3:'c'}),
                      (('aa', 'bb', 'a*', 'cc'), ('a*', 'b*', 'a*'), (('substance', 'aa'), ('substance', 'bb'), ('substance', 'a*')), {3:'cc'}),
                      ]
    
        for to_add, to_remove, res, exp in test_cases:
            c = cell.Cell()
            c.add_subs(to_add)
            success = c.find_tags(to_remove, substances = True, remove = True)
            #self.assertEqual(success, res)
            self.assertEqual(c.substances, exp)
            
    def test_only_find(self):
        test_cases = [(('a', 'b', 'a', 'c'), ('a', 'a', 'c'), ('a', 'a', 'c')),
                      (('a', 'b', 'a', 'c'), ('a', 'a', 'a'), False),
                      (('aa', 'bb', 'a*', 'cc'), ('a*', 'b*', 'a*'), ('aa', 'bb', 'a*')),
                      (('*', '*_*', '*__*', '*___*'), ('*','***','****','*****'), ('*', '*_*', '*__*', '*___*'))]
    
        for to_add, to_remove, res in test_cases:
            c = cell.Cell()
            c.add_subs(to_add)
            before = c.substances.items()
            success = c.find_tags(to_remove, substances = True, remove = False)
            #self.assertEqual(success, res)
            self.assertEqual(before, c.substances.items())
            
            
    def test_substance_conversion(self):
        test_cases = [(';V1[](a,b)(c,d)1;', ('a', 'b'), ['c', 'd']), #substance conversion
                      (';V2[](a)(c,d)1;', ('a', 'a'), ['c', 'd', 'c', 'd']),
                      (';V2/1[](a)(c,d)1;', ('a', 'a'), ['a', 'a']),
                      (';V1[]()(c,d)1;', (), ['c', 'd']), #Substance creation
                      (';V1[]()(*,***)1;', (), ['*', '***']), 
                      (';V2[](a,a)()1;', ('a', 'a', 'a', 'a'), []), #Substance destruction
                      (';V3[](a,a)()1;', ('a', 'a', 'a'), ['a']),
                      (';V2[](a,a)1;', ('a', 'a', 'a', 'a'), []),
                      (';V1[](a***,***d)1;', ('abcd', 'abcd'), []),
                      (';V3[](*,*,*)1;', ('*', 'a', 'd'), []), 
                      (';V1[](*,**,*_*)1;', ('a', 'aa', 'c_c'), []),
                      ]
    
        for seq, subs, exp in test_cases:
            c = cell.Cell(seq)
            c.add_subs(subs)
            c.read_seq()
            c.run_all()
            sub = c.substances.values()
            self.assertEqual(sub, exp)
                

class TestChangeSequenceLength(unittest.TestCase):
    def test_protection_change(self):
        test_cases = [(';Q1[]AA1TTAACC;', (8, 3), [(13,14)]),
                      (';Q1[]AA1TTAACC;', (8, 10), [(20,21)]),
                      (';Q1[]AA1TTAACC;', (9, -3), [(7,8)]),
                      (';Q1[]AA1TTAACC;', (12, 10), [(10,11)]),
                      (';Q1[]AA1TTAACC;', (13, -50), [(10,11)])]
        
        for test in test_cases:
            c = cell.Cell(test[0])
            c.read_seq(0)
            c.run_all()
            c.change_sequence_length(test[1][0], test[1][1])
            self.assertEqual(c.pr_pos, test[2])
            
    def test_protection_change_real(self):
        '''!!!STOP THIS TEST FROM FAILING!'''
        print '\n!!!FIX: test_protection_change_real'
        test_cases = [';Q1[]AA2R1[]T!T1TTAA;',
                      ';Q2[]AA2R1[]T!T1AATTAA;',
                      ';Q2[]AA2R1[]T!TT!T1AATTTAA;',
                      ';CAACQ3[]AA2L1[]T!T1T;;TAA;;CAAC;',
                      ';CAACQ3[]AA2L1[]T!T1T;;CAAC;;TAA;',
                      ';CAAC;;TAACCCC;;Q3[]AA2L1[]T!T1T;;CAAC;',
                      ]
        for test in test_cases:
            c = cell.Cell(test)
            c.read_seq(0)
            c.read_seq(13)
            c.run_all()
            pr_pos_before = set(c.pr_pos)
            c.reset_protected()
            pr_pos_next = set(c.pr_pos)
            #self.assertEqual(pr_pos_before, pr_pos_next)
            
class TestCellDisintegration(unittest.TestCase):
    def test_Z_unit(self):
        tests = [';Z1[]1;', ';Z150[]1;']
        
        for seq in tests:
            c = cell.Cell()
            ind = c.new_sub_cell(seq)
            c.cells[ind].read_seq()
            before = [unit for key, unit in c.cells[ind].units.items() if unit['life'] > 1]
            c.run_all()
            
            #only checking sequence for the moment:-(
            self.assertEqual(c.sequence, seq)
    
class TestEnergy(unittest.TestCase):
    def test_energy_consumption(self):
        test_cases = [
                      (';M3/1[]1AAA;', 97), #Check ordinary unit energy consumption
                      (';M1/99[]1AAA;', 1),
                      (';R3/45.001[]A!AT1AATAAT;', 10),
                      (';L3/45.001[]A!AT1AATA;;AT;', 55),
                      
                      (';N10/50[]1F1[]A!AT1SPAATAR1[]C!C1TT;', 50), #Check translation process energy consumption
                      (';N10/10[]1SP5R1[]A!A1;', 50),
                      (';N10/10[]1SP5R1[]A!A1L1[]A!A1;', 0),
                      (';N50/1[]1SP5R1[]A!A1R1[]A!A1R1[]A!A1R1[]A!A1R1[]A!A1R1[]A!A1R1[]A!A1R1[]A!A1R1[]A!A1R1[]A!A1SR1[]A!A1;', 50),#Ten units expressed
                      (';N50/1[]1SP5N2[]1;', 95), #Check that only previously created translation unit executes jobs
                      (';N10/99[]10F3[]A!A1SPAAQ1[]AA1;', 1),
                      (';N10/99[]10F3[]A!A1PAAQ1[]AA1;', 100), #Protected promoter does not get expressed
                      
                      ]
        
        for test_case in test_cases:
            c = cell.Cell(test_case[0])
            c.read_seq(0)
            c.change_energy(100)
            c.run_all()
            self.assertEqual(c.has_energy(total=True), test_case[1], 'Energy consumption failed {0} != {1}, Sequence at end: {2}'. format(c.has_energy(total=True), test_case[1], c.sequence))
    
    def test_energy_change(self):
        
        tests = [(100, -50, 50),
                 (0, 50, 50),
                 (10, -11, 10),
                 (0, 0, 0),
                 (10, -1, 9),
                 (10, -0, 10)]
        
        for start, diff, end in tests:
            c = cell.Cell()
            c.change_energy(start)
            en = c.has_energy(total = True)
            self.assertEqual(en, start)
            c.change_energy(diff)
            en = c.has_energy(total = True)
            self.assertEqual(en, end)

class TestActionList(unittest.TestCase):
    
    def test_get_action_list(self):
        test_cases = [(';R1[]A!A1AAA;', [[[0],'R', 0]]),
                      (';L2[]A!A1AAA;', [[[0],'L', 0], [[0],'L', 0]])
                      ]
        
        for sequence, action_list in test_cases:
            c = cell.Cell(sequence)
            c.read_seq(0)
            res = c.start_round()
            self.assertEqual(res, action_list)
            
    def test_requested_reads_negative(self):
        test_cases = [
                      (';N10[]3F1[]A!A2F-3[]A!A2SPAACCCCM0[]2;'),
                      (';N10[]3F1[]A!A2F-1[]A!A2SPAACCCCM0[]2;'),
                      (';N10[]3F1[]AC!2F-4.0[]A!C2SP3ACCCCM0[]2;'),
                      (';N10[]4F-5.0[]!AAA4SP5AAAN1.0[]2S;')
                      ]
        
        for seq in test_cases:
            c= cell.Cell(seq)
            c.read_seq(0)
            res = c.get_requested_reads()
            self.assertEqual(res, [])
            
    def test_requested_reads(self):
        test_cases = [
                      (';N10[]1F1[]GX!AT1SXXP4GCATXXXSR2[]C!C3;', [22, 22, 22, 22, 22])
                      ]
        
        for seq, exp_reads in test_cases:
            c= cell.Cell(seq)
            c.read_seq(0)
            res = c.get_requested_reads()
            self.assertEqual(res, exp_reads)

class TestProbabilistic(unittest.TestCase):
            
    def test_sample_probability(self):
        test_cases = range(1, 101, 13)
        test_cases += [0.025, 0.201, 4.75, 2.987, 15, 456.123, 0.909]
        test_cases += [-0.8, -4, -38.881, -10.01]
        
        rounds = 10
        
        c = cell.Cell()
        for test in test_cases:
            for i in range(rounds):
                ret = float(c.sample_probability(test))
                if test < 0:
                    self.assertEqual(ret <= 0, test <= 0)
                if test > 0:
                    self.assertEqual(ret >= 0, test >= 0)
            
    def test_sample_normal(self):
        factors = [0, 0.1, 0.25, 0.5, 1, 1.24, 2, 50, 100, 999999.999]
        test_cases = range(1, 101, 21)
        test_cases += [0.025, 0.201, 4.75, 2.987, 999]
        test_cases += [-0.8, -4, -38.881, -10.01]
        
        rounds = 10
        
        for factor in factors:
            c = cell.Cell()
            cell.sample_mode = 'normal'
            cell.normal_factor = factor
            for test in test_cases:
                for i in range(rounds):
                    ret = c.sample_probability(test)
                    if test < 0:
                        self.assertEqual(ret <= 0, test <= 0)
                    elif test > 0:
                        self.assertEqual(ret >= 0, test >= 0)
                    elif test == 0:
                        self.assertEqual(ret, 0)

class TestTranscriptionFactor(unittest.TestCase):
        
    def test_tf_bad(self):
        tests = [(';N10[]1F1[]G!AT1SXXXPGCTXXXR2[]CCC!TT!CCC3;', 0),
                 (';N10[]1F1[]G!AT1SXXPXGCTXXXR2[]C!C3;', 0),
                 (';N10[]1F1[]GX!AT1SXXP4GCATXXXSR2[]C!C3;', 0)]
        
        for seq, num_units in tests:
            c = cell.Cell(seq)
            c.read_seq()
            c.run_all()
            self.assertEqual(len(c.units), num_units)
 
class TestExportUnit(unittest.TestCase):
    '''Implements test methods for the sequence export unit'''
    
    def test_export(self):
        test_cases = [(';E1[]G!AT1SXXXPGATXXXM1[]1;XXTC;',';ATXXXM1[]1;'),        
                      (';E1[]CC!1GGCCTTAATT;',';TTAATT;'),
                      (';E5[]GX!AT1SXXXPGGATPGATN1[]1;',
                       ';ATPGATN1[]1;;ATPGATN1[]1;;ATPGATN1[]1;;ATPGATN1[]1;;ATPGATN1[]1;'),
                      (';E2[]GX!AT1SXXXPGGATPGA;3',';ATPGA;;ATPGA;')]
        
        for seq, exp in test_cases:
            c = cell.Cell()
            ind = c.new_sub_cell(seq)
            c.cells[ind].read_seq()
            c.run_all()
            self.assertEqual(c.sequence, exp, 'Result, Expected \n{0}\n{1}'.format(c.sequence, exp))

class TestSquenceExportImport(unittest.TestCase):
    
    def test_roundtrips(self):
        test_cases = [
                      (';E1[]G!A1I1[]AA2SGGGAAAA;', ';E1[]G!A1I1[]AA2SGGGAAAA;;AAAA;'),
                      (';E1[]G!A1I1[]TT2SGGGACCCCCCTT;', ';E1[]G!A1I1[]TT2SGGGACCCCCCTT;;ACCCCCCTT;'),
                      (';E1[]CC!1I1[]AA3GGCCTTAATT;', ';E1[]CC!1I1[]AA3GGCCTTAATT;;TTAATT;'),
                      (';E1[]G!A1I1[]CC1SGGGAAAA;',';E1[]G!A1I1[]CC1SGGGAAAA;'), #test bad roundtrips (expecting no change in sequence)
                      (';E1[]G!A1I1[]XG1SGGGAAAA;',';E1[]G!A1I1[]XG1SGGGAAAA;'),
                      (';E1[]C!A1I1[]AA1SGGGAAAA;',';E1[]C!A1I1[]AA1SGGGAAAA;')
                      ]
        for seq, exp in test_cases:
            c = cell.Cell()
            ind = c.new_sub_cell(seq)
            c.cells[ind].read_seq()
            c.run_all()
            c.run_all()
            self.assertEqual(c.cells[0].sequence, exp, 'Export import roundtrip failed.\nRes: {0}\nExp: {1}'.format(c.cells[0].sequence, exp))
        
class TestSequenceDegradationUnit(unittest.TestCase):
            
    def test_degradation_unit(self):
        test_cases = [
                      (';D3[]AAAACCCTTT2;;XAAAACCCTTTX;', ';D3[]AAAACCCTTT2;'),
                      (';D1[]XXCCC2;;CCCA;;GTCCCA;', ';D1[]XXCCC2;;CCCA;'),
                      (';D1[]XXCAC2;;CCCA;;GTCCCA;', ';D1[]XXCAC2;;CCCA;;GTCCCA;'),
                      (';D1[]TTTT3Q1[]GGGG3GGGGTTTT;', ''),
                      (';D1[]GT3Q1[]AGGG3AGGGGTTTT;', ''),
                      (';D1[]G3Q1[]GGGG3TTTT;', ';D1[]G3Q1[]GGGG3TTTT;')
                      ]
        
        for seq, exp in test_cases:
            c = cell.Cell(seq)
            c.read_seq()
            c.run_all()
            self.assertEqual(c.sequence, exp)

class TestProtectionUnit(unittest.TestCase):
    '''Testing the sequence protection unit.'''
    def test_pos_protected(self):
        '''All positions should be protected'''
        
        tests = [
                 (';Q1[]AACC3AACC;', range(10,14), True),
                 (';Q1[]AACC3P3AACC;', range(12,16), True),
                 (';Q3[]AGCC3AGCCAGCCAGCC;', range(10,22), True),
                 (';Q2[]AXCC3AACCAGCC;', range(10,18), True),
                 (';Q2[]AXCC3AACCAGCC;', range(10,18), True),
                 (';Q2[]AXCC3AACCAGCC;', range(0,10), False),
                 (';Q2[]AXCC3AACCAGCCXXXXXXX;', range(18,25), False),
                 (';Q2[]AXCC3AACCGGGGTTTC;', range(18,25), False),
                 (';Q1[]GGGGTTTT3GGGGTTTTO1[]G!GTT3;', range(14), False),
                 (';Q1[]GGGGTTTT3GGGGTTTTO1[]G!GTT3;', range(14, 22), True),
                 (';Q1[]GGGG3D1[]GTTT3GGGGTTTT;', range(19, 23), True),
                 (';Q1[]GGGG3D1[]GTTT3GGGGTTTT;', range(23, 26), False),
                 (';Q1[GGGGTTTT]GGGGT3;', range(4, 9), True),
                 (';Q1[GGGGTTTT]GGTTT3;', range(6, 11), True),
                 (';Q1[GGGGTTTT]GGTTT3;', range(6), False),
                 (';Q1[GGGGTTTT]GGTTT3;', range(11, 19), False)
                 ]
        
        for seq, test_pos, exp in tests:
            c = cell.Cell(seq)
            c.read_seq()
            c.run_all()
            for pos in test_pos:
                self.assertEqual(c.pos_protected(pos), exp)
    
    def test_protection_sequence(self):
        '''expecting no change in sequence'''
        
        tests = [';Q1[]AACC3R3[]AA!CC2;;AACC;', # restriction unit
                 ';Q1[]AACC3R3[]A!ACC2;;AACC;',
                 ';Q1[]AACC3R3[]AA!CC2;;AACC;',
                 ';Q1[]AACC3R3[]AA!CC!2;;AACC;',
                 ';Q1[]AACC3O3[]A!ACC2;;AACC;', # copy unit
                ';Q3[]XACC3O3[]X!ACC2;;AACC;;GACC;;TACC;',
                ';Q1[]GGGGTTTT3O1[]GGG!GTT3GGGGTTTT;',
                ';Q1[]GGGGTTTT3D1[]GGTT3GGGGTTTT;', # degradation unit
                ';Q1[]TTTT3D1[]GGTT3GGGGTTTT;',
                ';Q1[]TTTT3D1[]GGT3GGGGTTTT;',
                ';Q1[]GGGG3D1[]GTTTT3GGGGTTTT;'
                 ]
        
        for seq in tests:
            c = cell.Cell(seq)
            c.read_seq()
            c.run_all()
            self.assertEqual(c.sequence, seq)
            
    def test_protection_multiple_rounds(self):
        '''expecting no change in sequence'''
        tests = [(';Q1[]AACC5XXN10[]5SP1AACCTTCCR1[]TT!CC1;', 3), #translation stop
                 (';Q1[]AACC6XXN10[]5SP1R1[]AA!CC1AACCTTCCR1[]TT!CC1;', 5)]
        
        for seq, num_rounds in tests:
            c = cell.Cell(seq)
            c.read_seq()
            for i in range(num_rounds):
                c.run_all()
            self.assertEqual(c.sequence, seq)
            
class TestSequenceImportUnit(unittest.TestCase):
    '''Implements test methods for the sequence import unit'''
    
    def test_import(self):
        #Format: initial sequence, initial external sequence, expected sequence, expected external sequence
        test_cases = [(';I1.0[]GAT1;', ';GATTACA;', ';I1.0[]GAT1;;GATTACA;', ''), #test import good
                      (';I1[]GAT1XXAA;2', ';GATTACA;', ';I1[]GAT1XXAA;2;GATTACA;' , ''),
                      (';I1[]CC1XXAA;', ';TTTTCC;', ';I1[]CC1XXAA;;TTTTCC;' , ''),
                      (';I1[]CC1XXAA;', ';GGGTAC;;TTTTCC;', ';I1[]CC1XXAA;;TTTTCC;' , ';GGGTAC;'),
                      (';I1[]CC1XXAA;', ';GGGTAC;;TTTTCC;', ';I1[]CC1XXAA;;TTTTCC;', ';GGGTAC;'),
                      (';I1[]CC1XXAA;', ';GGGTAC;;TTTTCC;;GGGTAC;', ';I1[]CC1XXAA;;TTTTCC;', ';GGGTAC;;GGGTAC;'),
                      (';I1[]CC1XXAA;', ';GGGTAC;;TTTTCC;;TTXTAC;', ';I1[]CC1XXAA;;TTTTCC;', ';GGGTAC;;TTXTAC;'),
                      (';I1.0[]C!C1XXAA;', ';TTTTCC;', ';I1.0[]C!C1XXAA;;TTTTCC;', ''),         # Ignore the ! 
                      (';I1[]!TT1XXAA;', ';TTTTCC;', ';I1[]!TT1XXAA;;TTTTCC;', ''),
                      (';I1[]TC!1XXAA;', ';TTTTCC;', ';I1[]TC!1XXAA;;TTTTCC;', ''),
                      (';I1[]XCC1;', ';GATTACA;', ';I1[]XCC1;', ';GATTACA;'),                   #test bad imports
                      (';I1[]XG1XXAA;2', ';GATTACA;', ';I1[]XG1XXAA;2', ';GATTACA;'),
                      (';I1[]G1XXAA;', ';TTTTCC;', ';I1[]G1XXAA;', ';TTTTCC;'),
                      ]
        for start1, start2, end1, end2  in test_cases:
            c = cell.Cell(start2)
            ind = c.new_sub_cell(start1)
            c.cells[ind].read_seq()
            c.run_all()
            self.assertEqual(c.cells[ind].sequence, end1)
            self.assertEqual(c.sequence, end2)
    
class TestFindOccurence(unittest.TestCase):
    def test_occurence(self):
        c = cell.Cell(';CAAAAGTCTGACAAAAGTATGA;;GGTCA;3')
        
        tests = [('AAAA', 2, [2, 13]),
                 ('AGTC', 1, [5]),
                 ('A;', 3, [22, 22, 22]),
                 (';3', 1, [30]),
                 (';;', 1, [23])]
        
        for target, max, exp in tests:
            occ = c.find_occurences('own', 0, len(c.sequence), target, max)
            self.assertEqual(set(occ), set(exp))
    
    def test_find_occ_bad(self):
        test_cases = [(';E1[]CC!1I1[]AA3GGCCTTAATT;', 'CC', [18]),
                      (';E1.0[]!GG1I1[]AA3GGCCTTAATT;', 'GG', [18]),
                      (';E1[]!GG1I1[]AA3GGCCTTAATT;', 'AA', [22]),
                      (';E1.0[]!GG1I1[]AA3GGAACCAATTAATT;', 'AA', [20, 24, 28]),
                      (';Q3.0[]AGCC3AGCCAGCCAGCC;', 'AGCC', [12, 16, 20])
                      ]
        
        for seq, target, exp in test_cases:
            c = cell.Cell(seq)
            occ = c.find_occurences('own', 0, len(seq), target, 10)
            occ = set(occ)
            self.assertEqual(occ, set(exp), 'Result, Expected \n{0}\n{1}'.format(occ, set(exp)))
        
    def test_find_occ_prom(self):
        test_cases = [
                      (';P1AAAT;', 'AA', [3]),
                      (';P13.0AAGT;', 'XA', [6]),
                      (';ATGPCP1AAAT;', 'AA', [8]),
                      (';AP0.334XXXT;', 'XXX', [8]),
                      (';AP111.134XXXT;', 'XXX', [10]),
                      (';AAP2.223ACCC;', 'ACCC', [9])
                      ]
        
        for seq, target, exp in test_cases:
            c= cell.Cell(seq)
            occ = c.find_occurences('prom', 0, len(seq), target, 10)
            occ = set(occ)
            self.assertEqual(occ, set(exp), 'Result, Expected \n{0}\n{1}'.format(occ, set(exp)))

class TestRestriction(unittest.TestCase):
    
    def test_restriction(self):
        tests = [
                 #test blunt
                 (';R10[]A!G3AAAGGG;',';R10[]A!G3AAA;;GGG;'),                                   
                 (';R10[]AA!GG3AACGGGAAAGGGAAAGGG;',';R10[]AA!GG3AACGGGAAA;;GGGAAA;;GGG;'),
                 #test blunt multiple
                 (';R10[]A!G3AAAGGGAAAGGG;',';R10[]A!G3AAA;;GGGAAA;;GGG;'),
                 (';R10[]GGA!AGG1GGAAGGGGAAGGGGAAGGGGAAGGGGAAGGGGAAGGGGAAGGGGAAGGGGAAGG;',
                  ';R10[]GGA!AGG1GGA;;AGGGGA;;AGGGGA;;AGGGGA;;AGGGGA;;AGGGGA;;AGGGGA;;AGGGGA;;AGGGGA;;AGG;'),
                 #test overlap
                 (';R10[]A!AG!G3;AAGG;',';R10[]A!AG!G3;AAG;2;2AGG;'),
                 (';R10[]A!GG!3AAAAGGGG;',';R10[]A!GG!3AAAAGG;2;2GGGG;'),
                 (';R1[]A!AGX!G1AAGAGGGG;',';R1[]A!AGX!G1AAGA;3;3AGAGGGG;'),
                 (';R2[]A!A!GG3AAAGGGGAAAGGG;',';R2[]A!A!GG3AAA;1;1AGGGGAAA;1;1AGGG;'),
                 #test overlap bad
                 (';R10[]A!AG!G3;;3AAAGGG;',';R10[]A!AG!G3;;3AAAGGG;'),
                 (';R10[]!XAG!G3;TAAGGCG;4',';R10[]!XAG!G3;TAAGGCG;4'),
                 (';R20[]!XXXX!3;;3TAATGGCG;2',';R20[]!XXXX!3;;3TAATGGCG;2'),
                 ]
        
        for seq, exp in tests:
            c= cell.Cell(seq)
            c.read_seq()
            c.run_all()
            self.assertEqual(c.sequence, exp)
            
class TestLigation(unittest.TestCase):
    def test_ligation_sequence(self):
        tests = [
            # test blunt adjacent location
            (';L10[]AA!XA3AAA;;AAA;',';L10[]AA!XA3AAAAAA;'),
            (';L10[]G!X3AAG;;AAA;',';L10[]G!X3AAGAAA;'),
            (';L10[]G!A3AAG;;AAA;;2AATC;',';L10[]G!A3AAGAAA;;2AATC;'),
        
            #test_blunt_simple move
            (';L10[]G!C3AAG;;AAA;;CCC;',';L10[]G!C3AAGCCC;;AAA;'),
            (';L10[]AG!XG3AAG;;CCG;;GGG;',';L10[]AG!XG3AAGGGG;;CCG;'),
            (';L10[]AG!XG3AAG;;2CCT;;GGG;;ACTAA;3',';L10[]AG!XG3AAGGGG;;ACTAA;3;2CCT;'),
                 
            #Complex Move
            (';L10[]AG!CC3;;CCG;;TTT;;CAG;;TCT;',';L10[]AG!CC3;;TTT;;CAGCCG;;TCT;'),
            (';L10[]AXG!CXC3;;CGC;;XXG;;ACG;;AAC;',';L10[]AXG!CXC3;;XXG;;ACGCGC;;AAC;'),
            
            #Complex move mixed with blunt ends
            (';L10[]AXG!CXC3;;CGC;;2XXG;;2ACG;;AAC;2;2GTAAC;',
             ';L10[]AXG!CXC3;;2XXG;;2ACGCGC;;AAC;2;2GTAAC;'),
            
            #test_overlap_adj
            (';L10[]!XX!3AAA;2;2AAA;',';L10[]!XX!3AAAA;'),
            (';L10[]!GA!3CGA;2;2GAX;',';L10[]!GA!3CGAX;'),
            (';L10[]!GA!G3SGAACTCG;;2GATCA;;CCCGA;2;2GAG;',
             ';L10[]!GA!G3SGAACTCG;;2GATCA;;CCCGAG;'),
            (';L10[]CC!GA!G3SGAACTCG;;2GATCA;;CCCGA;2;2GAG;',
             ';L10[]CC!GA!G3SGAACTCG;;2GATCA;;CCCGAG;'),
            (';L10[]CC!GA!GG3SGAACTCG;;2GATCA;;CCCGA;2;2GAGG;',
             ';L10[]CC!GA!GG3SGAACTCG;;2GATCA;;CCCGAGG;'),
            
            #test_overlap_simple
            (';L10[]!GA!3CGA;2;GAX;;2GAX;',';L10[]!GA!3CGAX;;GAX;'),
            (';L10[]!GA!3CGA;2;GAA;;3GATTACA;;AAAX;;2GAX;',';L10[]!GA!3CGAX;;GAA;;3GATTACA;;AAAX;'),
            
            #test_overlap_complex
            (';L10[]!GA!3CCC;;2GATTAGA;;2CCCCC;2;GAGAGA;2',
             ';L10[]!GA!3CCC;;2CCCCC;2;GAGAGATTAGA;'),
            (';L10[]!CC!A3SAAAG;2;2CCAG;;3TTTTT;;GGCC;2;TTGA;',
             ';L10[]!CC!A3SAAAG;2;3TTTTT;;GGCCAG;;TTGA;'),
            (';L10[]!CC!A3S;;2AAAG;2;2CCAG;;3TT;;TTT;;1GGCC;2;TTGA;',
             ';L10[]!CC!A3S;;2AAAG;2;3TT;;TTT;;1GGCCAG;;TTGA;'),
            
            #test_multiple ligations in one round
            (';L10[]A!A3AAA;;AAA;;AAA;;AAA;;AAA;;AAA;;AAA;;AAA;;AAA;;AAA;;AAA;',
             ';L10[]A!A3AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA;'),
            (';L12[]GGA!XGG3;;AGGGGA;;AGGGGA;;AGGGGA;;AGGGGA;;AGGGGA;;AGGGGA;;AGGGGA;;AGGGGA;;AGGGGA;',
             ';L12[]GGA!XGG3;;AGGGGAAGGGGAAGGGGAAGGGGAAGGGGAAGGGGAAGGGGAAGGGGAAGGGGA;'),
            (';L10[]GG!AX!GG3;;2AAGGGGAA;2;2AAGGGGAA;2;2AAGGGGAA;2;2AAGGGGAA;2',
             ';L10[]GG!AX!GG3;;2AAGGGGAAGGGGAAGGGGAAGGGGAA;2'),
            ]
            
        for seq, exp in tests:
            c = cell.Cell(seq)
            c.read_seq()
            c.run_all()
            self.assertEqual(c.sequence, exp)

class TestConstProm(unittest.TestCase):
    '''Test for the constitutive promoter'''
    
    def test_promoter_repression(self):
        '''Repress promoter with transcription factors'''
        
        tests = [(';N10[]3F-5.0[]!AAA3SP5AAAN3[]7SR15[]GATTA!CA5GGTAC;', 2),
                 (';N10[]3F-10.3[]CCT!T3SP10CCTTN3[]7R15[]GATTA!CA5GGTAC;', 2)]
        for seq, len_units in tests:
            c = cell.Cell(seq)
            c.read_seq()
            c.run_all()
            self.assertEqual(len(c.units), len_units)
        
class TestCopyUnit(unittest.TestCase):
    
    def test_copy_simple(self):
        test_cases = [(';O1[]AAA!CCC2SAAACCCTTTGGG;', ';O1[]AAA!CCC2SAAACCCTTTGGG;;CCCTTTGGG;'),
                      (';O1[]AAA!CCC2SAAACCCTTTGGG;3', ';O1[]AAA!CCC2SAAACCCTTTGGG;3;CCCTTTGGG;'),
                      (';O1[]AX!ACCC2;2AAACCCTTTGGGCC;12', ';O1[]AX!ACCC2;2AAACCCTTTGGGCC;12;ACCCTTTGGGCC;'),
                      (';O3[]AA!CC2;AAAACCCC;', ';O3[]AA!CC2;AAAACCCC;;CCCC;;CCCC;;CCCC;'),
                      (';O5[]A!G1AGAAC;', ';O5[]A!G1AGAAC;;GAAC;;GAAC;;GAAC;;GAAC;;GAAC;')
                      ]
        
        for seq, exp in test_cases:
            c = cell.Cell(seq)
            c.read_seq()
            c.run_all()
            self.assertEqual(c.sequence, exp)
            
class TestUnitNeedsTag(unittest.TestCase):
    def test_needs_substances(self):
        tests = [(';R1[]A!A(a)1AA;', ('a',), ';R1[]A!A(a)1A;;A;'),
                 (';R1[]A!A(a)1AA;', ('b',), ';R1[]A!A(a)1AA;'),
                 (';R1[]A!A(a,b,glucose_)1AA;', ('b', 'a', 'glucose_'), ';R1[]A!A(a,b,glucose_)1A;;A;'),
                 (';R1[]A!A(_glucose_)1AA;', ('_glucose_'), ';R1[]A!A(_glucose_)1A;;A;'),
                 (';R1[]A!A(_glucose_)1AA;', ('glucose_', '***', 'mister_man'), ';R1[]A!A(_glucose_)1AA;')
                 ]
        
        for seq, subs, res in tests:
            c = cell.Cell(seq)
            c.add_subs(subs)
            before = c.substances
            c.read_seq()
            c.run_all()
            self.assertEqual(c.sequence, res)
            self.assertEqual(c.substances, before)
            
    def test_needs_unit(self):
        tests = [(';R1[]A!A(a)1N1(a)[]1AA;', ';R1[]A!A(a)1N1(a)[]1A;;A;'),
                 (';R1[]A!A(a)1N1(b)[]1AA;', ';R1[]A!A(a)1N1(b)[]1AA;'),
                 (';R1[]A!A(a,b)1N1(a,b)[]1AA;', ';R1[]A!A(a,b)1N1(a,b)[]1AA;')
                 ]
        
        for seq, res in tests:
            c = cell.Cell(seq)
            c.read_seq()
            c.run_all()
            self.assertEqual(c.sequence, res)
            
class TestTranslationUnit(unittest.TestCase):
    
    def test_translation_list(self):
        test_cases = [ (';N10[]1N5[]1;', {0: 10, 1: 5}),
                       (';N10[]1SP5M1[]1;', {0: 5}),
                       (';N10[]1SP10N1[]1;', {0: 0}),
                       (';N10[]1N10[]1N10[]1N10[]1N10[]1;', {0: 10, 1: 10, 2: 10, 3: 10, 4: 10})]
        
        for seq, trans_list in test_cases:
            c = cell.Cell(seq)
            c.read_seq()
            c.run_all()
            self.assertEqual(c.tr_count, trans_list)

class TestTransportUnit(unittest.TestCase):
    
    def test_unit_export(self):
        test_cases = [(';R1(res)[]A!A10U10[](***)2;', {0: {'type': 'R','action': 1,'energy': 0,'target': 'A!A','life': 9,'tags': ('res',)}}),
                      (';U1(self)[](se**)2;', {0: {'type': 'U','action': 1,'energy': 0,'target': ('se**',),'life': 1,'tags': ('self',)}}),
                      (';U1(other_tag,yet_another,self)[](self)2;', {0: {'type': 'U','action': 1,'energy': 0,'target': ('self',),'life': 1,'tags': ('other_tag', 'yet_another','self')}}),
                      (';L10(ligation_unit_one)[]A!A10L10(ligation_unit_one)[]A!A10U10[](ligation_unit_one)2;', 
                       {0: {'type': 'L','action': 10,'energy': 0,'target': 'A!A','life': 9,'tags': ('ligation_unit_one',)},
                        1: {'type': 'L','action': 10,'energy': 0,'target': 'A!A','life': 9,'tags': ('ligation_unit_one',)}})] 
        
        for seq, units in test_cases:
            c = cell.Cell()
            ind = c.new_sub_cell(seq)
            c.cells[ind].read_seq()
            c.run_all()
            self.assertEqual(c.units, units)
            
    def test_unit_import(self):
        test_cases = [(';R1(res)[]A!A10;', ';U-10[](res)1;', {1: {'type': 'R','action': 1,'energy': 0,'target': 'A!A','life': 9,'tags': ('res',)}}),
                      (';L10(ligation_unit_one)[]A!A10L10(ligation_unit_one)[]A!A10;', ';U-10[](ligation_unit****)1;', 
                       {1: {'type': 'L','action': 10,'energy': 0,'target': 'A!A','life': 9,'tags': ('ligation_unit_one',)},
                        2: {'type': 'L','action': 10,'energy': 0,'target': 'A!A','life': 9,'tags': ('ligation_unit_one',)}})] 
        
        for seq_ext, seq, units in test_cases:
            c = cell.Cell(seq_ext)
            ind = c.new_sub_cell(seq)
            c.read_seq()
            c.cells[ind].read_seq()
            c.run_all()
            self.assertEqual(c.cells[ind].units, units)
            
    def test_substance_export(self):
        test_cases = [(';U10[](a)2;', ('a', 'a', 'a', 'a', 'a'), ['a', 'a', 'a', 'a', 'a'], []),
                      (';U10[](b**,a)2;', ('a', 'a', 'a', 'a', 'a'), ['a', 'a', 'a', 'a', 'a'], []),
                      (';U4[](*)2;', ('*', '*', '*', '*', '*'), ['*', '*', '*', '*'], ['*']),
                      (';U4[](*,**)2;', ('*', '**', '*', '**','***','abcd'), ['**', '**', '*', '*'], ['***', 'abcd']),
                      ] 
        
        for seq, subs, res, res1 in test_cases:
            c = cell.Cell()
            ind = c.new_sub_cell(seq)
            c.cells[ind].read_seq()
            c.cells[ind].add_subs(subs)
            c.run_all()
            
            self.assertEqual(c.cells[ind].substances.values(), res1)
            self.assertEqual(Counter(c.substances.values()), Counter(res))
            
    def test_substance_import(self):
        test_cases = [(';U-6[](a)2;', ('a', 'a', 'a', 'a', 'a'), {0: 'a', 1: 'a', 2: 'a', 3: 'a', 4: 'a'}),
                      (';U-10.889[](b,a)2;', ('a', 'a', 'a', 'a', 'a'), {0: 'a', 1: 'a', 2: 'a', 3: 'a', 4: 'a'})] 
        
        for seq, subs, res in test_cases:
            c = cell.Cell()
            ind = c.new_sub_cell(seq)
            c.add_subs(subs)
            c.cells[ind].read_seq()
            c.run_all()
            self.assertEqual(c.substances, {})
            self.assertEqual(c.cells[ind].substances, res)
            
if __name__ == '__main__':
    unittest.main()