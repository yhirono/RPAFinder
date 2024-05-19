import numpy as np
import unittest
from rpa_finder.reaction_system import ReactionSystem, dimcoker


class TestReactionSystem(unittest.TestCase):
    def setUp(self):
        pass

    def test_dimcoker(self):
        mat = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0]])
        assert dimcoker(mat) == 1

    def test_find_reactions_to_add(self):
        system = ReactionSystem(
            [
                ['','"x"'],
                ['"x"','"y"'],
                ['"y"',''],
            ]
        )
        assert system.find_reactions_to_add([0],[]) == [1]
        assert len(system.find_reactions_to_add([0],[1]))==0
        assert system.find_reactions_to_add([0,1],[]) == [1,2]

    def test_compute_influence_index(self):        
        system1 = ReactionSystem(
            [
                ['0','"x"'],
                ['"x"','"y"'],
                ['"y"','"z"'],
                ['"z"', '0']
            ]
        )

        assert system1.compute_influence_index([[0],[1]]) == 0  
        assert system1.compute_influence_index([[],[]]) == 0
        assert system1.compute_influence_index([[0,1,2],[0,1,2,3]]) == 0

    def test_is_output_complete(self):
        system = ReactionSystem(
            [
                ['"A"+"B"','"B"'],
                ['"C"','"D"'],
                ['"A"','"D"'],
                ['2"A"+1"B"', '"D"']
            ]
        )
        assert system.is_output_complete([[0,1],[0,2,3]]) == True
        assert system.is_output_complete([[0,1],[0]])== False


    def test_enumerate_labeled_buffering_structures_1(self):
        network1 = [
                ['','"x"'],
                ['"x"','"y"'],
                ['"y"',''],
            ]
        assert ReactionSystem(network1).enumerate_labeled_buffering_structures() \
            == [[[0], [0, 1], [0, 1, 2], []], [[1], [0], [], [1]], [[2], [1], [], [2]]]


    def test_enumerate_labeled_buffering_structures_2(self):
        network2 = [
            ['', '"v1"'],
            ['','"v2"'],
            ['"v1"','"v3"'],
            ['"v2"','"v4"'],
            ['"v3"','"v4"'],
            ['"v3"','"v5"'],
            ['"v4"','"v6"'],
            ['"v6"','"v5"'],
            ['"v5"','"v7"'],
            ['"v6"','"v8"'],
            ['"v7"','"v3"'],
            ['"v7"+"v8"','"v9"'],
            ['"v6"',''],
            ['"v9"','']
        ]
        assert ReactionSystem(network2).enumerate_labeled_buffering_structures() \
            == [[[0], [0, 2, 3, 4, 5, 6, 7, 8], [0, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], []],
                [[1], [1, 2, 3, 4, 5, 6, 7, 8], [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], []],
                [[2], [0], [], [2]],
                [[3], [1], [], [3]],
                [[4], [2, 4, 6, 7], [5, 8, 10], [4, 11]],
                [[5], [4, 6, 7], [5, 8, 10], [11]],
                [[6], [3], [], [6]],
                [[7], [2, 3, 4, 6, 7], [4, 5, 6, 7, 8, 10], [11]],
                [[8], [4], [], [8]],
                [[9, 12], [2, 3, 4, 5, 6, 7, 8], [4, 5, 6, 7, 8, 9, 10, 11, 12, 13], []],
                [[10], [6, 7], [], [10, 11]],
                [[11], [7], [], [11]],
                [[13], [8], [], [13]]
                ]


    def test_enumerate_labeled_buffering_structures_3(self):
        network3 = [
                ['', '"v1"'],
                ['"v1"','"v2"'],
                ['"v2"','"v3"'],
                ['"v3"','"v1"'],
                ['"v2"','']
            ]        
        assert ReactionSystem(network3).enumerate_labeled_buffering_structures() \
            == [[[0], [0, 1, 2], [0, 1, 2, 3, 4], []],
                [[1], [0], [], [1]],
                [[2], [0, 2], [1, 2, 3], []],
                [[3], [2], [], [3]],
                [[4], [0, 1, 2], [1, 2, 3], [4]]]

    def test_enumerate_labeled_buffering_structures_4(self):
        network4 = [
            ['', '"v1"'],
            ['2"v1"','"v2"'],
            ['"v2"','"v3"'],
            ['"v3"','"v1"'],
            ['"v2"','']
        ]
        assert ReactionSystem(network4).enumerate_labeled_buffering_structures() \
            == [[[0], [0, 1, 2], [0, 1, 2, 3, 4], []],
                [[1], [0], [], [1]],
                [[2, 4], [0, 1, 2], [1, 2, 3, 4], []],
                [[3], [2], [], [3]]
            ]

    def test_enumerate_labeled_buffering_structures_5(self):
        network5 = [
            ['"z1"', '"z1"+"x"'],
            ['"z1"+"z2"',''],
            ['','"z1"'],
            ['2"x"','2"x"+"z2"'],
            ['"x"',''],
            ['"y"','"y"+"z1"'],
            ['','"y"'],
            ['"y"','']
        ]
        assert ReactionSystem(network5).enumerate_labeled_buffering_structures() \
            == [[[0], [2, 3], [], [0, 1]],
                [[1], [3], [], [1]],
                [[2], [0, 2, 3], [0, 1, 2, 3, 4], []],
                [[3], [0, 2, 3], [0, 4], [1, 3]],
                [[4], [2, 3], [0, 4], [1]],
                [[5], [0, 2, 3], [0, 1, 3, 4, 5], []],
                [[6], [0, 1, 2, 3], [0, 1, 3, 4, 5, 6, 7], []],
                [[7], [0, 1, 2, 3], [0, 1, 3, 4, 5], [7]]
                ]

    def test_enumerate_labeled_buffering_structures_6(self):
        network6 = [
            ['', '"v1"'],
            ['"v1"','"v2"'],
            ['"v2"',''],
            ['"v1"+"v2"','"v3"+"v4"'],
            ['"v3"+"v4"','"v1"+"v2"']
        ]
        assert ReactionSystem(network6).enumerate_labeled_buffering_structures() \
            == [[[0], [0, 1, 2, 3], [0, 1, 2, 3, 4], []],
                [[1], [0, 2, 3], [3, 4], [1]],
                [[2], [1, 2, 3], [3, 4], [2]],
                [[3], [2, 3], [3, 4], []],
                [[4, 5], [2, 3], [], [4]]
            ]

    def test_enumerate_labeled_buffering_structures_7(self):
        network7 = [
            ['"xms"', '"xm"'],
            ['"xm"','"xms"'],
            ['','"xm"'],
            ['','"xms"'],
            ['"xms"','']
        ]
        assert ReactionSystem(network7).enumerate_labeled_buffering_structures() \
            == [[[0], [0], [0, 1], []],
                [[1], [0], [], [1]],
                [[2], [0, 1], [0, 1, 2, 4], []],
                [[3], [0, 1], [0, 1, 3, 4], []],
                [[4], [0, 1], [0, 1], [4]]]

    def test_enumerate_labeled_buffering_structures_8(self):
        network8 = [
            ['"Glucose"', '"G6P"'],
            ['"G6P"', '"F6P"'],['"F6P"', '"G6P"'],
            ['"F6P"', '"F16P"'],
            ['"F16P"', '"G3P" + "DHAP"'],
            ['"DHAP"', '"G3P"'],
            ['"G3P"', '"PGP"'],
            ['"PGP"', '"3PG"'], ['"3PG"', '"PGP"'],
            ['"3PG"', '"2PG"'], ['"2PG"', '"3PG"'],
            ['"2PG"', '"PEP"'], ['"PEP"', '"2PG"'],
            ['"PEP"', '"PYR"'],
            ['"G6P"', '"PG6"'],
            ['"PG6"', '"Ru5P" + "CO2"'],
            ['"Ru5P"', '"X5P"'],
            ['"Ru5P"', '"R5P"'],
            ['"X5P" + "R5P"', '"G3P" + "S7P"'],
            ['"G3P" + "S7P"', '"X5P" + "R5P"'],
            ['"G3P" + "S7P"', '"F6P" + "E4P"'],
            ['"F6P" + "E4P"', '"G3P" + "S7P"'],
            ['"X5P" + "E4P"', '"F6P" + "G3P"'],
            ['"F6P" + "G3P"', '"X5P" + "E4P"'],
            ['"PG6"', '"G3P" + "PYR"'],
            ['"PYR"', '"Acetal" + "CO2"'],
            ['"Acetal"', '"Ethanol"'],
            ['"Ethanol"', '"Acetal"'],
            ['"R5P"', ''],
            ['"CO2"', ''],
            ['', '"Glucose"'],
            ['"Ethanol"', ''],
            ['"Acetal"', ''],
            ['"PYR"', '"Ala"'],
            ['"Ala"', '"PYR"'],
            ['"Ala"', ''],
        ]
        assert ReactionSystem(network8).enumerate_labeled_buffering_structures() \
            == [[[0], [12], [], [0]],
                [[1, 2, 3, 14, 15, 16, 17, 24],
                [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 20],
                [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,31,32,33,34,35],
                []],
                [[4], [8], [], [4]],
                [[5], [5], [], [5]],
                [[6], [6, 10, 19, 20], [18, 19, 20, 21, 22, 23], [6]],
                [[7], [15], [], [7]],
                [[8], [15], [7, 8], []],
                [[9], [1, 15], [7, 8], [9]],
                [[10], [1, 15], [7, 8, 9, 10], []],
                [[11], [0, 1, 15], [7, 8, 9, 10], [11]],
                [[12], [0, 1, 15], [7, 8, 9, 10, 11, 12], []],
                [[13], [0, 1, 13, 15], [7, 8, 9, 10, 11, 12], [13]],
                [[18, 19, 20, 21, 22], [6, 19, 20], [18, 19, 20, 21], [22]],
                [[23], [6, 19, 20], [18, 19, 20, 21, 22, 23], []],
                [[25, 33, 34, 35],
                [2, 3, 4, 7, 16],
                [25, 26, 27, 29, 31, 32, 33, 34, 35],
                []],
                [[26, 27, 31, 32], [2, 7], [26, 27, 31, 32], []],
                [[28], [6, 17, 19, 20], [18, 19, 20, 21], [22, 28]],
                [[29], [4], [], [29]],
                [[30],
                [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
                [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35],
                []]]



    def test_get_emergent_conserved_quantities(self):
        network1 = [
                ['','"x"'],
                ['"x"','"y"'],
                ['"y"',''],
            ]
        assert len(ReactionSystem(network1).get_emergent_conserved_quantities([[],[]])) == 0
        assert len(ReactionSystem(network1).get_emergent_conserved_quantities([[0],[]])) == 1
        assert len(ReactionSystem(network1).get_emergent_conserved_quantities([[0,1],[1]])) == 1
        assert len(ReactionSystem(network1).get_emergent_conserved_quantities([[0,1],[0,1,2]])) == 0

    def test_get_emergent_conserved_quantities_symb(self):
        network1 = [
                ['','"x"'],
                ['"x"','"y"'],
                ['"y"',''],
            ]
        assert ReactionSystem(network1).get_emergent_conserved_quantities_symb([[],[]]).shape[0] == 0
        assert ReactionSystem(network1).get_emergent_conserved_quantities_symb([[0],[]]).shape[0] == 1
        assert ReactionSystem(network1).get_emergent_conserved_quantities_symb([[0,1],[1]]).shape[0] == 1
        assert ReactionSystem(network1).get_emergent_conserved_quantities_symb([[0,1],[0,1,2]]).shape[0] == 0


    def test_get_lost_conserved_quantities(self):
        network1 = [
                ['"x"','"y"'],
                ['"y"','"z"'],
                ['"z"','"x"'],
            ]
    
        assert len(ReactionSystem(network1).get_lost_conserved_quantities([[],[]])) == 0
        assert len(ReactionSystem(network1).get_lost_conserved_quantities([[0],[]])) == 0
        assert len(ReactionSystem(network1).get_lost_conserved_quantities([[0],[0]])) == 1

    def test_get_lost_conserved_quantities_symb(self):
        network1 = [
                ['"x"','"y"'],
                ['"y"','"z"'],
                ['"z"','"x"'],
            ]
    
        sys1 = ReactionSystem(network1)

        assert sys1.get_lost_conserved_quantities_symb([[],[]]).shape[0] == 0
        assert sys1.get_lost_conserved_quantities_symb([[0],[]]).shape[0] == 0
        assert sys1.get_lost_conserved_quantities_symb([[0],[0]]).shape[0] == 1

    def test_get_emergent_cycles(self):
        network1 = [
                ['"x"','"y"'],
                ['"y"','"z"'],
                ['"z"','"w"']
            ]
        assert len(ReactionSystem(network1).get_emergent_cycles([[],[]])) == 0
        assert len(ReactionSystem(network1).get_emergent_cycles([[],[0]])) == 1
        assert len(ReactionSystem(network1).get_emergent_cycles([[1],[0,1]])) == 1
        assert len(ReactionSystem(network1).get_emergent_cycles([[0,1,2],[0,1,2]])) == 0

    def test_get_emergent_cycles_symb(self):
        network1 = [
                ['"x"','"y"'],
                ['"y"','"z"'],
                ['"z"','"w"']
            ]
        sys1 = ReactionSystem(network1)
        assert sys1.get_emergent_cycles_symb([[],[]]).shape[0] == 0
        assert sys1.get_emergent_cycles_symb([[],[0]]).shape[0] == 1
        assert sys1.get_emergent_cycles_symb([[1],[0,1]]).shape[0] == 1
        assert sys1.get_emergent_cycles_symb([[0,1,2],[0,1,2]]).shape[0] == 0


if __name__ == '__main__':
    unittest.main()