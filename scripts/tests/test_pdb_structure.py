import unittest
from pdb_structure import *


class TestPDBFile(unittest.TestCase):
    def setUp(self):
        pass

    def test_pdbfile_init(self):
        pdb_file = PDBFile("5udi.pdb")
        self.assertEqual(len(pdb_file.models), 1)
        self.assertEqual(len(pdb_file.models[0].atoms), 8387)

        # Check neighbors
        for a in pdb_file.models[0].atoms:
            if a.resn in ["HOH", "CA", "MG"]:
                self.assertEqual(len(a.neighbors), 0, msg="Atom "+str(a)+" should have no neighbors")
            else:
                self.assertGreater(len(a.neighbors), 0, msg="Atom "+str(a)+" should have neighbors")

    def test_clean_waters(self):
        pdb_file = PDBFile("5udi.pdb")

        self.assertTrue([ a for a in pdb_file.models[0].atoms if a.resn == "HOH"])
        near_water = pdb_file.models[0].get_atom("A", 29, "H")
        nearby = pdb_file.models[0].get_nearby(near_water.pos, 1.8)
        self.assertEqual(len(nearby), 3)  # Currently near itself, a water and CA

        pdb_file.clean_waters()
        self.assertEqual(len(pdb_file.models[0].atoms), 7993, msg="Should prune 394 waters")
        self.assertFalse([a for a in pdb_file.models[0].atoms if a.resn == "HOH"])

        nearby = pdb_file.models[0].get_nearby(near_water.pos, 1.8)
        self.assertEqual(len(nearby), 2)  # Now only near itself and CA

    def test_rebuild_atomids(self):
        pdb_file = PDBFile("5udi.pdb")
        first_atom = pdb_file.models[0].get_atom("A", 8, "N")
        last_atom = pdb_file.models[0].get_atom("B", 125, "O")
        self.assertEqual(first_atom.id, 1)
        self.assertEqual(last_atom.id, 8389)

        pdb_file.rebuild_atomids()

        self.assertEqual(first_atom.id, 1)
        self.assertEqual(last_atom.id, 8387)

        pdb_file.clean_waters()
        pdb_file.rebuild_atomids()

        self.assertEqual(first_atom.id, 1)
        self.assertEqual(last_atom.id, 8387)
        last_atom = pdb_file.models[0].atoms[-1]
        self.assertEqual(last_atom.id, 7993)


if __name__ == '__main__':
    unittest.main()
