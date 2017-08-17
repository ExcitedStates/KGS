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

    def test_rmsd(self):
        pdb1 = PDBFile("5udi.pdb").models[0]
        pdb2 = PDBFile("5udi_perturbed.pdb").models[0]

        rmsd = pdb1.rmsd(pdb2, names=['CA'])
        self.assertTrue(5.2 < rmsd < 5.4, msg="CA RMSD is "+str(rmsd));
        self.assertTrue(pdb1.rmsd(pdb1) < 0.01, msg="RMSD is "+str(rmsd));

    def test_align(self):
        pdb0 = PDBFile("5udi.pdb").models[0]
        pdb1 = PDBFile("5udi.pdb").models[0]
        pdb2 = PDBFile("5udi_perturbed.pdb").models[0]
        names = ['CA']

        pdb1.alignto(pdb2, names)

        rmsd01 = pdb0.rmsd_cur(pdb1, names)
        rmsd12 = pdb1.rmsd_cur(pdb2, names)
        rmsd02 = pdb0.rmsd_cur(pdb2, names)

        self.assertTrue(1.0 < rmsd01, msg="CA RMSD from 0 to 1 is "+str(rmsd01))
        self.assertTrue(rmsd12 < rmsd02, msg="CA RMSD should decrease")
        self.assertAlmostEqual(rmsd12, pdb0.rmsd(pdb2, names), msg="CA RMSD should match after alignment")

if __name__ == '__main__':
    unittest.main()
