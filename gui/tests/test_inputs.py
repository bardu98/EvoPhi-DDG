import unittest
from gui.inputs import single_input,batch_input,fasta_records

SEQUENCE='ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY'


class InputTests(unittest.TestCase):
    def test_multiple_variants_and_multisite_are_distinct(self):
        frame=single_input(SEQUENCE,'A1V\nC2G\nA1V_C2G')
        self.assertEqual(frame.MTS.tolist(),['A1V','C2G','A1V_C2G'])

    def test_batch_fasta_identifier_mapping(self):
        fasta='>protein_A description\n'+SEQUENCE+'\n>protein_B\n'+SEQUENCE[::-1]+'\n'
        csv='ID,MTS\nprotein_B,Y1A\nprotein_A,A1V_C2G\n'
        frame=batch_input(csv,fasta)
        self.assertEqual(frame.Sequence.tolist(),[SEQUENCE[::-1],SEQUENCE])

    def test_wrong_wt_names_and_duplicate_fasta_ids(self):
        with self.assertRaisesRegex(ValueError,'Row 1.*sequence residue 1 is A'):
            batch_input('ID,Sequence,MTS\np,'+SEQUENCE+',K1V\n')
        with self.assertRaisesRegex(ValueError,'absent from the FASTA'):
            batch_input('ID,MTS\nunknown,A1V\n','>protein\n'+SEQUENCE+'\n')
        with self.assertRaisesRegex(ValueError,'more than once'):
            fasta_records('>same\n'+SEQUENCE+'\n>same\n'+SEQUENCE+'\n')


if __name__=='__main__':unittest.main()
