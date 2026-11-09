"""Input safeguards, window membership and reuse of accepted remote MSA tickets."""
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch
from inference.inputs import variant,window


class InputTests(unittest.TestCase):
    def test_variants_are_applied_to_one_wt_background(self):
        wt='ACDEFGHIKLMNPQRSTVWY'*2
        mutant,text,changes=variant(wt,'A1V_L10A')
        self.assertEqual(mutant[0],'V');self.assertEqual(mutant[9],'A')
        self.assertEqual(sum(a!=b for a,b in zip(wt,mutant)),2)
        self.assertEqual(text,'A1V_L10A')

    def test_reject_wrong_coordinates_duplicate_and_wide_span(self):
        wt='A'*1200
        for mutation in ['K13A','A0V','A1V_A1G','A1A','A1V_A1200V']:
            with self.subTest(mutation=mutation),self.assertRaises(ValueError):variant(wt,mutation)

    def test_window_keeps_all_changed_sites(self):
        for changes in [[(0,'V')],[(3718,'V')],[(500,'V'),(1521,'G')]]:
            left,right=window(3719,changes)
            self.assertEqual(right-left,1022)
            self.assertTrue(all(left<=p<right for p,_ in changes))


class TicketTests(unittest.TestCase):
    def test_resume_uses_existing_ticket_without_submission(self):
        from inference.alignments import search
        sequence='ACDEFGHIKLMNPQRSTVWY'
        buffer=io.BytesIO()
        with tarfile.open(fileobj=buffer,mode='w:gz') as archive:
            payload=('>101\n'+sequence+'\n>hit\n'+sequence+'\n').encode()
            member=tarfile.TarInfo('uniref.a3m');member.size=len(payload)
            archive.addfile(member,io.BytesIO(payload))
        class Response:
            content=buffer.getvalue()
            def json(self):return dict(id='accepted-ticket',status='COMPLETE')
        calls=[]
        def request(session,method,url,**kwargs):
            calls.append((method,url));return Response()
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)
            (folder/'ticket.json').write_text(json.dumps(dict(id='accepted-ticket',status='PENDING')))
            with patch('inference.alignments.request',side_effect=request),patch('inference.alignments.time.sleep'):
                result=search(sequence,folder,'https://example.invalid',10,1)
            self.assertTrue(result.exists())
            self.assertEqual([m for m,_ in calls],['GET','GET'])
            self.assertIn('/ticket/accepted-ticket',calls[0][1])
            self.assertIn('/result/download/accepted-ticket',calls[1][1])
            with patch('inference.alignments.request',side_effect=AssertionError('Unexpected request')):
                self.assertEqual(search(sequence,folder,'https://example.invalid',10,1),result)


if __name__=='__main__':unittest.main()
