import csv
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
import variant_pipeline as p

FIELDS = ['Allele','Consequence','IMPACT','SYMBOL','Feature','CANONICAL','MANE_SELECT','HGVSc','HGVSp','gnomAD_AF','CADD_PHRED','SIFT','PolyPhen','REVEL','ClinVar_CLNSIG']

def ann(allele, feature, canonical='NO', impact='HIGH', af='0.000001', gene='KCNQ2'):
    consequence = {'HIGH':'splice_acceptor_variant','MODERATE':'missense_variant','LOW':'synonymous_variant','MODIFIER':'intron_variant'}[impact]
    values = dict(Allele=allele, Consequence=consequence, IMPACT=impact, SYMBOL=gene,
                  Feature=feature, CANONICAL=canonical, gnomAD_AF=af, CADD_PHRED='25', REVEL='.8')
    return '|'.join(values.get(k,'') for k in FIELDS)

class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root/'panel.txt').write_text('KCNQ2\n')

    def vcf(self, records, fields=FIELDS):
        path = self.root/'input.vcf'
        head = '##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description="VEP Format: '+ '|'.join(fields) +'">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
        path.write_text(head+''.join('\t'.join(map(str,r))+'\n' for r in records))
        return path

    def rec(self, pos, alt, annotations, gt='0/1', ref='A', qual=60):
        return ('1',pos,'.',ref,alt,qual,'PASS','CSQ='+','.join(annotations),'GT:DP:GQ',f'{gt}:40:60')

    def run_case(self, records, model='dominant', overrides=None):
        path = self.vcf(records)
        override_path = None
        if overrides:
            override_path = self.root/'overrides.tsv'
            with override_path.open('w') as fh:
                fh.write('chrom\tpos\tref\talt\ttranscript\n')
                for row in overrides:
                    fh.write('\t'.join(row)+'\n')
        return p.main(['--vcf',str(path),'--panel',str(self.root/'panel.txt'),
                       '--out-prefix',str(self.root/'out'/'candidates'),'--model',model] +
                      (['--transcript-override',str(override_path)] if override_path else []))

    def test_canonical_high_alt_low(self):
        v=self.run_case([self.rec(100,'G',[ann('G','ENST1','YES','HIGH'),ann('G','ENST2','NO','LOW')])])[0]
        self.assertEqual((v.status,v.csq['Feature']),('candidate','ENST1'))
        self.assertEqual(len(v.all_csq),2)

    def test_canonical_low_alt_high_is_review_not_lost(self):
        v=self.run_case([self.rec(101,'G',[ann('G','ENST1','YES','LOW'),ann('G','ENST2','NO','HIGH')])])[0]
        self.assertEqual(v.status,'candidate')
        self.assertIsNotNone(v.score)
        with open(self.root/'out'/'candidates.tsv') as fh:
            rows=list(csv.DictReader(fh,delimiter='\t'))
        self.assertEqual(len(rows),1)
        self.assertIn('HIGH', rows[0]['all_csq_json'])
        self.assertIn('LOW', rows[0]['all_csq_json'])
        self.assertIn('101', (self.root/'out'/'candidates.report.txt').read_text())

    def test_override_reprocesses_high(self):
        v=self.run_case([self.rec(101,'G',[ann('G','ENST1','YES','LOW'),ann('G','ENST2','NO','HIGH')])],
                        overrides=[('1','101','A','G','ENST2')])[0]
        self.assertEqual((v.status,v.csq['Feature'],v.transcript_status),('candidate','ENST2','manual'))
        self.assertIsNotNone(v.score)

    def test_missing_canonical_and_invalid_override(self):
        record=self.rec(102,'G',[ann('G','ENST1','NO','HIGH')])
        self.assertEqual(self.run_case([record])[0].status,'candidate')
        v=self.run_case([record],overrides=[('1','102','A','G','ENST_BAD')])[0]
        self.assertEqual((v.status,v.transcript_status),('candidate','manual_not_found'))

    def test_two_alts_no_cross_contamination_and_genotype(self):
        r=self.rec(103,'G,T',[ann('G','ENSTG','YES','LOW'),ann('T','ENSTT','YES','HIGH')],gt='0/2')
        vs=self.run_case([r]); self.assertEqual(len(vs),2)
        self.assertEqual((vs[0].status,vs[1].status),('excluded_consequence','candidate'))
        self.assertEqual([x['Feature'] for x in vs[1].all_csq],['ENSTT'])

    def test_recessive_and_af_cutoff(self):
        r=self.rec(104,'G',[ann('G','ENST1','YES','HIGH',af='.04')],gt='1/1')
        self.assertEqual(self.run_case([r],model='recessive')[0].status,'additional_candidate')
        self.assertEqual(self.run_case([r],model='dominant')[0].status,'additional_candidate')
        r2=self.rec(105,'G',[ann('G','ENST1','YES','HIGH',af='.05')],gt='1/1')
        self.assertEqual(self.run_case([r2],model='recessive')[0].status,'excluded_frequency')

    def test_inheritance_never_excludes_rare_alleles_and_outputs_are_separate(self):
        records = [self.rec(201,'G',[ann('G','ENST1',af='0.000001')],gt='1/1'),
                   self.rec(202,'G',[ann('G','ENST2',af='0.000001')],gt='./.'),
                   self.rec(203,'G',[ann('G','ENST3',af='0.02')],gt='0/1')]
        vs=self.run_case(records,model='dominant')
        self.assertEqual([v.status for v in vs],['candidate','candidate','additional_candidate'])
        self.assertEqual(vs[0].inheritance_compatibility,'other_dosage_review')
        self.assertEqual(vs[1].inheritance_compatibility,'unknown_genotype')
        with open(self.root/'out'/'candidates.tsv') as fh:
            candidates=list(csv.DictReader(fh,delimiter='\t'))
        with open(self.root/'out'/'candidates.excluded.tsv') as fh:
            excluded=list(csv.DictReader(fh,delimiter='\t'))
        self.assertEqual([r['pos'] for r in candidates],['201','202'])
        self.assertEqual([r['pos'] for r in excluded],[])
        with open(self.root/'out'/'candidates.additional_candidates.tsv') as fh:
            additional=list(csv.DictReader(fh,delimiter='\t'))
        self.assertEqual([r['pos'] for r in additional],['203'])
        self.assertIn('ENST1',candidates[0]['all_csq_json'])
        self.assertEqual(candidates[0]['alt_dosage'],'2')

    def test_recessive_mode_does_not_exclude_heterozygous_allele(self):
        v=self.run_case([self.rec(204,'G',[ann('G','ENST1',af='0.000001')],gt='0/1')],model='recessive')[0]
        self.assertEqual(v.status,'candidate')
        self.assertIn('compound_het_not_assessed',v.inheritance_compatibility)

    def test_unknown_conflicting_and_invalid_af(self):
        for af,status in [('', 'candidate'),('nan','review'),('1.3','review')]:
            r=self.rec(106,'G',[ann('G','ENST1','YES','HIGH',af=af)])
            self.assertEqual(self.run_case([r])[0].status,status)
        r=self.rec(107,'G',[ann('G','ENST1','YES','HIGH',af='0'),ann('G','ENST2','NO','HIGH',af='.02')])
        self.assertEqual(self.run_case([r])[0].status,'review')

    def test_header_order_and_missing_field(self):
        fields=list(reversed(FIELDS))
        values=dict(zip(FIELDS,ann('G','ENST1','YES','HIGH').split('|')))
        reversed_ann='|'.join(values.get(k,'') for k in fields)
        path=self.vcf([self.rec(108,'G',[reversed_ann])],fields)
        self.assertEqual(p.load_vcf(path)[0].csq['Feature'],'ENST1')
        bad=self.vcf([self.rec(108,'G',[reversed_ann])],fields[:-1])
        # The removed field is Allele, which is essential for ALT matching.
        with self.assertRaisesRegex(ValueError,'missing required'):
            p.load_vcf(bad)

    def test_missing_canonical_feature_still_filters(self):
        fields=[f for f in FIELDS if f not in ('CANONICAL','Feature')]
        values=dict(zip(FIELDS,ann('G','ENST1','YES','HIGH').split('|')))
        reduced='|'.join(values.get(k,'') for k in fields)
        path=self.vcf([self.rec(111,'G',[reduced])],fields)
        vs=p.main(['--vcf',str(path),'--panel',str(self.root/'panel.txt'),
                   '--out-prefix',str(self.root/'out'/'missing')])
        self.assertEqual(vs[0].status,'candidate')
        self.assertEqual(vs[0].transcript_status,'all_transcripts')

    def test_all_transcripts_scored_once(self):
        r=self.rec(112,'G',[ann('G','ENST1','YES','HIGH'),ann('G','ENST2','NO','HIGH')])
        one=self.run_case([self.rec(113,'G',[ann('G','ENST1','YES','HIGH')])])[0]
        two=self.run_case([r])[0]
        self.assertEqual(one.score,two.score)
        self.assertEqual(two.status,'candidate')

    def test_indel_allele_and_output_accounting(self):
        r=self.rec(109,'AT',[ann('T','ENST1','YES','HIGH')],ref='A')
        v=self.run_case([r])[0]
        self.assertEqual(v.status,'candidate')
        report=(self.root/'out'/'candidates.report.txt').read_text()
        self.assertIn('ALT alleles processed: 1',report)
        self.assertIn('Primary candidates: 1',report)

    def test_parallel_af_lanes_independent_of_model_and_dosage(self):
        records = [self.rec(301,'G',[ann('G','E1',af='0.000001')],gt='1/1'),
                   self.rec(302,'G',[ann('G','E2',af='0.001')],gt='0/1'),
                   self.rec(303,'G',[ann('G','E3',af='0.049999')],gt='1/1'),
                   self.rec(304,'G',[ann('G','E4',af='0.05')],gt='0/1'),
                   self.rec(305,'G',[ann('G','E5',af='0.00001')],gt='0/1')]
        for model in ('dominant','recessive'):
            vs=self.run_case(records,model=model)
            self.assertEqual([v.status for v in vs],
                ['candidate','additional_candidate','additional_candidate',
                 'excluded_frequency','additional_candidate'])
            out=self.root/'out'
            def positions(name):
                with open(out/name) as fh:
                    return [r['pos'] for r in csv.DictReader(fh,delimiter='\t')]
            self.assertEqual(positions('candidates.tsv'),['301'])
            self.assertEqual(set(positions('candidates.additional_candidates.tsv')),
                             {'302','303','305'})
            self.assertEqual(positions('candidates.excluded.tsv'),['304'])
            self.assertEqual(len({v.key for v in vs}),5)
            report=(out/'candidates.report.txt').read_text()
            self.assertIn('Primary candidates: 1',report)
            self.assertIn('Additional candidates: 3',report)

    def test_unknown_af_visible_and_not_falsely_called_zero(self):
        v=self.run_case([self.rec(306,'G',[ann('G','E6',af='')])])[0]
        self.assertEqual((v.status,v.frequency_status,v.af),('candidate','UNKNOWN',None))
        with open(self.root/'out'/'candidates.tsv') as fh:
            row=next(csv.DictReader(fh,delimiter='\t'))
        self.assertEqual(row['frequency_status'],'UNKNOWN')
        self.assertEqual(row['gnomAD_AF'],'')

    def test_unmatched_csq_is_review(self):
        r=self.rec(110,'G',[ann('T','ENST1','YES','HIGH')])
        self.assertEqual(self.run_case([r])[0].status,'review')

if __name__=='__main__':
    unittest.main(verbosity=2)
