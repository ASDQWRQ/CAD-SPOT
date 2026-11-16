import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import decision_audit as audit


class DecisionAuditTests(unittest.TestCase):
    def fixture(self):
        r = np.random.default_rng(27); n = 80
        return pd.DataFrame(dict(sample_id=['s1']*40+['s2']*40, barcode=[f'b{i:02}' for i in range(40)]*2,
            x=np.tile(np.arange(40),2), y=np.zeros(n), target=r.normal(size=n), reference=r.normal(size=n),
            depth=r.normal(size=n), mixture=r.normal(size=n)))

    def test_ties_use_barcode_without_labels(self):
        d = pd.DataFrame({'barcode':['z','a','b']})
        self.assertEqual(audit.budget_order(d, np.ones(3), np.arange(3)).tolist(), [1,2,0])

    def test_budget_validation(self):
        with self.assertRaises(ValueError):
            audit.audit_table(self.fixture(), 'target', 'reference', ['depth'], ['mixture'], budgets=[10])

    def test_membership_and_exclusions(self):
        d = self.fixture(); d.loc[0, 'x'] = np.nan
        out, summary = audit.audit_table(d, 'target', 'reference', ['depth'], ['mixture'])
        self.assertEqual(out.loc[0, 'model_agreement'], 'excluded')
        self.assertFalse(out.loc[0, 'composition_top25'])
        for k in [10,25,50]:
            counts = out[[f'{v}_top{k}' for v in audit.CONDITIONAL]].sum(axis=1)
            self.assertTrue((out[f'agreement_top{k}'].eq('all_three') == counts.eq(3)).all())
        self.assertEqual(len(summary), 2*3*8)

    def test_zero_scores_never_enter_lists(self):
        d = self.fixture(); d['target'] = 0; d['reference'] = 0
        out, _ = audit.audit_table(d, 'target', 'reference', ['depth'], ['mixture'])
        self.assertFalse(out.raw_top25.any())
        self.assertTrue(out.model_agreement.eq('not_selected').all())

    def test_reordering_does_not_change_top_lists(self):
        d = self.fixture(); a, _ = audit.audit_table(d, 'target', 'reference', ['depth'], ['mixture'])
        b, _ = audit.audit_table(d.sample(frac=1, random_state=5), 'target', 'reference', ['depth'], ['mixture'])
        keys = ['sample_id','barcode']
        a=a.set_index(keys).sort_index(); b=b.set_index(keys).sort_index()
        pd.testing.assert_series_equal(a.composition_top25, b.composition_top25)

    def test_complete_export(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory)/'audit'
            audit.export_audit(self.fixture(), prefix, 'target', 'reference', ['depth'], ['mixture'])
            self.assertEqual(len(list(Path(directory).iterdir())), 4)
            manifest = json.loads(Path(str(prefix)+'.manifest.json').read_text())
            self.assertTrue(manifest['positive_priority_required'])
            self.assertIn('not biological validation', Path(str(prefix)+'.report.html').read_text())


if __name__ == '__main__':
    unittest.main()
