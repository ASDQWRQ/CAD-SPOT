"""Ranking contracts and compatibility, not evidence of simulated advantage."""
from contextlib import redirect_stderr
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import cadspot as c
from test_cadspot import spots_fixture


def peak_fixture(isolated=False):
    d = spots_fixture().head(32).copy()
    if isolated:
        d.loc[0, ["x", "y"]] = [1000., 1000.]
    residuals = np.full(len(d), -1.)
    residuals[0] = 31.
    return d, residuals


class RankingMethodTests(unittest.TestCase):
    def test_default_is_positive_product_with_exact_formula(self):
        d = spots_fixture()
        _, a, b = c.fit_pair(d, "target", "reference", ["tech", "comp"])
        r = c.spatial_ranking(d, a, b)
        explicit = c.spatial_ranking(d, a, b, ranking_method="positive_product")
        pd.testing.assert_frame_equal(r, explicit)
        np.testing.assert_array_equal(r.priority_score,
            np.maximum(r.target_resid_z, 0) * np.maximum(r.reference_resid_z, 0))
        expected = r.spatial_status.eq("ok") & (r.target_resid_z > 0) & (r.reference_resid_z > 0)
        pd.testing.assert_series_equal(r.ranking_eligible, expected, check_names=False)
        pd.testing.assert_series_equal(r.cadspot_candidate,
            expected & (r.priority_robust_z >= 3), check_names=False)
        self.assertTrue(r.ranking_method.eq("positive_product").all())

    def test_product_candidates_do_not_require_positive_neighbor_support(self):
        for isolated in (False, True):
            with self.subTest(isolated=isolated):
                d, a = peak_fixture(isolated)
                r = c.spatial_ranking(d, a, a)
                self.assertLessEqual(r.loc[0, "neighbor_reference_resid_z"], 0)
                if isolated:
                    self.assertEqual(r.loc[0, "neighbor_reference_resid_z"], 0)
                self.assertTrue(r.loc[0, "positive_product_eligible"])
                self.assertFalse(r.loc[0, "positive_enrichment_eligible"])
                self.assertTrue(r.loc[0, "ranking_eligible"])
                self.assertTrue(r.loc[0, "cadspot_candidate"])
                self.assertEqual(r.loc[0, "priority_rank"], 1)
                self.assertEqual(r.loc[0, "neighbor_count"], 6)

    def test_legacy_and_positive_support_keep_support_filter(self):
        d, a = peak_fixture()
        for method in ("legacy", "positive_support"):
            with self.subTest(method=method):
                r = c.spatial_ranking(d, a, a, ranking_method=method)
                expected = r.positive_product_eligible & (r.neighbor_reference_resid_z > 0)
                pd.testing.assert_series_equal(r.ranking_eligible, expected, check_names=False)
                pd.testing.assert_series_equal(r.cadspot_candidate,
                    expected & (r.priority_robust_z >= 3), check_names=False)
                self.assertFalse(r.loc[0, "cadspot_candidate"])
                self.assertTrue(np.isnan(r.loc[0, "priority_rank"]))
                if method == "positive_support":
                    np.testing.assert_array_equal(r.priority_score,
                        r.baseline_positive_product * np.maximum(r.neighbor_reference_resid_z, 0))

    def test_legacy_numeric_regression_against_original_release(self):
        d = pd.DataFrame(dict(sample_id=["s"]*8, barcode=[f"b{i}" for i in range(8)],
                              x=np.arange(8) % 4, y=np.arange(8) // 4))
        a = np.array([-3., -1., 2., 4., -2., 3., 1., 5.])
        b = np.array([-2., 4., 3., -1., -3., 1., 2., 5.])
        r = c.spatial_ranking(d, a, b, robust_z_cutoff=0, ranking_method="legacy")
        # Frozen from the unmodified 20261115 source before this follow-up edit.
        np.testing.assert_allclose(r.priority_score, [
            -.5981281624369797, -.04228543349242864, .3266525866007573,
            .7045110169310198, -.7149943802619831, .25490706728654244,
            .12193269936451741, .3876691840932852], rtol=0, atol=1e-14)
        np.testing.assert_allclose(r.priority_robust_z, [
            -2.4677924927526798, -.7238373446463193, .4337047548885724,
            1.619234618850481, -2.8344599538934103, .20860337054704775,
            -.2086033705470477, .6251441743068712], rtol=0, atol=1e-14)
        np.testing.assert_allclose(r.priority_rank,
            [np.nan, np.nan, 1., np.nan, np.nan, np.nan, np.nan, np.nan])
        self.assertEqual(r.cadspot_candidate.tolist(), [False, False, True, False, False, False, False, False])

    def test_candidate_cutoff_is_inclusive_and_does_not_filter_ranks(self):
        d, a = peak_fixture()
        initial = c.spatial_ranking(d, a, a)
        cutoff = initial.loc[0, "priority_robust_z"]
        at = c.spatial_ranking(d, a, a, robust_z_cutoff=cutoff)
        above = c.spatial_ranking(d, a, a, robust_z_cutoff=np.nextafter(cutoff, np.inf))
        self.assertTrue(at.loc[0, "cadspot_candidate"])
        self.assertFalse(above.loc[0, "cadspot_candidate"])
        pd.testing.assert_series_equal(at.priority_rank, above.priority_rank)

    def test_constant_product_uses_finite_fallback_and_no_candidates(self):
        d, _ = peak_fixture()
        r = c.spatial_ranking(d, np.zeros(len(d)), np.zeros(len(d)))
        self.assertTrue(r.priority_score.eq(0).all())
        self.assertTrue(r.priority_robust_z.eq(0).all())
        self.assertFalse(r.cadspot_candidate.any())
        self.assertTrue(r.priority_rank.isna().all())

    def test_priority_scale_metadata_reports_actual_mad_sd_and_constant(self):
        d, a = peak_fixture()
        cases = [
            ("MAD", "legacy", d.target.to_numpy(), d.reference.to_numpy()),
            ("SD", "positive_product", a, a),
            ("constant", "positive_product", np.zeros(len(d)), np.zeros(len(d))),
        ]
        for expected_method, ranking_method, target, reference in cases:
            with self.subTest(scale_method=expected_method):
                r = c.spatial_ranking(d, target, reference, ranking_method=ranking_method)
                scores = r.priority_score
                center = scores.median()
                mad = 1.4826 * float(np.median(np.abs(scores - center)))
                if expected_method == "MAD":
                    self.assertGreater(mad, 0)
                    expected_scale = mad
                elif expected_method == "SD":
                    self.assertEqual(mad, 0)
                    self.assertGreater(scores.std(ddof=0), 0)
                    expected_scale = scores.std(ddof=0)
                else:
                    self.assertEqual(mad, 0)
                    self.assertEqual(scores.std(ddof=0), 0)
                    expected_scale = 1.
                self.assertTrue(r.priority_scale_method.eq(expected_method).all())
                self.assertTrue(r.priority_center.eq(center).all())
                self.assertTrue(r.priority_scale.eq(expected_scale).all())
                np.testing.assert_array_equal(r.priority_robust_z,
                    (r.priority_score - r.priority_center) / r.priority_scale)

    def test_spatial_exclusions_apply_to_every_ranking_method(self):
        d, a = peak_fixture()
        d.loc[0, "x"] = np.nan
        for method in c.RANKING_METHODS:
            with self.subTest(method=method):
                r = c.spatial_ranking(d, a, a, ranking_method=method)
                self.assertEqual(r.loc[0, "spatial_status"], "excluded_missing_coordinates")
                self.assertTrue(np.isnan(r.loc[0, "priority_score"]))
                self.assertTrue(np.isnan(r.loc[0, "priority_rank"]))
                self.assertFalse(r.loc[0, "ranking_eligible"])
                self.assertFalse(r.loc[0, "cadspot_candidate"])
                self.assertTrue(r.loc[0, ["priority_center", "priority_scale", "priority_scale_method"]].isna().all())
                small = c.spatial_ranking(d.head(6), a[:6], a[:6], ranking_method=method)
                self.assertFalse(small.ranking_eligible.any())
                self.assertFalse(small.cadspot_candidate.any())
                self.assertTrue(small.priority_score.isna().all())
                self.assertTrue(small[["priority_center", "priority_scale", "priority_scale_method"]].isna().all().all())

    def test_unknown_api_method_raises_value_error(self):
        d, a = peak_fixture()
        for method in ("unknown", "", "Positive_Product", None):
            with self.subTest(method=method), self.assertRaisesRegex(ValueError, "Unknown ranking_method"):
                c.spatial_ranking(d, a, a, ranking_method=method)

    def test_coordinate_changes_only_change_annotations_for_product(self):
        d, a = peak_fixture()
        isolated, _ = peak_fixture(isolated=True)
        first = c.spatial_ranking(d, a, a)
        second = c.spatial_ranking(isolated, a, a)
        for col in ("priority_score", "priority_robust_z", "priority_rank", "cadspot_candidate",
                    "positive_product_eligible", "ranking_eligible"):
            pd.testing.assert_series_equal(first[col], second[col])
        self.assertFalse(first.neighbor_reference_resid_z.equals(second.neighbor_reference_resid_z))

    def test_support_comparison_remains_a_separate_restricted_universe(self):
        d, a = peak_fixture()
        r = c.spatial_ranking(d, a, a)
        summary = c.ranking_comparison(r, top_k=1, top_fractions=(.05,))
        signed = summary[summary.scope.eq("pooled") & summary.comparator.eq("signed_residual_product")]
        positive = signed[signed.restriction.eq("positive_quadrant")]
        support = signed[signed.restriction.eq("cadspot_positive_support")]
        self.assertTrue(positive.n_comparison_universe.eq(1).all())
        self.assertTrue(positive.overlap_fraction.eq(1).all())
        self.assertTrue(support.n_comparison_universe.eq(0).all())

    def test_archived_eligibility_does_not_override_selected_method(self):
        d, a = peak_fixture()
        d["positive_product_eligible"] = False
        d["positive_enrichment_eligible"] = True
        d["ranking_eligible"] = False
        d["ranking_method"] = "legacy"
        d["priority_center"] = -999.
        d["priority_scale"] = 999.
        d["priority_scale_method"] = "old_value"
        r = c.spatial_ranking(d, a, a)
        self.assertTrue(r.loc[0, "cadspot_candidate"])
        self.assertFalse(r.loc[0, "positive_enrichment_eligible"])
        self.assertTrue(r.archived__positive_enrichment_eligible.all())
        self.assertTrue(r.archived__ranking_method.eq("legacy").all())
        self.assertTrue(r.ranking_method.eq("positive_product").all())
        self.assertTrue(r.archived__priority_center.eq(-999.).all())
        self.assertTrue(r.archived__priority_scale.eq(999.).all())
        self.assertTrue(r.archived__priority_scale_method.eq("old_value").all())
        self.assertTrue(r.priority_center.eq(0).all())
        self.assertTrue(r.priority_scale.eq(r.priority_score.std(ddof=0)).all())
        self.assertTrue(r.priority_scale_method.eq("SD").all())

    def test_cli_defaults_choices_metadata_and_associations(self):
        d = spots_fixture()
        with tempfile.TemporaryDirectory(dir=ROOT, prefix="test_ranking_tmp_") as folder:
            p = Path(folder)
            source = p / "input.tsv"
            d.to_csv(source, sep="\t", index=False)
            common = ["--input", str(source), "--target", "target", "--reference", "reference",
                      "--technical", "tech", "--composition", "comp", "--coordinate-units", "synthetic_units"]
            first_assoc = None
            for method in c.RANKING_METHODS:
                assoc_path = p / f"{method}.tsv"
                rank_path = p / f"{method}_rank.tsv"
                args = common + ["--output", str(assoc_path), "--coordinate-output", str(rank_path)]
                # The product run intentionally omits the CLI flag to test its default.
                if method != "positive_product":
                    args += ["--ranking-method", method]
                self.assertEqual(c.main(args), 0)
                exported = pd.read_csv(rank_path, sep="\t")
                self.assertTrue(exported.ranking_method.eq(method).all())
                self.assertTrue({"priority_center", "priority_scale", "priority_scale_method"} <= set(exported))
                manifest = json.loads(Path(str(assoc_path) + ".manifest.json").read_text())
                self.assertEqual(manifest["cadspot_version"], "2.2.0")
                self.assertEqual(c.__version__, "2.2.0")
                self.assertEqual(manifest["arguments"]["ranking_method"], method)
                assoc = pd.read_csv(assoc_path, sep="\t")
                if first_assoc is None:
                    first_assoc = assoc
                else:
                    pd.testing.assert_frame_equal(first_assoc, assoc)
                _, a, b = c.fit_pair(d, "target", "reference", ["tech", "comp"])
                expected = c.spatial_ranking(d, a, b, ranking_method=method)
                np.testing.assert_allclose(exported.priority_score, expected.priority_score)
                np.testing.assert_array_equal(exported.cadspot_candidate, expected.cadspot_candidate)
            with redirect_stderr(io.StringIO()) as error, self.assertRaises(SystemExit) as raised:
                c.main(common + ["--output", str(p / "invalid.tsv"), "--ranking-method", "unknown"])
            self.assertEqual(raised.exception.code, 2)
            self.assertIn("invalid choice", error.getvalue())
            self.assertFalse((p / "invalid.tsv").exists())


if __name__ == "__main__":
    unittest.main()
