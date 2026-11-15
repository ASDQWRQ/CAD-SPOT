"""Focused regression tests; no network, installations, or full-size bootstrap."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np
import pandas as pd
from scipy import sparse, stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import cadspot as c
import cadspot_null as n


def spots_fixture():
    rng = np.random.default_rng(472)
    rows = []
    for sec, shift in [("s1", 0), ("s2", 4)]:
        for i in range(32):
            cov = rng.normal()
            target = 2*cov + shift + rng.normal()
            reference = 1.5*cov - shift + .4*target + rng.normal()
            rows.append(dict(sample_id=sec, barcode=f"b{i:02}", x=float(i % 8), y=float(i // 8),
                             target=target, reference=reference, tech=cov, comp=rng.uniform()))
    return pd.DataFrame(rows)


def expression_fixture():
    spots = spots_fixture()
    rng = np.random.default_rng(20)
    genes = ["T1", "T2", "R1", "R2", "C1", "C2", "C3", "C4", "EX"]
    values = rng.lognormal(0, .7, (len(genes), len(spots)))
    expression = n.ExpressionData(sparse.csr_matrix(values), genes, spots)
    meta = pd.DataFrame(dict(gene=genes, mean_expression=np.ones(len(genes)),
                             detection_rate=np.ones(len(genes)), gene_length=np.full(len(genes), 1000)))
    return expression, meta


class TableAndEstimatorTests(unittest.TestCase):
    def test_validation_rejects_missing_target_numeric_and_duplicate_ids(self):
        d = spots_fixture()
        for col, value in [("target", np.nan), ("target", "bad"), ("comp", np.inf), ("barcode", "")]:
            bad = d.copy()
            if isinstance(value, str):
                bad[col] = bad[col].astype(object)
            bad.loc[0, col] = value
            with self.subTest(col=col, value=value), self.assertRaises(ValueError):
                c.validate_table(bad, "target", "reference", ["tech"], ["comp"])
        with self.assertRaises(ValueError):
            c.validate_table(pd.concat([d, d.iloc[[0]]]), "target", "reference", ["tech"])
        # The same barcode in different sections is legitimate.
        self.assertEqual(len(c.validate_table(d, "target", "reference", ["tech"])), 64)

    def test_primary_score_residual_and_partial_are_distinct_correct_estimators(self):
        d = spots_fixture()
        rho, a, b = c.fit_pair(d, "target", "reference", ["tech", "comp"])
        self.assertAlmostEqual(rho, stats.spearmanr(a, b).statistic, places=14)
        ranked = d.copy()
        for col in ["target", "reference", "tech", "comp"]:
            ranked[col] = stats.rankdata(ranked[col])
        X = c.design(ranked, ["tech", "comp"])
        manual = c.resid(ranked[["target", "reference"]].to_numpy(), X)
        partial = c.fit_pair(d, "target", "reference", ["tech", "comp"], method="partial_spearman")[0]
        self.assertAlmostEqual(partial, np.corrcoef(manual.T)[0, 1], places=13)
        self.assertNotAlmostEqual(rho, partial, places=4)

    def test_spline_skewed_ties_do_not_produce_duplicate_knot_failure(self):
        d = spots_fixture()
        d["tech"] = 0.0
        d.loc[d.index[-2:], "tech"] = [1., 2.]
        rho = c.fit_pair(d, "target", "reference", ["tech", "comp"], method=c.METHODS[2])[0]
        self.assertTrue(np.isfinite(rho))


class SpatialTests(unittest.TestCase):
    def test_missing_coordinates_and_insufficient_neighbors_explicit(self):
        d = spots_fixture()
        d.loc[0, "x"] = np.nan
        d = pd.concat([d[d.sample_id.eq("s1")], d[d.sample_id.eq("s2")].head(6)], ignore_index=True)
        d = c.validate_table(d, "target", "reference", ["tech"], spatial=True)
        r = c.spatial_ranking(d, d.target.to_numpy(), d.reference.to_numpy())
        self.assertEqual(r.loc[0, "spatial_status"], "excluded_missing_coordinates")
        self.assertTrue(np.isnan(r.loc[0, "priority_score"]))
        self.assertFalse(r.loc[0, "cadspot_candidate"])
        self.assertTrue(r.loc[r.sample_id.eq("s2"), "spatial_status"].eq("excluded_insufficient_neighbors").all())
        for _, row in r[r.spatial_status.eq("ok")].iterrows():
            neighbors = json.loads(row.neighbor_ids)
            self.assertEqual(len(neighbors), 6)
            self.assertEqual(len(set(neighbors)), 6)
            self.assertNotIn(row.barcode, neighbors)
            self.assertNotIn("b00", neighbors)

    def test_exact_six_neighbors_in_seven_spot_section_with_duplicate_positions(self):
        d = spots_fixture().head(7).copy()
        d[["x", "y"]] = 0.0
        r = c.spatial_ranking(d, d.target.to_numpy(), d.reference.to_numpy())
        self.assertTrue(r.neighbor_count.eq(6).all())
        self.assertTrue(r.spatial_status.eq("ok").all())

    def test_neighbors_never_cross_sections_and_ties_deterministic(self):
        d = spots_fixture()
        d["barcode"] = d.sample_id + "_" + d.barcode
        _, a, b = c.fit_pair(d, "target", "reference", ["tech"])
        r = c.spatial_ranking(d, a, b)
        for section, sub in r.groupby("sample_id"):
            local = set(sub.barcode)
            for encoded in sub.neighbor_ids:
                self.assertTrue(set(json.loads(encoded)) <= local)
        d["a"], d["b"] = a, b
        shuffled = d.sample(frac=1, random_state=10).reset_index(drop=True)
        rr = c.spatial_ranking(shuffled, shuffled.a.to_numpy(), shuffled.b.to_numpy())
        left = r.set_index(["sample_id", "barcode"]).neighbor_ids.sort_index()
        right = rr.set_index(["sample_id", "barcode"]).neighbor_ids.sort_index()
        pd.testing.assert_series_equal(left, right)

    def test_core_cli_final_outputs(self):
        d = spots_fixture()
        d["cadspot_candidate"] = True
        with tempfile.TemporaryDirectory(dir=ROOT, prefix="test_tmp_") as folder:
            p = Path(folder)
            d.to_csv(p/"input.tsv", sep="\t", index=False)
            self.assertEqual(c.main(["--input", str(p/"input.tsv"), "--target", "target", "--reference", "reference",
                "--technical", "tech", "--composition", "comp", "--coordinate-units", "synthetic_units",
                "--output", str(p/"association.tsv"), "--coordinate-output", str(p/"rank.tsv"),
                "--ranking-comparison-output", str(p/"baseline.tsv"), "--archived-comparison-output", str(p/"archive.tsv")]), 0)
            assoc = pd.read_csv(p/"association.tsv", sep="\t")
            self.assertEqual(set(assoc.method), set(c.METHODS))
            baseline = pd.read_csv(p/"baseline.tsv", sep="\t")
            exact = baseline[baseline.comparator.eq("signed_residual_product")]
            self.assertTrue(exact.n_cadspot_top.eq(exact.n_baseline_top).all())
            self.assertEqual(set(exact.scope), {"pooled", "section"})
            archive = pd.read_csv(p/"archive.tsv", sep="\t")
            self.assertEqual(archive.iloc[0].n_archived_candidates, len(d))
            ranked = pd.read_csv(p/"rank.tsv", sep="\t")
            np.testing.assert_allclose(ranked.baseline_signed_product, ranked.target_resid_z*ranked.reference_resid_z)

    def test_negative_negative_product_not_joint_high_and_requested_baseline_preserved(self):
        self.assertEqual(float(c.positive_product(-4., -3.)), 0.)
        self.assertEqual((-4.)*(-3.), 12.)
        d = spots_fixture()
        d["cadspot_candidate"] = True
        r = c.spatial_ranking(d, d.target.to_numpy(), d.reference.to_numpy(),
                              target_scores=d.target, reference_scores=d.reference)
        neg = (r.target_resid_z < 0) & (r.reference_resid_z < 0)
        self.assertFalse(r.loc[neg, "cadspot_candidate"].any())
        self.assertTrue(r.loc[neg, "priority_rank"].isna().all())
        self.assertTrue(r.archived__cadspot_candidate.all())
        baseline = c.ranking_comparison(r, top_k=10, top_fractions=(.1,))
        requested = baseline[baseline.comparator.eq("raw_score_z_signed_product") & baseline.restriction.eq("all_spots")]
        self.assertEqual(set(requested.scope), {"pooled", "section"})
        self.assertEqual(set(requested.selection), {"top_k", "top_fraction"})
        self.assertTrue((requested.negative_negative_top_count > 0).any())
        self.assertIn("overlap_fraction_of_baseline", baseline)
        exact = baseline[baseline.comparator.eq("signed_residual_product")]
        self.assertEqual(set(exact.restriction), {"all_spots", "positive_quadrant", "cadspot_positive_support"})
        self.assertTrue(exact.n_cadspot_top.eq(exact.n_baseline_top).all())
        pooled_positive = exact[exact.scope.eq("pooled") & exact.restriction.eq("positive_quadrant")].iloc[0]
        u = r[(r.target_resid_z > 0) & (r.reference_resid_z > 0) & r.spatial_status.eq("ok")]
        self.assertAlmostEqual(pooled_positive.score_rank_spearman, c.correlation(u.priority_score, u.baseline_signed_product), places=14)
        archived = c.archived_candidate_comparison(r)
        self.assertEqual(archived.iloc[0].n_archived_candidates, len(d))


class BootstrapTests(unittest.TestCase):
    def test_refits_and_retains_fixed_section_block_strata(self):
        d = spots_fixture()
        d.loc[0, "x"] = np.nan
        with mock.patch.object(c, "fit_pair", wraps=c.fit_pair) as fitted:
            ci, draws, audit, blocks = c.block_bootstrap(d, "target", "reference", ["tech"], ["comp"],
                                                        block_size=2., n_boot=5, seed=12)
        # Point plus each draw: two adjusted layers, pooled and both sections.
        self.assertEqual(fitted.call_count, (5+1)*2*3)
        self.assertEqual(set(audit.section), {"s1", "s2"})
        self.assertTrue(audit.n_blocks.eq(audit.n_drawn_blocks).all())
        self.assertEqual(audit.groupby("replicate").section.nunique().tolist(), [2]*5)
        self.assertTrue(ci.residualization.eq("refitted_in_each_resample").all())
        self.assertTrue(ci.n_missing_coordinate_excluded.eq(1).all())
        self.assertEqual(set(draws.scope), {"pooled", "section"})
        self.assertTrue((audit.n_unique_drawn_blocks < audit.n_blocks).any())
        # A cached residual bootstrap would not invoke fit_pair on changed rows.
        pooled_inputs = [call.args[0] for call in fitted.call_args_list if call.args[0].sample_id.nunique() == 2]
        self.assertTrue(any(not x.barcode.reset_index(drop=True).equals(blocks.barcode.reset_index(drop=True)) for x in pooled_inputs[2:]))
        repeat = c.block_bootstrap(d, "target", "reference", ["tech"], ["comp"], block_size=2., n_boot=5, seed=12)[1]
        pd.testing.assert_frame_equal(draws, repeat)

    def test_one_block_or_missing_whole_section_fails(self):
        d = spots_fixture()
        with self.assertRaisesRegex(ValueError, "fewer than two"):
            c.make_spatial_blocks(d, 100.)
        d.loc[d.sample_id.eq("s2"), "x"] = np.nan
        with self.assertRaisesRegex(ValueError, "no finite coordinates"):
            c.make_spatial_blocks(d, 2.)


class NullTests(unittest.TestCase):
    def test_strict_matching_unique_seeded_exclusions_and_infeasible_failure(self):
        expr, meta = expression_fixture()
        kw = dict(expression=expr, metadata=meta, target_genes=["T1", "T2"], reference_genes=["R1", "R2"],
                  universe=expr.genes, exclusions=["EX"], n_sets=6, seed=19,
                  mean_tolerance=0., detection_tolerance=0., length_tolerance=0.)
        result = n.generate_matched_sets(**kw)
        pd.testing.assert_frame_equal(result, n.generate_matched_sets(**kw))
        self.assertTrue(result.groupby("set_id").gene.nunique().eq(2).all())
        self.assertFalse(set(result.gene) & {"T1", "T2", "R1", "R2", "EX"})
        bad_meta = meta.copy()
        bad_meta.loc[bad_meta.gene.str.startswith("C"), "gene_length"] = 2000
        with self.assertRaisesRegex(ValueError, "not relaxed"):
            n.generate_matched_sets(**{**kw, "metadata": bad_meta})
        # Every slot has candidates, but no unique full assignment exists.
        hall_meta = meta.copy()
        hall_meta.loc[hall_meta.gene.eq("C2"), "gene_length"] = 2000
        with self.assertRaisesRegex(ValueError, "unique-gene full matching"):
            n.generate_matched_sets(**{**kw, "metadata": hall_meta, "exclusions": ["EX", "C3", "C4"]})

    def test_fresh_observed_and_null_use_identical_score_and_residual_model(self):
        expr, meta = expression_fixture()
        target, ref = ["T1", "T2"], ["R1", "R2"]
        matched = n.generate_matched_sets(expr, meta, target, ref, expr.genes, ["EX"], 4, 7)
        spots = expr.spots.copy()
        spots["target_score"] = 999.  # Archived values cannot enter fresh observed statistic.
        for method in c.METHODS:
            with self.subTest(method=method):
                summary, null, scored = n.calibrate_gene_sets(expr, spots, target, ref, matched,
                    ["tech"], ["comp"], method=method)
                observed = c.fit_pair(scored, "fresh_target_score", "fresh_reference_score", ["tech", "comp"], method=method)[0]
                self.assertAlmostEqual(summary.iloc[0].observed_rho, observed, places=14)
                group = matched[matched.set_id.eq(null.iloc[0].set_id)]
                manual = scored.copy()
                manual["fresh_target_score"] = n.score_gene_set(expr, group.gene.tolist())
                expected = c.fit_pair(manual, "fresh_target_score", "fresh_reference_score", ["tech", "comp"], method=method)[0]
                self.assertAlmostEqual(null.iloc[0].rho, expected, places=14)
                self.assertTrue(summary.iloc[0].observed_random_model_equivalent)

    def test_gene_audit_missing_gene_failure_and_no_sample_mixing(self):
        expr, _ = expression_fixture()
        audit = n.audit_gene_sets(expr, ["T1"], ["R1"], ["T1", "ABSENT"], ["R1"])
        self.assertEqual(audit["target"]["missing_from_expression"], ["ABSENT"])
        with self.assertRaisesRegex(ValueError, "no zero imputation"):
            n.score_gene_set(expr, ["ABSENT"])
        shuffled = expr.spots.sample(frac=1, random_state=22).reset_index(drop=True)
        aligned = expr.aligned(shuffled)
        score = n.score_gene_set(expr, ["T1"])
        expected = pd.Series(score, index=n._spot_keys(expr.spots, ["sample_id", "barcode"]))
        np.testing.assert_allclose(n.score_gene_set(aligned, ["T1"]), expected.reindex(n._spot_keys(shuffled, ["sample_id", "barcode"])).to_numpy())
        with self.assertRaisesRegex(ValueError, "do not join samples"):
            expr.aligned(expr.spots, ["barcode"])
        with self.assertRaisesRegex(ValueError, "exactly the same"):
            expr.aligned(expr.spots.iloc[:-1])

    def test_empirical_p_fails_on_nan_instead_of_dropping_controls(self):
        self.assertEqual(c.empirical_calibration(.8, [.1, .2, .9])["plus_one_p"], .5)
        with self.assertRaises(ValueError):
            c.empirical_calibration(.8, [.1, np.nan])

    def test_sparse_cli_generation_and_calibration_end_to_end(self):
        expr, meta = expression_fixture()
        with tempfile.TemporaryDirectory(dir=ROOT, prefix="test_tmp_") as folder:
            p = Path(folder)
            sparse.save_npz(p/"x.npz", expr.values.T)
            expr.spots.to_csv(p/"spots.tsv", sep="\t", index=False)
            meta.to_csv(p/"meta.tsv", sep="\t", index=False)
            for name, genes in [("genes", expr.genes), ("universe", expr.genes),
                                ("target", ["T1", "T2"]), ("reference", ["R1", "R2"])]:
                pd.DataFrame({"gene": genes}).to_csv(p/f"{name}.tsv", sep="\t", index=False)
            common = ["--expression", str(p/"x.npz"), "--genes", str(p/"genes.tsv"),
                      "--expression-spots", str(p/"spots.tsv"), "--matrix-orientation", "spots-by-genes",
                      "--expression-scale", "synthetic_log_expression", "--target-genes", str(p/"target.tsv"),
                      "--reference-genes", str(p/"reference.tsv")]
            self.assertEqual(n.main(["generate"]+common+["--metadata", str(p/"meta.tsv"), "--universe", str(p/"universe.tsv"),
                "--mean-tolerance", "0", "--detection-tolerance", "0", "--length-tolerance", "0", "--seed", "5", "--n-sets", "3", "--output", str(p/"sets.tsv")]), 0)
            # Compatibility with early generation code's hardcoded composite key.
            generation_path = p/"sets.tsv.manifest.json"
            generation = json.loads(generation_path.read_text())
            del generation["arguments"]["spot_key_cols"]
            generation_path.write_text(json.dumps(generation))
            self.assertEqual(n.main(["calibrate"]+common+["--spot-table", str(p/"spots.tsv"), "--matched-sets", str(p/"sets.tsv"),
                "--technical", "tech", "--composition", "comp", "--null-output", str(p/"rho.tsv"), "--output", str(p/"summary.tsv")]), 0)
            self.assertEqual(len(pd.read_csv(p/"rho.tsv", sep="\t")), 3)
            manifest = json.loads((p/"summary.tsv.manifest.json").read_text())
            self.assertTrue(manifest["observed_random_model_equivalent"])
            self.assertFalse(manifest["spatial_inference"])
            self.assertTrue(manifest["generation_spot_key_default_inferred"])


if __name__ == "__main__":
    unittest.main()
