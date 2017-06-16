#! /usr/bin/env python

import unittest
import os
import sys
import logging

from sumcoevolity import posterior
from sumcoevolity.test.support import package_paths

_LOG = logging.getLogger(__name__)


class PosteriorSummaryTestCase(unittest.TestCase):

    def setUp(self):
        self.posterior_paths = (
                package_paths.data_path(
                        'simcoevolity-sim-0048-config-state-run-1.log'),
                package_paths.data_path(
                        'simcoevolity-sim-0048-config-state-run-2.log'))

    def test_with_burnin(self):
        ps = posterior.PosteriorSummary(self.posterior_paths, burnin = 11)
        self.assertEqual(ps.number_of_samples, 20)

        self.assertAlmostEqual(ps.parameter_summaries["concentration"]["mean"], 1.41422)
        self.assertAlmostEqual(ps.parameter_summaries["concentration"]["variance"], 0.0)
        self.assertEqual(ps.parameter_summaries["concentration"]["n"], 20)

        self.assertAlmostEqual(ps.parameter_summaries["number_of_events"]["mean"], 2.0)
        self.assertAlmostEqual(ps.parameter_summaries["number_of_events"]["variance"], 0.0)
        self.assertEqual(ps.parameter_summaries["number_of_events"]["n"], 20)

        self.assertAlmostEqual(ps.parameter_summaries["root_height_c1sp1"]["mean"], 0.008030407181410518)
        self.assertAlmostEqual(ps.parameter_summaries["root_height_c1sp1"]["variance"], 1.8561946270985903e-07)
        self.assertEqual(ps.parameter_summaries["root_height_c1sp1"]["n"], 20)

        nevents = ps.get_number_of_events()
        self.assertEqual(nevents, [(2, 1.0)])

        models = ps.get_models()
        self.assertEqual(models, [((0,0,1), 1.0)])

        self.assertTrue(ps.number_of_events_in_credibility_set(2, 0.95))
        self.assertFalse(ps.number_of_events_in_credibility_set(1, 0.95))
        self.assertFalse(ps.number_of_events_in_credibility_set(3, 0.95))

        self.assertTrue(ps.model_in_credibility_set((0, 0, 1), 0.95))
        self.assertFalse(ps.model_in_credibility_set((0, 1, 1), 0.95))

        self.assertEqual(ps.get_number_of_events_credibility_level(2), 0.0)
        self.assertEqual(ps.get_number_of_events_credibility_level(1), 1.0)
        self.assertEqual(ps.get_number_of_events_credibility_level(3), 1.0)

        self.assertEqual(ps.get_model_credibility_level((0, 0, 1)), 0.0)
        self.assertEqual(ps.get_model_credibility_level((0, 1, 1)), 1.0)


class PosteriorSampleTestCase(unittest.TestCase):

    def setUp(self):
        self.posterior_paths = (
                package_paths.data_path(
                        'simcoevolity-sim-0048-config-state-run-1.log'),
                package_paths.data_path(
                        'simcoevolity-sim-0048-config-state-run-2.log'))

    def test_with_burnin(self):
        ps = posterior.PosteriorSample(self.posterior_paths, burnin = 11)
        self.assertEqual(ps.number_of_samples, 20)

        expected_heights = (
                0.00798437475094069052,
                0.00782100759966990709,
                0.00796695487409240441,
                0.00803601657282466963,
                0.0078943048304134858,
                0.00782662502377271148,
                0.00755171798012292499,
                0.00776181733270465881,
                0.00772720362285364559,
                0.00767089531258386553,
                0.00907459399687917355,
                0.00826956209009048482,
                0.00848745864353430725,
                0.00831266083681090479,
                0.00876436416307290085,
                0.00756437101368213037,
                0.00842274661692436234,
                0.00834449429441110058,
                0.0074808667866104973,
                0.00764610728621554787,
        )

        self.assertAlmostEqual(ps.parameter_samples["root_height_c1sp1"], expected_heights)

        self.assertEqual(ps.get_rank("root_height_c1sp1", 0.005), 0.0)
        self.assertEqual(ps.get_rank("root_height_c1sp1", 0.01), 1.0)
        self.assertAlmostEqual(ps.get_rank("root_height_c1sp1", 0.008), 12.0/20.0)

        nevents = ps.get_number_of_events()
        self.assertEqual(nevents, [(2, 1.0)])

        models = ps.get_models()
        self.assertEqual(models, [((0,0,1), 1.0)])

        self.assertTrue(ps.number_of_events_in_credibility_set(2, 0.95))
        self.assertFalse(ps.number_of_events_in_credibility_set(1, 0.95))
        self.assertFalse(ps.number_of_events_in_credibility_set(3, 0.95))

        self.assertTrue(ps.model_in_credibility_set((0, 0, 1), 0.95))
        self.assertFalse(ps.model_in_credibility_set((0, 1, 1), 0.95))

        self.assertEqual(ps.get_number_of_events_credibility_level(2), 0.0)
        self.assertEqual(ps.get_number_of_events_credibility_level(1), 1.0)
        self.assertEqual(ps.get_number_of_events_credibility_level(3), 1.0)

        self.assertEqual(ps.get_model_credibility_level((0, 0, 1)), 0.0)
        self.assertEqual(ps.get_model_credibility_level((0, 1, 1)), 1.0)
        

if __name__ == '__main__':
    unittest.main()

