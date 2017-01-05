#! /usr/bin/env python

import unittest
import os
import math
import logging
import random

from sumcoevolity.utils.stats import *
from sumcoevolity.test import TestLevel
from sumcoevolity.test.support.sumcoevolity_test_case import SumcoevolityTestCase

_LOG = logging.getLogger(__name__)
GLOBAL_RNG = random.Random()

class SampleSummarizerTestCase(SumcoevolityTestCase):

    def test_init(self):
        ss = SampleSummarizer(tag='test')
        self.assertEqual(ss.tag, 'test')
        self.assertEqual(ss.minimum, None)
        self.assertEqual(ss.maximum, None)
        self.assertEqual(ss.mean, None)
        self.assertEqual(ss.variance, None)
        self.assertEqual(ss.std_deviation, None)
        self.assertEqual(ss.pop_variance, None)

    def test_add_one_sample(self):
        ss = SampleSummarizer(tag='test')
        ss.add_sample(1)
        self.assertEqual(ss.tag, 'test')
        self.assertEqual(ss.minimum, 1)
        self.assertEqual(ss.maximum, 1)
        self.assertApproxEqual(ss.mean, 1.0, 1e-9)
        self.assertEqual(ss.variance, float('inf'))
        self.assertEqual(ss.std_deviation, float('inf'))
        self.assertEqual(ss.pop_variance, 0)

        ss = SampleSummarizer(tag='test')
        ss.add_sample(3.45)
        self.assertEqual(ss.tag, 'test')
        self.assertEqual(ss.minimum, 3.45)
        self.assertEqual(ss.maximum, 3.45)
        self.assertApproxEqual(ss.mean, 3.45, 1e-9)
        self.assertEqual(ss.variance, float('inf'))
        self.assertEqual(ss.std_deviation, float('inf'))
        self.assertEqual(ss.pop_variance, 0)

    def test_update_samples(self):
        ss = SampleSummarizer(tag='test')
        ss.update_samples([1.0, 2.0, 3.0])
        self.assertEqual(ss.tag, 'test')
        self.assertEqual(ss.minimum, 1.0)
        self.assertEqual(ss.maximum, 3.0)
        self.assertApproxEqual(ss.mean, 2.0, 1e-9)
        self.assertApproxEqual(ss.variance, 1.0, 1e-9)
        self.assertEqual(ss.std_deviation, math.sqrt(1.0), 1e-9)
        self.assertApproxEqual(ss.pop_variance, 2/float(3), 1e-9)

    def test_init_with_samples(self):
        ss = SampleSummarizer([1.0, 2.0, 3.0])
        self.assertEqual(ss.minimum, 1.0)
        self.assertEqual(ss.maximum, 3.0)
        self.assertApproxEqual(ss.mean, 2.0, 1e-9)
        self.assertApproxEqual(ss.variance, 1.0, 1e-9)
        self.assertEqual(ss.std_deviation, math.sqrt(1.0), 1e-9)
        self.assertApproxEqual(ss.pop_variance, 2/float(3), 1e-9)

class MedianTestCase(unittest.TestCase):
    def test_empty(self):
        samples = []
        self.assertRaises(ValueError, median, samples)
    
    def test_sample_size_1(self):
        samples = [1.3]
        med = median(samples)
        self.assertEqual(samples[0], med)

    def test_sample_size_even(self):
        samples = [1.1, 1.2, 1.3, 1.4]
        med = median(samples)
        self.assertAlmostEqual(med, 1.25)

    def test_sample_size_odd(self):
        samples = [1.1, 1.2, 1.3, 1.4, 1.5]
        med = median(samples)
        self.assertAlmostEqual(med, 1.3)

class ModeListTestCase(unittest.TestCase):
    def test_empty(self):
        samples = []
        self.assertRaises(ValueError, mode_list, samples)

    def test_ints(self):
        samples = [1,2,3,4,5]
        md = mode_list(samples)
        self.assertEqual(md, samples)

        samples = [1,2,2,3,4,5]
        md = mode_list(samples)
        self.assertEqual(md, [2])
        md = mode_list(samples, bin_width=None)
        self.assertEqual(md, [2])
        md = mode_list(samples, bin_width='a')
        self.assertEqual(md, [2])

        samples = [1,2,2,3,4,5,5]
        md = mode_list(samples)
        self.assertEqual(sorted(md), sorted([2, 5]))

    def test_strings(self):
        samples = ['a', 'b', 'b', 'c', 'd']
        md = mode_list(samples)
        self.assertEqual(md, ['b'])

    def test_floats_no_binning(self):
        samples = [1.1,2.1,2.1,3.1,4.1,5.1]
        md = mode_list(samples, bin_width=None)
        self.assertEqual(md, [2.1])
        md = mode_list(samples, bin_width='auto')
        self.assertNotEqual(md, [2.1])

    def test_floats(self):
        samples = [1.111, 1.112, 1.115, 1.16, 1.121]
        md = mode_list(samples, bin_width = 0.01, zero_value = 'b')
        self.assertEqual(sorted(md), sorted([(1.11, 1.12)]))

class IntervalTestCase(unittest.TestCase):
    def setUp(self):
        self.samples = [GLOBAL_RNG.normalvariate(0, 1) for i in range(100000)]
        self.exp_samples = [GLOBAL_RNG.expovariate(1) for i in range(100000)]

    def test_standard_normal_hpd(self):
        if not TestLevel.test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            return
        hpdi = get_hpd_interval(self.samples, 0.95)
        self.assertAlmostEqual(hpdi[0], -1.96, places=1)
        self.assertAlmostEqual(hpdi[1], 1.96, places=1)

    def test_standard_normal_quantile(self):
        if not TestLevel.test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            return
        quants = quantile_95(self.samples)
        q025 = quantile(self.samples, p=0.025)
        q975 = quantile(self.samples, p=0.975)
        self.assertAlmostEqual(q025, quants[0])
        self.assertAlmostEqual(q975, quants[1])
        self.assertAlmostEqual(quants[0], -1.96, places=1)
        self.assertAlmostEqual(quants[1], 1.96, places=1)

    def test_exp_hpd(self):
        if not TestLevel.test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            return
        hpdi = get_hpd_interval(self.exp_samples, 0.95)
        self.assertAlmostEqual(hpdi[0], 0.0, places=1)
        self.assertAlmostEqual(hpdi[1], 2.9957, places=1)

    def test_exp_quantile(self):
        if not TestLevel.test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            return
        quants = quantile_95(self.exp_samples)
        q025 = quantile(self.exp_samples, p=0.025)
        q975 = quantile(self.exp_samples, p=0.975)
        self.assertAlmostEqual(q025, quants[0])
        self.assertAlmostEqual(q975, quants[1])
        self.assertAlmostEqual(quants[0], 0.0253, places=1)
        self.assertAlmostEqual(quants[1], 3.6889, places=1)

class GetSummaryTestCase(unittest.TestCase):
    def setUp(self):
        self.samples = [GLOBAL_RNG.normalvariate(0, 1) for i in range(100000)]

    def test_standard_normal(self):
        if not TestLevel.test_enabled(
                level = TestLevel.EXHAUSTIVE,
                log = _LOG,
                module_name = '.'.join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            return
        d = get_summary(self.samples)
        self.assertEqual(d['n'], len(self.samples))
        self.assertEqual(d['range'][0], min(self.samples))
        self.assertEqual(d['range'][1], max(self.samples))
        self.assertAlmostEqual(d['mean'], 0.0, places=1)
        self.assertAlmostEqual(d['median'], 0.0, places=1)
        self.assertEqual(len(d['modes'][0]), 2)
        self.assertAlmostEqual(d['modes'][0][0], 0.0, places=0)
        self.assertAlmostEqual(d['modes'][0][1], 0.0, places=0)
        self.assertAlmostEqual(d['variance'], 1.0, places=1)
        self.assertAlmostEqual(d['qi_95'][0], -1.96, places=1)
        self.assertAlmostEqual(d['qi_95'][1], 1.96, places=1)
        self.assertAlmostEqual(d['hpdi_95'][0], -1.96, places=1)
        self.assertAlmostEqual(d['hpdi_95'][1], 1.96, places=1)
        
class GetCountsTestCase(unittest.TestCase):

    def test_get_counts(self):
        x = [0,0,0,1,1,1,1,2,3,4]
        expected = {0: 3, 1: 4, 2: 1, 3: 1, 4: 1}
        counts = get_counts(x)
        self.assertEqual(counts, expected)

class GetFreqsTestCase(unittest.TestCase):

    def test_get_counts(self):
        x = [0,0,0,1,1,1,1,2,3,4]
        expected = {0: 0.3, 1: 0.4, 2: 0.1, 3: 0.1, 4: 0.1}
        freqs = get_freqs(x)
        self.assertAlmostEqual(sum(freqs.values()), 1.0)
        for k, v in freqs.iteritems():
            self.assertAlmostEqual(v, expected[k])

class FreqLessThanTestCase(unittest.TestCase):

    def test_estimate_prob_zero(self):
        x = [0.0045, 0.00021, 0.00012, 0.009999, 0.001, 0.01, 0.010001, 0.9,
                0.09, 1.3]
        self.assertAlmostEqual(freq_less_than(x, 0.01), 0.5)
        self.assertAlmostEqual(freq_less_than(x, 2.0), 1.0)
        self.assertAlmostEqual(freq_less_than(x, 1.3), 0.9)

class MeanSquaredErrorTestCase(unittest.TestCase):
    def test_zero(self):
        x = [-1.0, 2.0, 4.0]
        y = [-1.0, 2.0, 4.0]
        mse = mean_squared_error(x,y)
        self.assertAlmostEqual(mse, 0.0)

    def test_one(self):
        x = [1.0, 2.0, 3.0]
        y = [2.0, 1.0, 4.0]
        mse = mean_squared_error(x,y)
        self.assertAlmostEqual(mse, 1.0)

    def test_simple(self):
        x = [-1.0, 5.5, 10.1, 1016.3]
        y = [-2.0, 8.5, 12.1, 1012.3]
        mse = mean_squared_error(x,y)
        self.assertAlmostEqual(mse, 30/float(4))

class RootMeanSquaredErrorTestCase(unittest.TestCase):
    def test_zero(self):
        x = [-1.0, 2.0, 4.0]
        y = [-1.0, 2.0, 4.0]
        rmse = root_mean_square_error(x,y)
        self.assertAlmostEqual(rmse, 0.0)

    def test_one(self):
        x = [1.0, 2.0, 3.0]
        y = [2.0, 1.0, 4.0]
        rmse = root_mean_square_error(x,y)
        self.assertAlmostEqual(rmse, 1.0)

    def test_simple(self):
        x = [-1.0, 5.5, 10.1, 1016.3]
        y = [-2.0, 8.5, 12.1, 1012.3]
        rmse = root_mean_square_error(x,y)
        self.assertAlmostEqual(rmse, math.sqrt(30/float(4)))
        

if __name__ == '__main__':
    unittest.main()

