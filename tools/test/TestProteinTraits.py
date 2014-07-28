'''
Created on 20/02/2014

@author: victor
'''
import unittest
import numpy
from tools.proteinTraits import pca, project_points_and_get_lengths_in_axis


class Test(unittest.TestCase):


    def test_pca(self):
        data = numpy.array([[3.7, 1.7], [4.1, 3.8], [4.7, 2.9], [5.2, 2.8], [6.0, 4.0], [6.3, 3.6], [9.7, 6.3], [10.0, 4.9], [11.0, 3.6], [12.5, 6.4]])
        expected_axis = [[ 0.92849112,  0.37135459], [-0.37135459,  0.92849112]]
        axis = pca(data)
        numpy.testing.assert_array_almost_equal(axis, expected_axis,8)

    def test_project(self):
        data = numpy.array([[3.7, 1.7], [4.1, 3.8], [4.7, 2.9], [5.2, 2.8], [6.0, 4.0], [6.3, 3.6], [9.7, 6.3], [10.0, 4.9], [11.0, 3.6], [12.5, 6.4]])
        axis = [[ 0.92849112,  0.37135459], [-0.37135459,  0.92849112]]
        expected_sizes = [9.9160884714759483, 2.98968700380644]
        sizes = project_points_and_get_lengths_in_axis(data,axis)
        numpy.testing.assert_array_almost_equal(sizes, expected_sizes,12)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_pca']
    unittest.main()