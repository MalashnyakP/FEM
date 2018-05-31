using System;
using System.Collections.Generic;
using System.Text;

namespace FEM
{
    class SquareGenerator
    {
        private double[] value;
        private double[,,][] matrix;
        private Dictionary<int, int[]> adapter = Globals.magicDictionary;

        public SquareGenerator(double lower, double middle, double higher)
        {
            value = new double[] { lower, middle, higher };
        }

        public double[,,][] getMatrix()
        {
            matrix = new double[3, 3, 3][];
            for (int z = 0; z < 3; z++)
            {
                for (int y = 0; y < 3; y++)
                {
                    for (int x = 0; x < 3; x++)
                    {
                        matrix[x, y, z] = new double[] { value[x], value[2 - y], value[z] };                    }
                }
            }

            return matrix;
        }

        public double[][] getItInSingleDArray()
        {
            getMatrix();
            double[][] oneline = new double[20][];
            int[] point;
            for (int i = 0; i < 20; i++)
            {
                point = adapter[i];
                oneline[i] = matrix[point[0], point[1], point[2]];
            }
            return oneline;
        }
    }
}

