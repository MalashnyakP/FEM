using System;
using System.Collections.Generic;
using System.Text;

namespace FEM
{
    static class VertexFill
    {
        public static KeyValuePair<double[,,][],int> fillMatrixWithMainVertexes(double[,,][] matrix,int nqp,int k,int n,int m,double SCALE_X, double SCALE_Y, double SCALE_Z)
        {
            for (int deltaZ = 0; deltaZ < k + 1; deltaZ++)
            {
                for (int deltaY = 0; deltaY < n + 1; deltaY++)
                {
                    for (int deltaX = 0; deltaX < m + 1; deltaX++)
                    {
                        matrix[deltaX * 2, deltaY * 2, deltaZ * 2] = new double[] { SCALE_X * deltaX, SCALE_Y * deltaY, SCALE_Z * deltaZ };
                        nqp++;
                    }
                }
            }
            KeyValuePair<double[,,][], int> res = new KeyValuePair<double[,,][], int>(matrix, nqp);

            return res;
        }

        public static KeyValuePair<double[,,][], int> fillMatrixWithIntermidiateVertexes(double[,,][] matrix, int nqp, int k, int n, int m, double SCALE_X, double SCALE_Y, double SCALE_Z)
        {
            double[] current;
            for (int deltaZ = 0; deltaZ < k + 1; deltaZ++)
            {
                for (int deltaY = 0; deltaY < n + 1; deltaY++)
                {
                    for (int deltaX = 0; deltaX < m + 1; deltaX++)
                    {
                        current = matrix[deltaX * 2, deltaY * 2, deltaZ * 2];
                        if (deltaX != m)
                        {
                            matrix[deltaX * 2 + 1, deltaY * 2, deltaZ * 2] = new double[] { current[0] + SCALE_X / 2, current[1], current[2] };
                            nqp++;
                        }

                        if (deltaY != n)
                        {
                            matrix[deltaX * 2, deltaY * 2 + 1, deltaZ * 2] = new double[] { current[0], current[1] + SCALE_Y / 2, current[2] };
                            nqp++;
                        }

                        if (deltaZ != k)
                        {
                            matrix[deltaX * 2, deltaY * 2, deltaZ * 2 + 1] = new double[] { current[0], current[1], current[2] + SCALE_Z / 2 };
                            nqp++;
                        }
                    }
                }
            }
            KeyValuePair<double[,,][], int> res = new KeyValuePair<double[,,][], int>(matrix,nqp);
            
            return res;
        }
    }
}
