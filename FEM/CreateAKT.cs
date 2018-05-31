using System.Collections.Generic;
using System.Linq;

using System;
using System.Collections.Generic;
using System.Text;

namespace FEM
{
    static class CreateAKT
    {
        public static double[][] createAKT(int nqp,int k, int n,int m,double[,,][] matrix,int[,,] local_to_global)
        {
            double[][] AKT = new double[nqp][];
            int counter = 0;
            for (int deltaZ = 0; deltaZ < k * 2 + 1; deltaZ++)
            {
                for (int deltaY = 0; deltaY < n * 2 + 1; deltaY++)
                {
                    for (int deltaX = 0; deltaX < m * 2 + 1; deltaX++)
                    {
                        if (matrix[deltaX, deltaY, deltaZ] != null)
                        {
                            AKT[counter] = matrix[deltaX, deltaY, deltaZ];
                            local_to_global[deltaX, deltaY, deltaZ] = counter;
                            counter++;
                        }
                    }
                }
            }
            return AKT;
        }

        public static int[] createZU(double[][] AKT)
        {
            int i = 0;
            while (AKT[i][2] == 0)
            {
                i++;
            }
            int[] ZU = Enumerable.Range(0, i).ToArray();
            return ZU;
        }

        public static double[,] createZP(int m,int n,int nel)
        {
            int loadElementsCount = m * n;
            double[,] ZP = new double[loadElementsCount, 3];
            int firstOne = nel - loadElementsCount;
            for (int i = firstOne, counter = 0; i < nel; i++, counter++)
            {
                ZP[counter, 0] = i;
                ZP[counter, 1] = 5;
                ZP[counter, 2] = 10;
            }
            return ZP;
        }
    }
}
