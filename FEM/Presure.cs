using System;
using System.Collections.Generic;
using System.Text;

namespace FEM
{
    static class Presure
    {
        public static double[,] createPSI()
        {
            double[,] PSIET = new double[8, 9];

            double[] values;
            double[][] nodes = Globals.GaussNodes9;

            for (int i = 0; i < 8; i++)
            {
                for (int j = 0; j < 9; j++)
                {
                    values = nodes[j];
                    PSIET[i, j] = PSI.getPsi(i, values[0], values[1]);
                }
            }

            return PSIET;
        }

        public static double[] createF(int nqp, int nel, int m, int n, double[][] AKT, int[][] NT, double[,] PSIET, double[] c)
        {
            double[,,] DXYZET;
            double[] F = new double[3 * nqp];
            double[,,] DPSITE = Globals.DPSITE;
            int[][] PAdapter = new int[6][] {
            new int[8] { 0, 1, 5, 4, 8, 13, 16, 12},
            new int[8] { 1, 2, 6, 5, 9, 14, 17, 13},
            new int[8] { 2, 3, 7, 6, 10, 15, 18, 14 },
            new int[8] { 3, 0, 4, 7, 11, 12, 19, 15},
            new int[8] { 0, 1, 2, 3, 8, 9, 10, 11},
            new int[8] { 4, 5, 6, 7, 16, 17, 18, 19}
        };
            int site = 5;

            int loadElementsCount = m * n;
            int start = nel - loadElementsCount;
            for (int number = start; number < nel; number++)
            {
                DXYZET = new double[3, 2, 9];

                int[] coordinates = NT[number];


                double globalCoordinate = 0;
                double diPsi = 0;
                double sum = 0;


                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        for (int k = 0; k < 9; k++)
                        {
                            sum = 0;
                            for (int l = 0; l < 8; l++)
                            {
                                globalCoordinate = AKT[coordinates[PAdapter[site][l]]][i];
                                diPsi = DPSITE[k, j, l];
                                sum += globalCoordinate * diPsi;
                            }
                            DXYZET[i, j, k] = sum;
                        }
                    }
                }


                double presure = -0.3;

                double[] f2 = new double[8];

                for (int i = 0; i < 8; i++)
                {
                    sum = 0;
                    int counter = 0;
                    for (int k = 0; m < 3; m++)
                    {
                        for (int p = 0; n < 3; n++)
                        {
                            sum += presure *
                                (DXYZET[0, 0, counter] * DXYZET[1, 1, counter] - DXYZET[1, 0, counter] * DXYZET[0, 1, counter]) *
                                PSIET[i, counter]
                                * c[p] * c[k];
                            ++counter;
                        }
                    }
                    f2[i] = sum;
                }

                for (int i = 0; i < 8; i++)
                {
                    F[coordinates[PAdapter[site][i]] * 3 + 2] += f2[i];
                }
            }
            return F;
        }
        
    }
}
