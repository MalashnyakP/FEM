using System;
using System.Collections.Generic;
using System.Text;

namespace FEM
{
    static class createMG
    {
        public static double[,] getMG(int nel,int[][] NT,double[][] AKT,double[,,] DFIABG, double lam, double v, double mu, double[] c,int nqp)
        {
            double[,,] dxyzabg = new double[3, 3, 27];
            double[] dj = new double[27];
            double[,,] dfixyz = new double[27, 20, 3];
            double[,] MG = new double[3 * nqp, 3 * nqp];

            for (int number = 0; number < nel; number++)
            {
                int[] coordinates = NT[number];

                double globalCoordinate = 0;
                double diFi = 0;
                double sum = 0;

                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        for (int k = 0; k < 27; k++)
                        {
                            sum = 0;
                            for (int l = 0; l < 20; l++)
                            {
                                globalCoordinate = AKT[coordinates[l]][i];
                                diFi = DFIABG[k, j, l];
                                sum += globalCoordinate * diFi;
                            }
                            dxyzabg[i, j, k] = sum;
                        }
                    }
                }

                double[,] jak;
                for (int i = 0; i < 27; i++)
                {
                    jak = new double[3, 3] {
                    { dxyzabg[0,0,i], dxyzabg[1,0,i], dxyzabg[2,0,i] },
                    { dxyzabg[0,1,i], dxyzabg[1,1,i], dxyzabg[2,1,i] },
                    { dxyzabg[0,2,i], dxyzabg[1,2,i], dxyzabg[2,2,i] }
                };
                    dj[i] = (
                                jak[0, 0] * jak[1, 1] * jak[2, 2] +
                                jak[0, 1] * jak[1, 2] * jak[2, 0] +
                                jak[0, 2] * jak[1, 0] * jak[2, 1]
                            ) -
                            (
                                jak[0, 2] * jak[1, 1] * jak[2, 0] +
                                jak[0, 1] * jak[1, 0] * jak[2, 2] +
                                jak[0, 0] * jak[1, 2] * jak[2, 1]
                            );
                }


                double[] col = new double[3];

                for (int i = 0; i < 27; i++)
                {
                    for (int phi = 0; phi < 20; phi++)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            col[k] = DFIABG[i, k, phi];
                        }
                        double[,] matrix = new double[3, 3] {
                            { dxyzabg[0,0,i], dxyzabg[1,0,i], dxyzabg[2,0,i] },
                            { dxyzabg[0,1,i], dxyzabg[1,1,i], dxyzabg[2,1,i] },
                            { dxyzabg[0,2,i], dxyzabg[1,2,i], dxyzabg[2,2,i] }
                        };
                        double[] gaussianSolve = Gauss.Solve(matrix, col);

                        for (int k = 0; k < 3; k++)
                        {
                            dfixyz[i, phi, k] = gaussianSolve[k];
                        }
                    }
                }

                double[,][,] mge = new double[3, 3][,];

                mge[0, 0] = one_one(dfixyz, dj,  lam,  v,  mu,  c);
                mge[1, 1] = two_two(dfixyz, dj, lam, v, mu, c);
                mge[2, 2] = three_three(dfixyz, dj, lam, v, mu, c);

                mge[0, 1] = one_two(dfixyz, dj, lam, v, mu, c);
                mge[0, 2] = one_three(dfixyz, dj, lam, v, mu, c);
                mge[1, 2] = two_three(dfixyz, dj, lam, v, mu, c);

                mge[1, 0] = rotate(mge[0, 1]);
                mge[2, 0] = rotate(mge[0, 2]);
                mge[2, 1] = rotate(mge[1, 2]);

                int x, y, localX, localY, globalX, globalY;
                
                for (int i = 0; i < 60; i++)
                {
                    for (int j = 0; j < 60; j++)
                    {
                        x = i / 20;
                        y = j / 20;

                        localX = i % 20;
                        localY = j % 20;

                        globalX = (NT[number][localX]) * 3 + x;
                        globalY = (NT[number][localY]) * 3 + y;

                        MG[globalX, globalY] += mge[x, y][localX, localY];
                    }
                }
            }
            return MG;
        }

        private static double[,] one_one(double[,,] dfixyz, double[] dj,double lam,double v,double mu,double[] c)
        {
            double[,] res = new double[20, 20];
            for (int i = 0; i < 20; i++)
            {
                for (int j = 0; j < 20; j++)
                {
                    if (i > j)
                    {
                        res[i, j] = res[j, i];
                    }
                    else
                    {
                        double sum = 0;
                        int counter = 0;
                        for (int k = 0; k < 3; k++)
                        {
                            for (int l = 0; l < 3; l++)
                            {
                                for (int m = 0; m < 3; m++)
                                {
                                    sum += (
                                            (lam * (1 - v) * (dfixyz[counter, i, 0] * dfixyz[counter, j, 0]))
                                            +
                                            (mu * (dfixyz[counter, i, 1] * dfixyz[counter, j, 1] + dfixyz[counter, i, 2] * dfixyz[counter, j, 2]))
                                        ) * Math.Abs(dj[counter]) * c[m] * c[l] * c[k];
                                    ++counter;
                                }
                            }
                        }
                        res[i, j] = sum;
                    }
                }
            }
            return res;
        }
        private static double[,] two_two(double[,,] dfixyz, double[] dj, double lam, double v, double mu, double[] c)
        {
            double[,] res = new double[20, 20];
            for (int i = 0; i < 20; i++)
            {
                for (int j = 0; j < 20; j++)
                {
                    if (i > j)
                    {
                        res[i, j] = res[j, i];
                    }
                    else
                    {
                        double sum = 0;
                        int counter = 0;
                        for (int k = 0; k < 3; k++)
                        {
                            for (int l = 0; l < 3; l++)
                            {
                                for (int m = 0; m < 3; m++)
                                {
                                    sum += (
                                            (lam * (1 - v) * (dfixyz[counter, i, 1] * dfixyz[counter, j, 1]))
                                            +
                                            (mu * (dfixyz[counter, i, 0] * dfixyz[counter, j, 0] + dfixyz[counter, i, 2] * dfixyz[counter, j, 2]))
                                        ) * Math.Abs(dj[counter]) * c[m] * c[l] * c[k];
                                    ++counter;

                                }
                            }
                        }
                        res[i, j] = sum;
                    }
                }
            }
            return res;
        }
        private static double[,] three_three(double[,,] dfixyz, double[] dj, double lam, double v, double mu, double[] c)
        {
            double[,] res = new double[20, 20];
            for (int i = 0; i < 20; i++)
            {
                for (int j = 0; j < 20; j++)
                {
                    if (i > j)
                    {
                        res[i, j] = res[j, i];
                    }
                    else
                    {
                        double sum = 0;
                        int counter = 0;
                        for (int k = 0; k < 3; k++)
                        {
                            for (int l = 0; l < 3; l++)
                            {
                                for (int m = 0; m < 3; m++)
                                {
                                    sum += (
                                            (lam * (1 - v) * (dfixyz[counter, i, 2] * dfixyz[counter, j, 2]))
                                            +
                                            (mu * (dfixyz[counter, i, 0] * dfixyz[counter, j, 0] + dfixyz[counter, i, 1] * dfixyz[counter, j, 1]))
                                        ) * Math.Abs(dj[counter]) * c[m] * c[l] * c[k];
                                    ++counter;
                                }
                            }
                        }
                        res[i, j] = sum;
                    }
                }
            }
            return res;
        }

        private static double[,] one_two(double[,,] dfixyz, double[] dj, double lam, double v, double mu, double[] c)
        {
            double[,] res = new double[20, 20];
            for (int i = 0; i < 20; i++)
            {
                for (int j = 0; j < 20; j++)
                {
                    double sum = 0;
                    int counter = 0;
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            for (int m = 0; m < 3; m++)
                            {
                                sum += (
                                    (lam * v * (dfixyz[counter, i, 0] * dfixyz[counter, j, 1]))
                                      +
                                    (mu * (dfixyz[counter, i, 1] * dfixyz[counter, j, 0]))
                                    ) * Math.Abs(dj[counter]) * c[m] * c[l] * c[k];
                                ++counter;
                            }
                        }
                    }
                    res[i, j] = sum;

                }
            }
            return res;
        }
        private static double[,] one_three(double[,,] dfixyz, double[] dj, double lam, double v, double mu, double[] c)
        {
            double[,] res = new double[20, 20];
            for (int i = 0; i < 20; i++)
            {
                for (int j = 0; j < 20; j++)
                {
                    double sum = 0;
                    int counter = 0;
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            for (int m = 0; m < 3; m++)
                            {
                                sum += (
                                    (lam * v * (dfixyz[counter, i, 0] * dfixyz[counter, j, 2]))
                                      +
                                    (mu * (dfixyz[counter, i, 2] * dfixyz[counter, j, 0]))
                                    ) * Math.Abs(dj[counter]) * c[m] * c[l] * c[k];
                                ++counter;
                            }
                        }
                    }
                    res[i, j] = sum;

                }
            }
            return res;
        }
        private static double[,] two_three(double[,,] dfixyz, double[] dj, double lam, double v, double mu, double[] c)
        {
            double[,] res = new double[20, 20];
            for (int i = 0; i < 20; i++)
            {
                for (int j = 0; j < 20; j++)
                {
                    double sum = 0;
                    int counter = 0;
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            for (int m = 0; m < 3; m++)
                            {
                                sum += (
                                    (lam * v * (dfixyz[counter, i, 1] * dfixyz[counter, j, 2]))
                                      +
                                    (mu * (dfixyz[counter, i, 2] * dfixyz[counter, j, 1]))
                                    ) * Math.Abs(dj[counter]) * c[m] * c[l] * c[k];
                                ++counter;
                            }
                        }
                    }
                    res[i, j] = sum;

                }
            }
            return res;
        }

        private static double[,] rotate(double[,] toRotate)
        {
            double[,] res = new double[20, 20];
            for (int i = 0; i < 20; i++)
            {
                for (int j = 0; j < 20; j++)
                {
                    res[i, j] = toRotate[j, i];
                }
            }
            return res;
        }

        public static double[,] improveMG(int[] ZU,double[,] MG)
        {
            int index;
            for (int i = 0; i < ZU.Length; i++)
            {
                index = ZU[i] * 3;
                for (int j = 0; j < 3; j++)
                {
                    MG[index + j, index + j] = Globals.BIG_NUMBER;
                }
            }
            return MG;
        }
    }
}
