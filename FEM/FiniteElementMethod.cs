using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using Newtonsoft.Json;

namespace FEM
{
    class FiniteElementMethod
    {
        public double[,,][] matrix;

        public int x;
        public int y;
        public int z;

        public int m;
        public int n;
        public int k;

        public double E = 1.0;
        public double v;
        public double lam;
        public double mu;



        public int nqp;

        public double[][] AKT;

        public int[,,] local_to_global;
        public Dictionary<int, int[]> adapter = Globals.magicDictionary;
        public int[][] PAdapter = new int[6][] {
            new int[8] { 0, 1, 5, 4, 8, 13, 16, 12},
            new int[8] { 1, 2, 6, 5, 9, 14, 17, 13},
            new int[8] { 2, 3, 7, 6, 10, 15, 18, 14 },
            new int[8] { 3, 0, 4, 7, 11, 12, 19, 15},
            new int[8] { 0, 1, 2, 3, 8, 9, 10, 11},
            new int[8] { 4, 5, 6, 7, 16, 17, 18, 19}
        };

        public double[][] GaussNodes = Globals.GaussNodes;
        public double[,,] DFIABG = Globals.DFIABG;
        public double[,,] DFIABG_P = Globals.DFIABG_P;

        public double[] c = new double[3] { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

        int nel;
        public int[][] NT;

        public double[,] MG;
        public double[] F;

        public int[] ZU;
        public double[,] ZP;

        public double[,,] DPSITE = Globals.DPSITE;
        public double[,] PSIET = new double[8, 9];

        private double SCALE_X;
        private double SCALE_Y;
        private double SCALE_Z;

        private double[] U;

        public FiniteElementMethod(int _x, int _y, int _z, int _m, int _n, int _k, double _v)
        {
            x = _x;
            y = _y;
            z = _z;

            m = _m;
            n = _n;
            k = _k;

            v = _v;
            lam = E / ((1 + v) * (1 - 2 * v));
            mu = E / (2 * (1 + v));

            SCALE_X = (double)x / m;
            SCALE_Y = (double)y / n;
            SCALE_Z = (double)z / k;

            matrix = new double[m * 2 + 1, n * 2 + 1, k * 2 + 1][];
            nqp = 0;
           
            local_to_global = new int[m * 2 + 1, n * 2 + 1, k * 2 + 1];

            nel = m * n * k;
            NT = new int[nel][];
        }

        public Tuple<double[][], double[][]> Start()
        {
            KeyValuePair<double[,,][], int> mainVert = VertexFill.fillMatrixWithMainVertexes(matrix, nqp, k, n, m, SCALE_X, SCALE_Y, SCALE_Z);
            matrix = mainVert.Key;
            nqp = mainVert.Value;
            KeyValuePair<double[,,][], int> interVert = VertexFill.fillMatrixWithIntermidiateVertexes(matrix, nqp, k, n, m, SCALE_X, SCALE_Y, SCALE_Z);
            matrix = interVert.Key;
            nqp = interVert.Value;

            MG = new double[3 * nqp, 3 * nqp];
            F = new double[3 * nqp];

            AKT = CreateAKT.createAKT(nqp, k, n, m, matrix, local_to_global);
            ZU = CreateAKT.createZU(AKT);
            ZP = CreateAKT.createZP(m, n, nel);
            createNT();
            MG = createMG.getMG(nel, NT, AKT, DFIABG, lam, v, mu, c, nqp);
            MG = createMG.improveMG(ZU, MG);
            PSIET = Presure.createPSI();
            F = Presure.createF(nqp, nel, m, n, AKT, NT, PSIET, c);

            U = Gauss.Solve(MG, F);

            double[][] AKTres = new double[nqp][];
            for (int i = 0; i < nqp; i++)
            {
                double[] prev = AKT[i];
                double[] point = U.Skip(i * 3).Take(3).ToArray();
                AKTres[i] = new double[3] { Math.Round(prev[0] + point[0], 4), Math.Round(prev[1] + point[1], 4), Math.Round(prev[2] + point[2], 4) };
            }

            using (StreamWriter sw = new StreamWriter("C:\\FEM\\frontend\\FEMpoints.txt", false, System.Text.Encoding.Default))
            {
                sw.WriteLine(JsonConvert.SerializeObject((from a in AKTres select new { x = a[0], y = a[1], z = a[2], })));
            }

            using (StreamWriter sw = new StreamWriter("C:\\FEM\\frontend\\start.txt", false, System.Text.Encoding.Default))
            {
                sw.WriteLine(JsonConvert.SerializeObject((from a in AKT select new { x = a[0], y = a[1], z = a[2], })));
            }

            return new Tuple<double[][], double[][]>(AKT, AKTres);
        }

        private void createNT()
        {
            for (int figure = 0; figure < nel; figure++)
            {
                int current_element = figure;
                int Z_COOORDINATE = (int)(figure / (m * n));
                current_element %= (m * n);
                int Y_COORDINATE = (int)(current_element / m);
                current_element %= (m);
                int X_COORDINATE = current_element;

                createIndependentElements(X_COORDINATE * 2, Y_COORDINATE * 2, Z_COOORDINATE * 2, figure);
            }
        }

        private void createIndependentElements(int x, int y, int z, int belongElementNumber)
        {
            int[] globalCoordinates = new int[20];
            int[] delta;
            for (int i = 0; i < 20; i++)
            {
                delta = adapter[i];
                globalCoordinates[i] = local_to_global[x + delta[0], y + delta[1], z + delta[2]];
            }
            NT[belongElementNumber] = globalCoordinates;
        }
    }
}
