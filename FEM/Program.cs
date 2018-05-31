using System;

namespace FEM
{
    class Program
    {
        static void Main(string[] args)
        {
            FiniteElementMethod solve = new FiniteElementMethod(50, 50, 50, 3, 3, 3, 0.3);
            solve.Start();
        }
    }
}
