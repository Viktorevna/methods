using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace least_square_method
{
    class Program
    {
        public static Random random = new Random();
        public static double function(double x)
        {
            return 3 * x - Math.Cos(x + 1);
        }
        public static double inaccurate_function(double x)
        {
            return function(x) + Math.Pow(-1, random.Next())*0.01*random.NextDouble()*10;
        }
        public static double[,] create_table(Func<double,double> func, double a, double b, int pointsnum)
        {
            double[,] table = new double[pointsnum, 2];
            double x = a;
            for (int i = 0; i < pointsnum; i++)
            {
                table[i, 0] = x;
                table[i, 1] = func(x);
                x += (b - a) / (pointsnum-1);
            }
            return table;
        }
        public static Matrix table_to_X(double[,] table)
        {
            double[] Main = new double[table.GetLength(0)];
            for (int i = 0; i < Main.GetLength(0); i++)
                Main[i] = table[i, 0];
            return new Matrix(Main);
        }
        public static Matrix table_to_Y(double[,] table)
        {
            double[] Main = new double[table.GetLength(0)];
            for (int i = 0; i < Main.GetLength(0); i++)
                Main[i] = table[i, 1];
            return new Matrix(Main);
        }
        public static double fifunction(double x, int number)
        {
            return Math.Pow(x, number);
        }
        public static Matrix get_matrix_Q(Func<double,int,double> fi, int polynomdegree, Matrix X)
        {
            double[,] Main = new double[X.Height, polynomdegree+1];
            for (int i = 0; i < Main.GetLength(0); i++)
                for (int j = 0; j < Main.GetLength(1); j++)
                    Main[i, j] = fi(X.Main[0, i],j);
            return new Matrix(Main);
        }
        public static double MNKPolynomial(double x, Matrix A)
        {
            //Matrix H = Q.Transposed() * Q;
            //Matrix B = Q.Transposed() * Y;
            //Matrix A = Matrix.GausMethod(H, B);
            double[] fiMain = new double[A.Height];
            for (int i = 0; i < fiMain.GetLength(0); i++)
                fiMain[i] = fifunction(x, i);
            return Matrix.EuclideanNorm(A, new Matrix(fiMain));
        }
        static void Main(string[] args)
        {
            int polynomdegree = 10;
            int pointsnum = 100;
            double a = -1;
            double b = 1;
            output(a, b, polynomdegree, pointsnum, 40);
        }
        public static void output(double a, double b, int polynomdegree, int pointsnum, int accuracy)
        {
            double[,] table = create_table(inaccurate_function, a, b, pointsnum);
            Matrix X = table_to_X(table);
            Matrix Y = table_to_Y(table);
            Matrix Q = get_matrix_Q(fifunction, polynomdegree, X.Transposed());
            Console.WriteLine("LEAST SQUARE METHOD");
            Console.WriteLine("###################");
            Console.WriteLine("x\t\tf(x)\t\tMNK(x)\t\t|f(x)-MNK(x)|");
            double x = a;
            int i = 0;
            Matrix H = Q.Transposed() * Q;
            Matrix B = Q.Transposed() * Y;
            Matrix A = Matrix.GausMethod(H, B);
            do
            {
                i++;
                double pol = MNKPolynomial(x,A);
                var func = function(x);
                Console.WriteLine($"{x:0.###}\t\t{func:0.###}\t\t{pol:0.###}\t\t{(Math.Abs(func - pol)):0.###}");
                x += (b - a) / accuracy;
            } while (i <= accuracy);
        }
    }
}
