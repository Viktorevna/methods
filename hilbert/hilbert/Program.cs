using System;

namespace hilbert
{
    class Program
    {
        public static double Function(double x)
        {
            return 3 * x - Math.Cos(x + 1);
        }
        public static double Factorial(int n)
        {
            int result = 1;
            for (int i = 1; i < n+1; i++)
            {
                result *= i;
            }
            return result;
        }
        public static double Integrate(Func<double, double> f, double a, double b, int accuracy)
        {
            double result = 0;
            double dx = (b - a) / accuracy;
            double x = a + dx;
            for (int i=0; i <= accuracy; i++)
            {
                result += f(x - dx / 2) * dx;
                x += dx;
                //Console.WriteLine(x);
            }
            //Console.WriteLine(result);
            return result;
        }
        public static double Derivative(Func<double, double> f, double x, int number)
        {
            double dx = 0.0001;
            double[,] DerArray = new double[number + 1, number + 1];
            //if (number == 0)
                //return f(x);
            for (int i = 0; i < DerArray.GetLength(0); i++)
                DerArray[i, 0] = f(x + i * dx);
            for (int j = 1; j < DerArray.GetLength(1); j++)
                for (int i = 0; i < DerArray.GetLength(0) - j; i++)
                    DerArray[i, j] = (DerArray[i, j - 1] - DerArray[i + 1, j - 1]) / dx;
            return DerArray[0, number];
        }
        public static double LezhandrPolynomial(double x, int number)
        { 
            return 1 / (Factorial(number) * Math.Pow(2, number)) * 
                Derivative((double a) => Math.Pow((1 - a * a), number), x, number);
        }
        public static double Coefficient(double a, double b, double x, int number, int accuracy)
        {
            return Integrate((double z) => Function(z) * LezhandrPolynomial(z, number), a, b, accuracy) /
                Integrate((double z) => LezhandrPolynomial(z, number) * LezhandrPolynomial(z, number), a, b, accuracy);
        }
        public static double HPolynomial(double x, int n, double a, double b, int accuracy)
        {
            double result = 0;
            for (int i = 0; i <= n; i++)
            {
                result += Coefficient(a, b, x, i, accuracy) * LezhandrPolynomial(x, i);
            }
            return result;
        }
        static void Main(string[] args)
        {
            int polynomdegree = 1;
            int pointsnum = 40;
            int accuracy = 10000;
            double a = -1;
            double b = 1;
            output(a, b, polynomdegree, pointsnum, accuracy);
            Console.ReadLine();
        }
        public static void output(double a, double b, int polynomdegree, int pointsnum, int accuracy)
        {
            Console.WriteLine("Approximation in Hilbert spaces");
            Console.WriteLine("--------------------------------------------------------------");
            Console.WriteLine("x\t\tf(x)\t\tH(x)\t\t|f(x)-H(x)|");
            double x = a;
            int i = 0;
            do
            {
                i++;
                double pol = HPolynomial(x, polynomdegree, a, b, accuracy);
                var f = Function(x);
                Console.WriteLine($"{x:0.###}\t\t{f:0.###}\t\t{pol:0.###}\t\t{(Math.Abs(f - pol)):0.###}");
                x += (b - a) / pointsnum;
            } while (i <= pointsnum);
        }
    }
}
