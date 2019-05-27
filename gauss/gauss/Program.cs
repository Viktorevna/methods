using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NewtonCotesMethod
{
    static class Mathematics
    {
        public static double[] GausMethod(double[,] matrix, double[] b)
        {
            int width = matrix.GetLength(1);
            int height = matrix.GetLength(0);
            for (int j = 0; j < height; j++)
            {
                double cur_diag = matrix[j, j];
                if (cur_diag != 0)
                {
                    for (int i = 0; i < width; i++)
                        matrix[j, i] /= cur_diag;
                    b[j] /= cur_diag;
                }
                for (int k = 0; k < height; k++)
                {
                    if (k != j)
                    {
                        double cur = matrix[k, j];
                        for (int i = 0; i < width; i++)
                            matrix[k, i] -= matrix[j, i] * cur;
                        b[k] -= b[j] * cur;
                    }
                }
            }
            for (int j = height - 1; j > -1; j--)
            {
                double cur_diag = matrix[j, j];
                if (cur_diag != 0)
                {
                    for (int i = 0; i < width; i++)
                        matrix[j, i] /= cur_diag;
                    b[j] /= cur_diag;
                }
                for (int k = height - 1; k > -1; k--)
                {
                    if (k != j)
                    {
                        double cur = matrix[k, j];
                        for (int i = 0; i < width; i++)
                            matrix[k, i] -= matrix[j, i] * cur;
                        b[k] -= b[j] * cur;
                    }
                }
            }
            return b;
        }

        public static double Factorial(double x)
        {
            if (x == 0 || x == 1)
                return 1;
            return x * Factorial(x - 1);
        }
        public static double C(int n, int k)
        {
            return Factorial(n) / (Factorial(k) * Factorial(n - k));
        }

        public static double[] GetMuArray(int n, double a, double b, double alpha, double beta)
        {
            double[] Mu = new double[n];
            double c;
            if (alpha == 0)
                c = b;
            else
                c = a;
            for (int k = 0; k < n; k++)
            {
                double sum = 0;
                for (int i = 1; i <= k; i++)
                    sum += C(k, i) * Math.Pow(-1, i + 1) * Math.Pow(c, i) * Mu[k - i];
                Mu[k] = 1 / (k + 1 - alpha - beta) * Math.Pow(b - a, k + 1 - alpha - beta) + sum;
            }

            return Mu;
        }
    }

    static class MidpointMethod
    {
        public static double Integrate(Func<double, double> func,
            Func<double, double, double, double, double, double> P,
            double a, double b, double alpha, double beta, int n)
        {
            double sum = 0;
            double h = (b - a) / (n - 1);
            double x = a;
            for (int i = 0; i < n - 1; i++)
            {
                sum += func(x + 0.5 * h) * P(x + 0.5 * h, a, b, alpha, beta) * h;
                x += h;
            }
            return sum;
        }
    }

    static class NewtonCotesMethod
    {
        public static double Integrate(Func<double, double> func,
            double a, double b, double alpha, double beta, int n)
        {
            var Mu = Mathematics.GetMuArray(n, a, b, alpha, beta);
            double[] X = new double[n];
            double h = (b - a) / (n - 1);
            for (int i = 0; i < n; i++)
                X[i] = a + h * i;
            double[,] matrix = new double[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    matrix[i, j] = Math.Pow(X[j], i);
            var A = Mathematics.GausMethod(matrix, Mu);

            double integralSum = 0;
            for (int i = 0; i < n; i++)
            {
                integralSum += func(X[i]) * A[i];
            }
            return integralSum;
        }
    }

    static class GaussMethod
    {
        public static double CalculatePolynomial(double[] coefs, double x)
        {
            double result = 0;
            for (int i = 0; i < coefs.Length; i++)
                result += coefs[coefs.Length - 1 - i] * Math.Pow(x, i);
            return result;
        }

        public static double[] GetSquarePolynomialRoots(double[] coefs)
        {
            double D = coefs[1] * coefs[1] - 4 * coefs[0] * coefs[2];
            if (D < 0) D = 0;
            var roots = new double[2];
            roots[0] = (-coefs[1] + Math.Sqrt(D)) / (2 * coefs[0]);
            roots[1] = (-coefs[1] - Math.Sqrt(D)) / (2 * coefs[0]);
            return roots;
        }
        public static double Integrate(Func<double, double> func,
            double a, double b, double alpha, double beta, int n)
        {
            var Mu = Mathematics.GetMuArray(2 * n, a, b, alpha, beta);
            var SLAE = new double[n, n];
            var bSLAE = new double[n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    SLAE[i, j] = Mu[j + i];
                }
                bSLAE[i] = -Mu[n + i];
            }
            var acoefs = Mathematics.GausMethod(SLAE, bSLAE);
            var polcoefs = new double[acoefs.Length + 1];
            for (int i = 0; i < polcoefs.Length - 1; i++)
                polcoefs[i] = acoefs[i];
            polcoefs[polcoefs.Length - 1] = 1;

            var nodes = GetSquarePolynomialRoots(polcoefs);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    SLAE[i, j] = Math.Pow(nodes[j], i);
                }
                bSLAE[i] = Mu[i];
            }
            var intcoefs = Mathematics.GausMethod(SLAE, bSLAE);
            double result = 0;
            for (int i = 0; i < n; i++)
                result += intcoefs[i] * func(nodes[i]);
            return result;
        }
    }
    class Program
    {
        public static double Function(double x)
        {
            //return 1.3 * Math.Cos(3.5 * x) * Math.Exp(2 * x / 3) + 6 * Math.Sin(4.5 * x) * Math.Exp(-x / 8) + 5 * x;
            return 0.5 * Math.Cos(2 * x) * Math.Exp(0.4 * x) + 2.4 * Math.Sin(1.5 * x) * Math.Exp(-6 * x) + 6 * x;

        }

        public static double P(double x, double a, double b, double alpha, double beta)
        {
            //return 21;
            return 1 / (Math.Pow(x - a, alpha) * Math.Pow(b - x, beta));
        }

        static void Main(string[] args)
        {
            Stopwatch SW = new Stopwatch();
            double a = 1.1;
            double b = 2.5;
            double alpha = 0.4;
            double beta = 0;
            SW.Start();
            SW.Stop();
            SW.Restart();
            var n = 2;
            var integralSum = GaussMethod.Integrate(Function, a, b, alpha, beta, n);
            SW.Stop();
            output("Метод Гаусса", integralSum, n, SW.ElapsedTicks);
            SW.Restart();
            n = 3;
            integralSum = NewtonCotesMethod.Integrate(Function, a, b, alpha, beta, n);
            SW.Stop();
            output("Метод Ньютона-Котса", integralSum, n, SW.ElapsedTicks);
            for (int i = 10; i <= 10000; i *= 10)
            {
                SW.Restart();
                integralSum = MidpointMethod.Integrate(Function, P, a, b, alpha, beta, i);
                SW.Stop();
                output("Метод средней точки", integralSum, i, SW.ElapsedTicks);
            }
            Console.ReadLine();
        }

        static void output(string title, double result, int dimension, long time)
        {
            Console.WriteLine(title + ": " + "n = " + dimension + ", результат = " + result + ", время = " + time);
        }
    }
}
