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
        public static double[] GausMethod(double[,] matrix, double[] B)
        {
            int width = matrix.GetLength(1);
            int height = matrix.GetLength(0);
            double[] b = new double[B.Length];
            B.CopyTo(b,0);
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
        public static double[] GetMuArray(int n, double a, double b, double start, double end, double alpha, double beta)
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
                Mu[k] = 1 / (k + 1 - alpha - beta) * (Math.Pow(end - a, k + 1 - alpha - beta) - Math.Pow(start - a, k + 1 - alpha - beta)) + sum;
            }

            return Mu;
        }
        public static int GetGoodNum(double start, double end, double step)
        {
            return Convert.ToInt32(Math.Ceiling((end - start) / step))+1;
        }
        public static double AccuracyRunge(double intSum1, double intSum2, int ada)
        {
            return (intSum2 - intSum1) / (1 - Math.Pow(2, -ada));
        }
        public static double GetOptStep(double a, double b, int segmNum, double accuracy, int ada, double epsilon)
        {
            return 0.95 * ((b - a) / (segmNum)) * Math.Pow(epsilon / Math.Abs(accuracy), 1/Convert.ToDouble(ada));
        }
        public static double GetOptStepRichardson(Func<Func<double, double>,
            double, double, double, double, int, int, double> method, Func<double, double> func,
            double a, double b, double alpha, double beta, int n, int r, int ada, double epsilon) 
        {
            double startStep = (b - a);
            double[,] A = new double[r + 1, r + 1];
            double[] B = new double[r + 1];
            double[] X = new double[r + 1];
            do
            {
                for (int i = 0; i < A.GetLength(0); i++)
                {
                    double step = startStep / Math.Pow(2, i);
                    for (int j = 0; j < A.GetLength(0); j++)
                    {
                        if (j == 0) A[i, j] = -1;
                        else
                            A[i, j] = Math.Pow(step, ada + j - 1);
                    }
                    B[i] = -method(func, a, b, alpha, beta, n, Mathematics.GetGoodNum(a, b, step));
                }
                X = Mathematics.GausMethod(A, B);
                startStep /= 2;
            } while (Math.Abs(X[0]+B[r])>epsilon);

            return startStep/Math.Pow(2,r-1);
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
                sum += func(x + 0.5d * h) * P(x + 0.5d * h, a, b, alpha, beta) * h;
                x += h;
            }
            return sum;
        }
    }

    static class NewtonCotesMethod
    {
        public static double Integrate(Func<double, double> func,
            double a, double b, double start, double end, double alpha, double beta, int n)
        {
            var Mu = Mathematics.GetMuArray(n, a, b, start, end, alpha, beta);
            double[] X = new double[n];
            double h = (end - start) / (n - 1);
            for (int i = 0; i < n; i++)
                X[i] = start + h * i;
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
        public static double CompoundIntegrate(Func<double, double> func,
            double a, double b, double alpha, double beta, int n, int segmNum)
        {
            double result = 0;
            double h = (b - a) / (segmNum - 1);
            for (int i = 0; i < segmNum - 1; i++)
            {
                var intSum = Integrate(func, a, b, a + i * h, a + (i + 1) * h, alpha, beta, n);
                result += intSum;
            }
            return result;
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
            roots[0] = (-coefs[1] + Math.Sqrt(D)) / (2 * coefs[2]);
            roots[1] = (-coefs[1] - Math.Sqrt(D)) / (2 * coefs[2]);
            return roots;
        }
        public static double Integrate(Func<double, double> func,
            double a, double b, double start, double end, double alpha, double beta, int n)
        {
            var Mu = Mathematics.GetMuArray(2 * n, a, b, start, end, alpha, beta);
            var SLAE = new double[n, n];
            var bSLAE = new double[n];
            for (int s = 0; s < n; s++)
            {
                for (int j = 0; j < n; j++)
                {
                    SLAE[s, j] = Mu[j + s];
                }
                bSLAE[s] = -Mu[n + s];
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
        public static double CompoundIntegrate(Func<double, double> func,
            double a, double b, double alpha, double beta, int n, int segmNum)
        {
            double result = 0;
            double h = (b - a) / (segmNum - 1);
            for (int i = 0; i < segmNum - 1; i++)
                result += Integrate(func, a, b, a + i * h, a + (i + 1) * h, alpha, beta, n);
            return result;
        }

    }

    class Program
    {
        public static double Function(double x)
        {
            return 0.5d * Math.Cos(2d * x) * Math.Exp(0.4d * x) + 2.4d * Math.Sin(1.5d * x) * Math.Exp(-6d * x) + 6d * x;

        }

        public static double P(double x, double a, double b, double alpha, double beta)
        {
            return 1d / (Math.Pow(x - a, alpha) * Math.Pow(b - x, beta));
        }

        static void Main(string[] args)
        {
            Stopwatch SW = new Stopwatch();

            double accRunge;
            double richOptStep;

            int segmNum = 100;
            double integralSum;
            int n;
            double a = 1.1d;
            double b = 2.5d;
            double alpha = 0.4d;
            double beta = 0;
            SW.Start();
            SW.Stop();
            //SW.Restart();
            //n = 2;
            //integralSum = GaussMethod.CompoundIntegrate(Function, a, b, alpha, beta, n, segmNum);
            //SW.Stop();
            //output("CQF (Gauss)", integralSum, n, SW.ElapsedTicks);

            //richOptStep = Mathematics.GetOptStepRichardson(GaussMethod.CompoundIntegrate, Function, a, b, alpha, beta, n, 2, n - 1, 1e-6);

            //Console.WriteLine("Optimal step (Richardson): " + richOptStep);
            //Console.WriteLine(" => Optimal n: " + Mathematics.GetGoodNum(a, b, richOptStep));
            SW.Restart();
            n = 3;
            integralSum = NewtonCotesMethod.CompoundIntegrate(Function, a, b, alpha, beta, n, segmNum);
            SW.Stop();
            output("CQF (Newton-Cotes)", integralSum, n, SW.ElapsedTicks);

            accRunge = Mathematics.AccuracyRunge(integralSum, NewtonCotesMethod.CompoundIntegrate(Function, a, b, alpha, beta, n, segmNum*2), n - 1);

            var optStep = Mathematics.GetOptStep(a, b, segmNum, accRunge, n - 1, 1e-6);
            Console.WriteLine("Optimal step (Runge): " + optStep);
            Console.WriteLine(" => Optimal n: " + Mathematics.GetGoodNum(a, b, optStep));

            richOptStep = Mathematics.GetOptStepRichardson(NewtonCotesMethod.CompoundIntegrate, Function, a, b, alpha, beta, n, 2, n - 1, 1e-6);

            Console.WriteLine("Optimal step (Richardson): " + richOptStep);
            Console.WriteLine(" => Optimal n: " + Mathematics.GetGoodNum(a, b, richOptStep));
            Console.WriteLine("Error (Runge): " + accRunge);
        }

        static void output(string title, double result, int dimension, long time)
        {
            Console.WriteLine("\n"+title + ": " + "n = " + dimension + ", result = " + result + ", time = " + time);
        }
    }
}
