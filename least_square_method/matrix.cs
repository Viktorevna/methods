using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace least_square_method
{
    class Matrix
    {
        public double[,] Main;
        public int Width { get; }
        public int Height { get; }
        public Matrix(int Width, int Height)
        {
            Main = new double[Width, Height];
            this.Width = Width;
            this.Height = Height;
        }
        public Matrix()
        {
            Main = null;
            this.Width = 0;
            this.Height = 0;
        }
        public Matrix(double[,] A)
        {
            Main = A;
            this.Width = A.GetLength(0);
            this.Height = A.GetLength(1);
        }
        public Matrix(double[] A)
        {
            this.Width = A.GetLength(0);
            this.Height = 1;
            Main = new double[Width, Height];
            for (int i = 0; i < Width; i++)
                Main[i, 0] = A[i];
        }
        public void Write()
        {
            for (int i = 0; i < Width; i++)
            {
                if (i != 0) Console.WriteLine();
                for (int j = 0; j < Height; j++)
                    Console.Write(String.Format("{0:0.000000}", Main[i, j]) + "\t");
            }
            Console.WriteLine();

        }
        public static double Norm1(Matrix A)
        {
            double Result = 0;
            for (int j = 0; j < A.Height; j++)
            {
                double Sum = 0;
                for (int i = 0; i < A.Width; i++)
                    Sum += Math.Abs(A.Main[i, j]);
                if (Result < Sum) Result = Sum;
            }
            return Result;
        }
        public static double EuclideanNorm(Matrix A, Matrix B)
        {
            Matrix R = A * B;
            return OneOnematrixToDouble(R);
        }

        public static double NormV(Matrix A)
        {
            double Result = 0;
            for (int j = 0; j < A.Height; j++)
            {
                double Sum = 0;
                for (int i = 0; i < A.Width; i++)
                    Sum += A.Main[i, j] * A.Main[i, j];
                Result = Math.Sqrt(Sum);
            }
            return Result;
        }
        public static double OneOnematrixToDouble(Matrix A)
        {
            if (A.Width == 1 && A.Height == 1) return A.Main[0, 0];
            else
            {
                Console.WriteLine("Need 1X1 dimenstion!");
                return 0;
            }
        }
        public Matrix Transposed()
        {
            Matrix Result = new Matrix(Height, Width);
            for (int i = 0; i < Width; i++)
                for (int j = 0; j < Height; j++)
                    Result.Main[j, i] = Main[i, j];
            return Result;
        }
        public static Matrix operator +(Matrix A, Matrix B)
        {
            if ((A.Width == B.Width) && (A.Height == B.Height))
            {
                Matrix Result = new Matrix(A.Width, A.Height);
                for (int i = 0; i < A.Width; i++)
                    for (int j = 0; j < A.Height; j++)
                        Result.Main[i, j] = A.Main[i, j] + B.Main[i, j];
                return Result;
            }
            else
            {
                Console.WriteLine("Not correct dimensions!");
                return new Matrix();
            }
        }
        public static Matrix operator *(double b, Matrix A)
        {

            Matrix Result = new Matrix(A.Width, A.Height);
            for (int i = 0; i < A.Width; i++)
                for (int j = 0; j < A.Height; j++)
                    Result.Main[i, j] = A.Main[i, j] * b;
            return Result;
        }

        public static Matrix operator *(Matrix A, Matrix B)
        {
            if (A.Height == B.Width)
            {
                Matrix Result = new Matrix(A.Width, B.Height);
                for (int i = 0; i < A.Width; i++)
                    for (int j = 0; j < B.Height; j++)
                    {
                        double Sum = 0;
                        for (int r = 0; r < A.Height; r++)
                            Sum += A.Main[i, r] * B.Main[r, j];
                        Result.Main[i, j] = Sum;
                    }
                return Result;
            }
            else
            {
                throw new Exception("Not correct dimensions! [*] " + A.Height + " !=" + B.Width);
            }
        }
        public static double[,] CloneArray(double[,] A)
        {
            double[,] Result = new double[A.GetLength(0), A.GetLength(1)];
            for (int i = 0; i < A.GetLength(0); i++)
                for (int j = 0; j < A.GetLength(1); j++)
                    Result[i, j] = A[i, j];
            return Result;
        }
        public static Matrix GausMethod(Matrix A, Matrix B)
        {
            var matrix = CloneArray(A.Main);
            var b = CloneArray(B.Main);
            int width = matrix.GetLength(1);
            int height = matrix.GetLength(0);
            for (int j = 0; j < height; j++)
            {
                double cur_diag = matrix[j, j];
                if (cur_diag != 0)
                {
                    for (int i = 0; i < width; i++)
                        matrix[j, i] /= cur_diag;
                    b[j,0] /= cur_diag;
                }

                //MatrixPrint(matrix);
                //MatrixPrint(b);
                for (int k = 0; k < height; k++)
                {
                    if (k != j)
                    {
                        double cur = matrix[k, j];
                        for (int i = 0; i < width; i++)
                            matrix[k, i] -= matrix[j, i] * cur;
                        //MatrixPrint(matrix);
                        b[k,0] -= b[j,0] * cur;
                        //MatrixPrint(b);
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
                    b[j,0] /= cur_diag;
                }

                //MatrixPrint(matrix);
                //MatrixPrint(b);
                for (int k = height - 1; k > -1; k--)
                {
                    if (k != j)
                    {
                        double cur = matrix[k, j];
                        for (int i = 0; i < width; i++)
                            matrix[k, i] -= matrix[j, i] * cur;
                        b[k,0] -= b[j,0] * cur;
                        //MatrixPrint(matrix);
                        //MatrixPrint(b);
                    }
                }
            }
            return new Matrix(b).Transposed();
        }
    }
}
