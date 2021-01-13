using System;
namespace betablack
{
    public static class DoubleUtils
    {
        public static bool ApproximatelyEqual(double a, double b, double epsilon)
        {
            return Math.Abs(a - b) <= epsilon * Math.Max(Math.Abs(a), Math.Abs(b)) || Math.Abs(a - b) <= epsilon;
        }

        public static double DoubleEpsilon = 2.22044604925031E-16;
    }
}
