using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace betablack
{
    public static class Polynomial
    {
        public static double Evaluate(double[] coeffs, double x)
        {
            var num = 0.0;
            for (var i = coeffs.Length - 1; i >= 0; i--)
                num = coeffs[i] + num * x;
            return num;
        }
    }
}
