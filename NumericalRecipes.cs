using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace betablack
{
    public static class NumericalRecipes
    {
        public static bool zbrac(Func<double, double> f, ref double x1, ref double x2, double num2 = 1.6)
        {
            var num = 50;
            if (x1 == x2)
                throw FS.FailWith($"Bad initial range in zbrac x1 = x2 with {x1}");
            var num3 = 0;
            var num4 = f(x1);
            var num5 = f(x2);            
            while (num3++ < num)
            {
                if (num4 * num5 <= 0.0)
                {
                    return true;
                }
                
                if (Math.Abs(num4) < Math.Abs(num5))
                {
                    x1 /= num2;
                    num4 = f(x1);
                }
                else
                {
                    x2 *= num2;
                    num5 = f(x2);
                }
            }
            return false;
        }
    }
}
