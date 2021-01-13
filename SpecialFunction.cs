using System;
using System.Runtime.CompilerServices;
namespace betablack
{
    public static class SpecialFunction
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double wapr(double x)
        {
            var em = -0.36787944117144233;
            var c23 = 0.66666666666666663;
            var e12 = 5.43656365691809;
            var x0 = 0.0012303916502879625;
            var x1 = -0.36766871970031223;
            var an3 = 2.6666666666666665;
            var an4 = 1.6265060240963856;
            var an5 = 4.2564102564102564;
            var an6 = 0.8923640462102;
            var s2 = 1.4142135623730951;
            var s21 = -0.17157287525380971;
            var s22 = -0.24264068711928566;
            var s23 = -0.58578643762690485;
            var xx = x;
            var delx = xx - em;
            var w = 0.0;
            var itflag = 0;
            if (Math.Abs(xx) <= x0)
            {
                w = xx / (1.0 + xx / (1.0 + xx / (2.0 + xx / (0.6 + 0.34 * xx))));
                itflag = 1;
            }
            else if (xx <= x1)
            {
                var reta = Math.Sqrt(e12 * delx);
                w = reta / (1 + reta / (3 + reta / (reta / (an4 + reta / (reta * an6 + an5)) + an3))) - 1.0;
                itflag = 1;
            }
            else if (xx <= 20)
            {
                // sqrt is ok
                var reta = s2 * Math.Sqrt(1 - xx / em);
                var an2 = 4.612634277343749 * Math.Sqrt(Math.Sqrt(reta + 1.09556884765625));
                w = reta / (1.0 + reta / (3.0 + (s21 * an2 + s22) * reta / (s23 * (an2 + reta)))) - 1.0;
            }
            else
            {
                return lambert_w_boost(xx);
                // bottle neck here. too many log, exp, pow evaluations ! use boost 
                // var zl = Math.Log(xx);
                // w = Math.Log(xx / Math.Log(xx / Math.Pow(zl, Math.Exp(-1.124491989777808 / (0.4225028202459761 + zl)))));
            }

            if (itflag == 0.0)
            {
                var zn = Math.Log(xx / w) - w;
                var temp = 1 + w;
                var temp2 = temp + c23 * zn;
                temp2 = 2 * temp * temp2;
                w = w * (1.0 + (zn / temp) * (temp2 - zn) / (temp2 - 2 * zn));
            }
            return w;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double wapr_b(double x, ref double x_b, double res_b)
        {
            var res = wapr(x);
            x_b = res_b * res / (x * (res + 1));
            return res_b;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double lambert_w_boost(double z)
        {
            if (z < 6.0) throw FS.FailWith("require(z > 6)");
            else if (z < 18)
            {
                // 6 < z < 18
                // Max error in interpolated form: 1.985e-19
                const double offset = 1.80937194824218750e+00;
                double[] P =
                {
                -1.80690935424793635e+00,
                -3.66995929380314602e+00,
                -1.93842957940149781e+00,
                -2.94269984375794040e-01,
                1.81224710627677778e-03,
                2.48166798603547447e-03,
         1.15806592415397245e-04,
         1.43105573216815533e-06,
         3.47281483428369604e-09
                };
                double[] Q = {
         1.00000000000000000e+00,
         2.57319080723908597e+00,
         1.96724528442680658e+00,
         5.84501352882650722e-01,
         7.37152837939206240e-02,
         3.97368430940416778e-03,
         8.54941838187085088e-05,
         6.05713225608426678e-07,
         8.17517283816615732e-10
                    };
                return offset + boost_math_tools_evaluate_rational(P, Q, z);
            }
            else if (z < 9897.12905874)  // 2.8 < log(z) < 9.2
            {
                // Max error in interpolated form: 1.195e-18
                double Y = -1.40297317504882812e+00;
                double[] P = {
         1.97011826279311924e+00,
         1.05639945701546704e+00,
         3.33434529073196304e-01,
         3.34619153200386816e-02,
         -5.36238353781326675e-03,
         -2.43901294871308604e-03,
         -2.13762095619085404e-04,
         -4.85531936495542274e-06,
         -2.02473518491905386e-08,
      };
                double[] Q = {
         1.00000000000000000e+00,
         8.60107275833921618e-01,
         4.10420467985504373e-01,
         1.18444884081994841e-01,
         2.16966505556021046e-02,
         2.24529766630769097e-03,
         9.82045090226437614e-05,
         1.36363515125489502e-06,
         3.44200749053237945e-09,
      };
                var log_w = Math.Log(z);
                return log_w + Y + boost_math_tools_evaluate_rational(P, Q, log_w);
            }
            else if (z < 7.896296e+13)  // 9.2 < log(z) <= 32
            {
                // Max error in interpolated form: 6.529e-18
                double Y = -2.73572921752929688e+00;
                double[] P = {
         3.30547638424076217e+00,
         1.64050071277550167e+00,
         4.57149576470736039e-01,
         4.03821227745424840e-02,
         -4.99664976882514362e-04,
         -1.28527893803052956e-04,
         -2.95470325373338738e-06,
         -1.76662025550202762e-08,
         -1.98721972463709290e-11,
      };
                double[] Q = {
         1.00000000000000000e+00,
         6.91472559412458759e-01,
         2.48154578891676774e-01,
         4.60893578284335263e-02,
         3.60207838982301946e-03,
         1.13001153242430471e-04,
         1.33690948263488455e-06,
         4.97253225968548872e-09,
         3.39460723731970550e-12,
      };
                var log_w = Math.Log(z);
                return log_w + Y + boost_math_tools_evaluate_rational(P, Q, log_w);
            }
            else if (z < 2.6881171e+43) // 32 < log(z) < 100
            {
                // Max error in interpolated form: 2.015e-18
                double Y = -4.01286315917968750e+00;
                double[] P = {
         5.07714858354309672e+00,
         -3.32994414518701458e+00,
         -8.61170416909864451e-01,
         -4.01139705309486142e-02,
         -1.85374201771834585e-04,
         1.08824145844270666e-05,
         1.17216905810452396e-07,
         2.97998248101385990e-10,
         1.42294856434176682e-13,
      };
                double[] Q = {
         1.00000000000000000e+00,
         -4.85840770639861485e-01,
         -3.18714850604827580e-01,
         -3.20966129264610534e-02,
         -1.06276178044267895e-03,
         -1.33597828642644955e-05,
         -6.27900905346219472e-08,
         -9.35271498075378319e-11,
         -2.60648331090076845e-14,
      };
                var log_w = Math.Log(z);
                return log_w + Y + boost_math_tools_evaluate_rational(P, Q, log_w);
            }
            else // 100 < log(z) < 710
            {
                // Max error in interpolated form: 5.277e-18
                double Y = -5.70115661621093750e+00;
                double[] P = {
         6.42275660145116698e+00,
         1.33047964073367945e+00,
         6.72008923401652816e-02,
         1.16444069958125895e-03,
         7.06966760237470501e-06,
         5.48974896149039165e-09,
         -7.00379652018853621e-11,
         -1.89247635913659556e-13,
         -1.55898770790170598e-16,
         -4.06109208815303157e-20,
         -2.21552699006496737e-24,
      };
                double[] Q = {
         1.00000000000000000e+00,
         3.34498588416632854e-01,
         2.51519862456384983e-02,
         6.81223810622416254e-04,
         7.94450897106903537e-06,
         4.30675039872881342e-08,
         1.10667669458467617e-10,
         1.31012240694192289e-13,
         6.53282047177727125e-17,
         1.11775518708172009e-20,
         3.78250395617836059e-25,
      };
                var log_w = Math.Log(z);
                return log_w + Y + boost_math_tools_evaluate_rational(P, Q, log_w);
            }
        }


        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double boost_math_tools_evaluate_rational(double[] P, double[] Q, double z)
        {
            var num = Polynomial.Evaluate(P, z);
            var den = Polynomial.Evaluate(Q, z);
            return num / den;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double pdf_normal(double x)
        {
            return 1.0 / Math.Sqrt(6.2831853071795862) * Math.Exp(-x * x / 2.0);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double cdf_normal(double x) // hart78
        {
            var num = Math.Abs(x);
            var num2 = 0.0;
            if (num < 37.0)
            {
                var num3 = Math.Exp(-0.5 * num * num);
                if (num < 7.07106781186547)
                {
                    var num4 = 0.0352624965998911 * num + 0.700383064443688;
                    num4 = num4 * num + 6.37396220353165;
                    num4 = num4 * num + 33.912866078383;
                    num4 = num4 * num + 112.079291497871;
                    num4 = num4 * num + 221.213596169931;
                    num4 = num4 * num + 220.206867912376;
                    num2 = num3 * num4;
                    num4 = 0.0883883476483184 * num + 1.75566716318264;
                    num4 = num4 * num + 16.064177579207;
                    num4 = num4 * num + 86.7807322029461;
                    num4 = num4 * num + 296.564248779674;
                    num4 = num4 * num + 637.333633378831;
                    num4 = num4 * num + 793.826512519948;
                    num4 = num4 * num + 440.413735824752;
                    num2 /= num4;
                }
                else
                {
                    var num4 = num + 0.65;
                    num4 = num + 4.0 / num4;
                    num4 = num + 3.0 / num4;
                    num4 = num + 2.0 / num4;
                    num4 = num + 1.0 / num4;
                    num2 = num3 / num4 / 2.506628274631;
                }
            }
            if (x > 0.0)
            {
                num2 = 1.0 - num2;
            }
            return num2;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]

        public static double pow2(double x)
        {
            return x * x;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double cdf_normal_b(double x, ref double x_b, double res_b)
        {
            var num = Math.Abs(x);
            var num2 = 0.0;
            if (num < 37.0)
            {
                var num3 = Math.Exp(-0.5 * num * num);
                if (num < 7.07106781186547)
                {
                    var num4 = 0.0352624965998911 * num + 0.700383064443688;
                    num4 = num4 * num + 6.37396220353165;
                    num4 = num4 * num + 33.912866078383;
                    num4 = num4 * num + 112.079291497871;
                    num4 = num4 * num + 221.213596169931;
                    num4 = num4 * num + 220.206867912376;
                    num2 = num3 * num4;
                    num4 = 0.0883883476483184 * num + 1.75566716318264;
                    num4 = num4 * num + 16.064177579207;
                    num4 = num4 * num + 86.7807322029461;
                    num4 = num4 * num + 296.564248779674;
                    num4 = num4 * num + 637.333633378831;
                    num4 = num4 * num + 793.826512519948;
                    num4 = num4 * num + 440.413735824752;
                    num2 /= num4;
                }
                else
                {
                    var num4 = num + 0.65;
                    num4 = num + 4.0 / num4;

                    num4 = num + 3.0 / num4;
                    num4 = num + 2.0 / num4;
                    num4 = num + 1.0 / num4;
                    num2 = num3 / num4 / 2.506628274631;
                }
                // This is regardless of the sign of x
                x_b += res_b * num3 / Math.Sqrt(6.2831853071795862);
            }
            if (x > 0.0)
            {
                num2 = 1.0 - num2;
            }
            return num2;
        }

        static readonly double[] a = {
        3.3871328727963666080,     1.3314166789178437745e+2,
        1.9715909503065514427e+3,  1.3731693765509461125e+4,
        4.5921953931549871457e+4,  6.7265770927008700853e+4,
        3.3430575583588128105e+4,  2.5090809287301226727e+3 };
        static readonly double[] b = {
        1.0,                       4.2313330701600911252e+1,
        6.8718700749205790830e+2,  5.3941960214247511077e+3,
        2.1213794301586595867e+4,  3.9307895800092710610e+4,
        2.8729085735721942674e+4,  5.2264952788528545610e+3 };
        static readonly double[] c = {
        1.42343711074968357734,     4.63033784615654529590,
        5.76949722146069140550,     3.64784832476320460504,
        1.27045825245236838258,     2.41780725177450611770e-1,
        2.27238449892691845833e-2,  7.74545014278341407640e-4 };
        static readonly double[] d = {
        1.0,                        2.05319162663775882187,
        1.67638483018380384940,     6.89767334985100004550e-1,
        1.48103976427480074590e-1,  1.51986665636164571966e-2,
        5.47593808499534494600e-4,  1.05075007164441684324e-9 };
        static readonly double[] e = {
        6.65790464350110377720,     5.46378491116411436990,
        1.78482653991729133580,     2.96560571828504891230e-1,
        2.65321895265761230930e-2,  1.24266094738807843860e-3,
        2.71155556874348757815e-5,  2.01033439929228813265e-7 };
        static readonly double[] f = {
        1.0,                        5.99832206555887937690e-1,
        1.36929880922735805310e-1,  1.48753612908506148525e-2,
        7.86869131145613259100e-4,  1.84631831751005468180e-5,
        1.42151175831644588870e-7,  2.04426310338993978564e-15 };

        static double R8_HUGE = 1e30;

        public static double cdfinv_normal(double p)

        /******************************************************************************/
        /*
          Purpose:

            R8_NORMAL_01_CDF_INVERSE inverts the standard normal CDF.

          Discussion:

            The result is accurate to about 1 part in 10^16.

          Licensing:

            This code is distributed under the GNU LGPL license. 

          Modified:

            19 March 2010

          Author:

            Original FORTRAN77 version by Michael Wichura.
            C version by John Burkardt.

          Reference:

            Michael Wichura,
            The Percentage Points of the Normal Distribution,
            Algorithm AS 241,
            Applied Statistics,
            Volume 37, Number 3, pages 477-484, 1988.

          Parameters:

            Input, double P, the value of the cumulative probability 
            densitity function.  0 < P < 1.  If P is outside this range, an "infinite"
            value is returned.

            Output, double R8_NORMAL_01_CDF_INVERSE, the normal deviate value 
            with the property that the probability of a standard normal deviate being 
            less than or equal to this value is P.
        */
        {
            var const1 = 0.180625;
            var const2 = 1.6;

            double q;
            double r;
            var split1 = 0.425;
            var split2 = 5.0;
            double value;

            if (p <= 0.0)
            {
                if (p < 0) throw FS.FailWith($"p is lower than 0, with {p}");
                value = -R8_HUGE;//r8_huge ( );
                return value;
            }

            if (1.0 <= p)
            {
                if (p > 1) throw FS.FailWith($"p is higher than 1 with {p}");
                value = R8_HUGE; // r8_huge ( );

                return value;
            }

            q = p - 0.5;

            //if ( r8_abs ( q ) <= split1 )
            if (Math.Abs(q) <= split1)
            {
                r = const1 - q * q;
                value = q * (((((((a[7] * r + a[6]) * r + a[5]) * r + a[4]) * r + a[3]) * r + a[2]) * r + a[1]) * r + a[0]) /
                            (((((((b[7] * r + b[6]) * r + b[5]) * r + b[4]) * r + b[3]) * r + b[2]) * r + b[1]) * r + b[0]);



            }
            else
            {
                if (q < 0.0)
                {
                    r = p;
                }
                else
                {
                    r = 1.0 - p;
                }

                if (r <= 0.0)
                {
                    value = -1.0;
                    throw FS.FailWith($"Error, cdfinv with p = {p}");
                    //exit ( 1 );// exception??
                }

                r = Math.Sqrt(-Math.Log(r));

                if (r <= split2)
                {
                    r = r - const2;
                    value = (((((((c[7] * r + c[6]) * r + c[5]) * r + c[4]) * r + c[3]) * r + c[2]) * r + c[1]) * r + c[0]) /
                              (((((((d[7] * r + d[6]) * r + d[5]) * r + d[4]) * r + d[3]) * r + d[2]) * r + d[1]) * r + d[0]);
                }
                else
                {
                    r = r - split2;
                    value = (((((((e[7] * r + e[6]) * r + e[5]) * r + e[4]) * r + e[3]) * r + e[2]) * r + e[1]) * r + e[0]) /
                             (((((((f[7] * r + f[6]) * r + f[5]) * r + f[4]) * r + f[3]) * r + f[2]) * r + f[1]) * r + f[0]);
                }

                if (q < 0.0)
                {
                    value = -value;
                }
            }
            return value;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double exp_x_N_y(double x, double y) // hart78
        {
            var num = Math.Abs(y);
            var num2 = 0.0;
            if (num < 37.0)
            {
                if (num < 7.07106781186547)
                {
                    var num4 = 0.0352624965998911 * num + 0.700383064443688;
                    num4 = num4 * num + 6.37396220353165;
                    num4 = num4 * num + 33.912866078383;
                    num4 = num4 * num + 112.079291497871;
                    num4 = num4 * num + 221.213596169931;
                    num4 = num4 * num + 220.206867912376;
                    num2 = num4;
                    num4 = 0.0883883476483184 * num + 1.75566716318264;
                    num4 = num4 * num + 16.064177579207;
                    num4 = num4 * num + 86.7807322029461;
                    num4 = num4 * num + 296.564248779674;
                    num4 = num4 * num + 637.333633378831;
                    num4 = num4 * num + 793.826512519948;
                    num4 = num4 * num + 440.413735824752;
                    num2 /= num4;
                }
                else
                {
                    var num4 = num + 0.65;
                    num4 = num + 4.0 / num4;

                    num4 = num + 3.0 / num4;
                    num4 = num + 2.0 / num4;
                    num4 = num + 1.0 / num4;
                    num2 = 1.0 / num4 / 2.506628274631;
                }
            }
            else // x > 37
            {
                //! N(-37) = 0.0 in most implementation, but for x negative we know the asymptotic behaviour.
                if (y < 0.0)
                {
                    var x2 = y * y;
                    var x4 = x2 * x2;
                    var x6 = x4 * x2;
                    var expansion = -1.0 / Math.Sqrt(2.0 * Math.PI) / y *
                        (1.0
                        - 1.0 / x2
                        + 3.0 / x4
                        - 15.0 / x6
                        //+ 105.0 / Math.Pow(x, 8)
                        //- 945.0 / Math.Pow(x, 10)
                        //+ 10395.0 / Math.Pow(x, 12)
                        //- 135135 / Math.Pow(x, 14)
                        //+ 2027025 / Math.Pow(x, 16)
                        //- 34459425 / Math.Pow(x, 18)
                        );
                    return Math.Exp(x - 0.5 * y * y) * expansion;
                }
                else // > 37, N = 1
                    return Math.Exp(x);
            }

            if (y > 0.0)
            {
                // no saving in this case
                return Math.Exp(x) * (1.0 - Math.Exp(-0.5 * num * num) * num2);
            }
            else return Math.Exp(x - 0.5 * num * num) * num2;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double exp_x_N_y_b(double x, double y, ref double x_b, ref double y_b, double res_b)
        {
            var num = Math.Abs(y);
            var num2 = 0.0;
            if (num < 37.0)
            {
                if (num < 7.07106781186547)
                {
                    var num4 = 0.0352624965998911 * num + 0.700383064443688;
                    num4 = num4 * num + 6.37396220353165;
                    num4 = num4 * num + 33.912866078383;
                    num4 = num4 * num + 112.079291497871;
                    num4 = num4 * num + 221.213596169931;
                    num4 = num4 * num + 220.206867912376;
                    num2 = num4;
                    num4 = 0.0883883476483184 * num + 1.75566716318264;
                    num4 = num4 * num + 16.064177579207;
                    num4 = num4 * num + 86.7807322029461;
                    num4 = num4 * num + 296.564248779674;
                    num4 = num4 * num + 637.333633378831;
                    num4 = num4 * num + 793.826512519948;
                    num4 = num4 * num + 440.413735824752;
                    num2 /= num4;
                }
                else
                {
                    var num4 = num + 0.65;
                    num4 = num + 4.0 / num4;

                    num4 = num + 3.0 / num4;
                    num4 = num + 2.0 / num4;
                    num4 = num + 1.0 / num4;
                    num2 = 1.0 / num4 / 2.506628274631;
                }
            }
            else // x > 37
            {
                //! N(-37) = 0.0 in most implementation, but for x negative we know the asymptotic behaviour.
                if (y < 0.0)
                {
                    var x2 = y * y;
                    var x4 = x2 * x2;
                    var x6 = x4 * x2;
                    var expansion = -1.0 / Math.Sqrt(2.0 * Math.PI) / y *
                        (1.0
                        - 1.0 / x2
                        + 3.0 / x4
                        - 15.0 / x6
                        //+ 105.0 / Math.Pow(x, 8)
                        //- 945.0 / Math.Pow(x, 10)
                        //+ 10395.0 / Math.Pow(x, 12)
                        //- 135135 / Math.Pow(x, 14)
                        //+ 2027025 / Math.Pow(x, 16)
                        //- 34459425 / Math.Pow(x, 18)
                        );

                    var tmp = Math.Exp(x - 0.5 * y * y);
                    var res = tmp * expansion;
                    x_b += res_b * res;
                    y_b += res_b * tmp / Math.Sqrt(6.2831853071795862);
                    return res;
                }
                else // > 37, N = 1
                {
                    var res = Math.Exp(x);
                    x_b += res_b * res;
                    return res;
                }
            }

            if (y > 0.0)
            {
                var e1 = Math.Exp(x);
                var e2 = Math.Exp(-0.5 * num * num);
                var res = e1 * (1.0 - e2 * num2);
                x_b += res_b * res;
                y_b += res_b * e1 * e2 / Math.Sqrt(6.2831853071795862);

                return res;
            }
            else
            {
                var tmp = Math.Exp(x - 0.5 * num * num);
                var res = tmp * num2;
                x_b += res_b * res;
                y_b += res_b * tmp / Math.Sqrt(6.2831853071795862);
                return res;
            }
        }

    }
}
