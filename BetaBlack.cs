using System;
using MathNet.Numerics.RootFinding;
using System.Runtime.CompilerServices;

namespace betablack
{
    
    public static class BetaBlack
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double g(double beta, double x) // beta >= 0 && beta < 1.0
        {
            if (beta == 0.0 || x >= 0.0) return Math.Exp(x);
            var oneMinustBeta = 1.0 - beta;
            return beta * oneMinustBeta + beta * (1 + x) + oneMinustBeta * oneMinustBeta * Math.Exp(x / oneMinustBeta);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double g_b(double beta, double x, ref double beta_b, ref double x_b, double res_b)
        {
            if (beta == 0.0 || x >= 0.0)
            {
                var res = Math.Exp(x);
                x_b += res_b * res;
                return res;
            }
            else
            {
                var oneMinustBeta = 1.0 - beta;
                var theExp = Math.Exp(x / oneMinustBeta);
                var res = beta * oneMinustBeta + beta * (1 + x) + oneMinustBeta * oneMinustBeta * theExp;
                beta_b += res_b * (theExp * (x + 2 * beta - 2) + x - 2 * beta + 2);
                x_b += res_b * (beta + oneMinustBeta * theExp);
                return res;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double inv_g(double beta, double y)
        {
            if (beta == 0.0 || y >= 1.0) return Math.Log(y);
            var expo = -(beta * beta - 2.0 * beta + y) / (beta * (beta - 1.0));

            if (expo > 700.0)
            {
                //! This happens only for small beta
                //! Need solution here - default to newton iteration? 
                //! For the time being we use an expansion
                var l1 = expo + Math.Log((1 - beta) / beta); // no risk here, beta is small
                var l12 = l1 * l1;
                var l13 = l12 * l1;
                var l14 = l12 * l12;
                var l15 = l14 * l12;
                var l2 = Math.Log(l1);
                var l22 = l2 * l2;
                var l23 = l22 * l2;
                var l24 = l23 * l2;
                var wwb = l1 - l2 + l2 / l1 + l2 * (-2.0 + l2) / (2.0 * l12)
                    + l2 * (6.0 - 9.0 * l2 + 2.0 * l22) / (6.0 * l13)
                    + l2 * (-12.0 + 36.0 * l2 - 22.0 * l22 + 3.0 * l23) / (12.0 * l14)
                    + l2 * (60.0 - 300.0 * l2 + 350.0 * l22 - 125.0 * l23 + 12.0 * l24) / (60.0 * l15);
                return beta * wwb - wwb + beta - 2.0 + y / beta;
            }

            var aw = Math.Exp(expo) * (1.0 - beta) / beta;
            var wb = SpecialFunction.wapr(aw);
            return beta * wb - wb + beta - 2.0 + y / beta;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double inv_g_b(double beta, double y, ref double beta_b, ref double y_b, double res_b)
        {
            //! This is neat
            var x = inv_g(beta, y);
            var beta_loc_b = 0.0;
            var x_loc_b = 0.0;
            g_b(beta, x, ref beta_loc_b, ref x_loc_b, 1.0);
            beta_b += res_b * (-beta_loc_b / x_loc_b);
            y_b += res_b * (1.0 / x_loc_b);
            return x;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double dg(double beta, double x)
        {
            if (beta == 0.0 || x >= 0.0) return Math.Exp(x);
            var oneMinusBeta = 1.0 - beta;
            return beta + oneMinusBeta * Math.Exp(x / oneMinusBeta);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double d2g(double beta, double x)
        {
            if (beta == 0.0 || x >= 0.0) return Math.Exp(x);
            return Math.Exp(x / (1.0 - beta));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double d3g(double beta, double x)
        {
            if (beta == 0.0 || x >= 0.0) return Math.Exp(x);
            var aa = 1.0 / (1.0 - beta);
            return aa * Math.Exp(aa * x);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double d4g(double beta, double x)
        {
            if (beta == 0.0 || x >= 0.0) return Math.Exp(x);
            var aa = 1.0 / (1.0 - beta);
            return aa * aa * Math.Exp(aa * x);
        }

        //! Density
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double beta_black_pdf(double beta, double mu, double sigma, double x)
        {
            var gm1 = inv_g(beta, x);
            var den = dg(beta, gm1);
            return SpecialFunction.pdf_normal((gm1 - mu) / sigma) / den / sigma;
        }

        //! Cdf
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double beta_black_cdf(double beta, double mu, double sigma, double x)
        {
            return SpecialFunction.cdf_normal((inv_g(beta, x) - mu) / sigma);
        }

        //! Quantile
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double beta_black_cdfinv(double beta, double mu, double sigma, double p)
        {
            return g(beta, SpecialFunction.cdfinv_normal(p) * sigma + mu);
        }

        //! Next one to be differentiated
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double beta_black_trunc(double beta, double mu, double sigma, double x)
        {
            var betaMinusOne = beta - 1.0;
            var betaMinusOneSq = betaMinusOne * betaMinusOne;
            var sigma2 = sigma * sigma;
            var sqrt2 = Math.Sqrt(2.0);
            var sqrtPi = Math.Sqrt(Math.PI);
            var t = -mu / sigma;
            var one_over_two_betaMinusOneSq = 0.5 / betaMinusOneSq;
            var f0 = -1.0 / (2.0 * sqrtPi);
            if (x < t)
            {
                var x2 = x * x;
                var e2 = Math.Exp(-0.5 * x2);
                var a4 = x + sigma / betaMinusOne;
                var group4 = SpecialFunction.exp_x_N_y(-beta * mu / betaMinusOneSq + (sigma2 + 2.0 * mu) * one_over_two_betaMinusOneSq, a4);
                return f0 * sqrt2 * e2 * beta * sigma +
                        +betaMinusOneSq * group4
                        + beta * SpecialFunction.cdf_normal(x) * (mu - beta + 2.0);
            }
            else
            {
                var t2 = t * t;
                var e4 = Math.Exp(-0.5 * t2);
                var a1 = t + sigma / betaMinusOne;
                var group1 = -SpecialFunction.exp_x_N_y(-beta * mu / betaMinusOneSq + (sigma2 + 2.0 * mu) * one_over_two_betaMinusOneSq, a1);
                var a2 = sigma - t;
                var a3 = x - sigma;
                var group2 = -SpecialFunction.exp_x_N_y(mu + 0.5 * sigma2, -a2);
                var group3 = SpecialFunction.exp_x_N_y(mu + 0.5 * sigma2, a3);
                return f0 * sqrt2 * e4 * beta * sigma
                    - betaMinusOneSq * group1
                    + (group2 + group3)
                    - beta * (mu - beta + 2) * (SpecialFunction.cdf_normal(-t) - 1.0);
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double beta_black_trunc_b(double beta, double mu, double sigma, double x, ref double beta_b, ref double mu_b, ref double sigma_b, ref double x_b, double res_b)
        {
            var betaMinusOne = beta - 1.0;
            var betaMinusOneSq = betaMinusOne * betaMinusOne;
            var sigma2 = sigma * sigma;
            var sqrt2 = Math.Sqrt(2.0);
            var sqrtPi = Math.Sqrt(Math.PI);
            var t = -mu / sigma;
            var one_over_two_betaMinusOneSq = 0.5 / betaMinusOneSq;
            var f0 = -1.0 / (2.0 * sqrtPi);
            var res = 0.0;
            var betaMinusOne_b = 0.0;
            var betaMinusOneSq_b = 0.0;
            var one_over_two_betaMinusOneSq_b = 0.0;
            var x2_b = 0.0;
            var t_b = 0.0;
            if (x < t)
            {
                var x2 = x * x;
                var ar2 = -0.5 * x2;
                var e2 = Math.Exp(ar2);
                var a4 = x + sigma / betaMinusOne;
                var ar4 = -beta * mu / betaMinusOneSq + (sigma2 + 2.0 * mu) * one_over_two_betaMinusOneSq;
                var dgroup4_dar4 = 0.0;
                var dgroup4_da4 = 0.0;
                var group4 = SpecialFunction.exp_x_N_y_b(ar4, a4, ref dgroup4_dar4, ref dgroup4_da4, 1.0);
                var dNx_dx = 0.0;
                var Nx = SpecialFunction.cdf_normal_b(x, ref dNx_dx, 1.0);
                res = f0 * sqrt2 * e2 * beta * sigma
                            + betaMinusOneSq * group4
                            + beta * Nx * (mu - beta + 2.0);
                var e2_b = res_b * f0 * sqrt2 * beta * sigma;
                beta_b += res_b * (f0 * sqrt2 * e2 * sigma + Nx * (mu - 2.0 * beta + 2.0));
                sigma_b += res_b * f0 * sqrt2 * e2 * beta;
                betaMinusOneSq_b += res_b * group4;
                var group4_b = res_b * betaMinusOneSq;
                var Nx_b = res_b * beta * (mu - beta + 2.0);
                mu_b += res_b * beta * Nx;
                x_b += Nx_b * dNx_dx;
                var ar4_b = group4_b * dgroup4_dar4;
                var a4_b = group4_b * dgroup4_da4;
                beta_b += ar4_b * (-mu / betaMinusOneSq);
                mu_b += ar4_b * (-beta / betaMinusOneSq + 2.0 * one_over_two_betaMinusOneSq);
                betaMinusOneSq_b += ar4_b * (beta * mu / betaMinusOneSq / betaMinusOneSq);
                sigma_b += ar4_b * 2.0 * sigma * one_over_two_betaMinusOneSq;
                one_over_two_betaMinusOneSq_b += ar4_b * (sigma2 + 2.0 * mu);
                x_b += a4_b;
                sigma_b += a4_b / betaMinusOne;
                betaMinusOne_b += a4_b * (-sigma / betaMinusOne / betaMinusOne);
                var ar2_b = e2_b * e2;
                x2_b += -0.5 * ar2_b;
                x_b += x2_b * 2.0 * x;
            }
            else
            {
                var t2 = t * t;
                var ar4 = -0.5 * t2;
                var e4 = Math.Exp(ar4);
                var a1 = t + sigma / betaMinusOne;
                var ar1 = -beta * mu / betaMinusOneSq + (sigma2 + 2.0 * mu) * one_over_two_betaMinusOneSq;
                var dmgroup1_dar1 = 0.0;
                var dmgroup1_da1 = 0.0;
                var group1 = -SpecialFunction.exp_x_N_y_b(ar1, a1, ref dmgroup1_dar1, ref dmgroup1_da1, 1.0);
                var a2 = t - sigma;
                var a3 = x - sigma;
                var ar2 = mu + 0.5 * sigma2;
                var dmgroup2_dar2 = 0.0;
                var dmgroup2_da2 = 0.0;
                var group2 = -SpecialFunction.exp_x_N_y_b(ar2, a2, ref dmgroup2_dar2, ref dmgroup2_da2, 1.0);
                var dgroup3_dar2 = 0.0;
                var dgroup3_da3 = 0.0;
                var group3 = SpecialFunction.exp_x_N_y_b(ar2, a3, ref dgroup3_dar2, ref dgroup3_da3, 1.0);
                var dnmt_dmt = 0.0;
                var nmt = SpecialFunction.cdf_normal_b(-t, ref dnmt_dmt, 1.0);
                res = f0 * sqrt2 * e4 * beta * sigma
                            - betaMinusOneSq * group1
                            + group2 + group3
                            - beta * (mu - beta + 2) * (nmt - 1.0);
                var e4_b = res_b * f0 * sqrt2 * beta * sigma;
                beta_b += res_b * (f0 * sqrt2 * e4 * sigma - (mu - 2.0 * beta + 2.0) * (nmt - 1.0));
                sigma_b += res_b * f0 * sqrt2 * e4 * beta;
                betaMinusOneSq_b += -res_b * group1;
                var group1_b = -res_b * betaMinusOneSq;
                var group2_b = res_b;
                var group3_b = res_b;
                mu_b += -res_b * beta * (nmt - 1.0);
                var nmt_b = res_b * -beta * (mu - beta + 2);
                t_b += -nmt_b * dnmt_dmt;
                var ar2_b = group3_b * dgroup3_dar2;
                var a3_b = group3_b * dgroup3_da3;
                var a2_b = -group2_b * dmgroup2_da2;
                ar2_b += -group2_b * dmgroup2_dar2;
                mu_b += ar2_b;
                sigma_b += ar2_b * sigma;
                x_b += a3_b;
                sigma_b += -a3_b;
                sigma_b += -a2_b;
                t_b += a2_b;
                var ar1_b = -group1_b * dmgroup1_dar1;
                var a1_b = -group1_b * dmgroup1_da1;
                beta_b += ar1_b * (-mu / betaMinusOneSq);
                mu_b += ar1_b * (-beta / betaMinusOneSq + 2.0 * one_over_two_betaMinusOneSq);
                betaMinusOneSq_b += ar1_b * (beta * mu / betaMinusOneSq / betaMinusOneSq);
                sigma_b += ar1_b * 2.0 * sigma * one_over_two_betaMinusOneSq;
                one_over_two_betaMinusOneSq_b += ar1_b * (sigma2 + 2.0 * mu);
                t_b += a1_b;
                sigma_b += a1_b / betaMinusOne;
                betaMinusOne_b += -a1_b * sigma / betaMinusOne / betaMinusOne;
                var ar4_b = e4_b * e4;
                var t2_b = -0.5 * ar4_b;
                t_b += t2_b * 2.0 * t;
            }
            betaMinusOneSq_b += -one_over_two_betaMinusOneSq_b * 0.5 / betaMinusOneSq / betaMinusOneSq;
            mu_b += -t_b / sigma;
            sigma_b += t_b * mu / sigma / sigma;
            betaMinusOne_b += betaMinusOneSq_b * 2.0 * betaMinusOne;
            beta_b += betaMinusOne_b;
            return res;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double m1_beta_black(double beta, double mu, double sigma)
        {
            var betaMinusOne = beta - 1.0;
            var betaMinusOneSq = betaMinusOne * betaMinusOne;
            var sigma2 = sigma * sigma;
            var sqrt2 = Math.Sqrt(2.0);
            var sqrtPi = Math.Sqrt(Math.PI);
            var t = -mu / sigma;
            var t2 = t * t;
            var one_over_two_betaMinusOneSq = 0.5 / betaMinusOneSq;
            var e4 = Math.Exp(-0.5 * t2);
            var a1 = t + sigma / betaMinusOne;
            var group1 = -SpecialFunction.exp_x_N_y(-beta * mu / betaMinusOneSq + (sigma2 + 2.0 * mu) * one_over_two_betaMinusOneSq, a1);
            var a2 = sigma - t;
            var group2 = SpecialFunction.exp_x_N_y(mu + 0.5 * sigma2, a2);
            var f0 = -1.0 / (2.0 * sqrtPi);
            return f0 * sqrt2 * e4 * beta * sigma
                - betaMinusOneSq * group1
                + group2
                - beta * (mu - beta + 2) * (SpecialFunction.cdf_normal(-t) - 1.0);
        }

        // Checked.
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double m1_beta_black_b(double beta, double mu, double sigma, ref double beta_b, ref double mu_b, ref double sigma_b, double res_b)
        {
            var betaMinusOne = beta - 1.0;
            var betaMinusOneSq = betaMinusOne * betaMinusOne;
            var sigma2 = sigma * sigma;
            var sqrt2 = Math.Sqrt(2.0);
            var sqrtPi = Math.Sqrt(Math.PI);
            var t = -mu / sigma;
            var t2 = t * t;
            var one_over_two_betaMinusOneSq = 0.5 / betaMinusOneSq;
            var e4 = Math.Exp(-0.5 * t2);
            var a1 = t + sigma / betaMinusOne;
            var b1 = -beta * mu / betaMinusOneSq + (sigma2 + 2.0 * mu) * one_over_two_betaMinusOneSq;
            var dmgroup1_db1 = 0.0;
            var dmgroup1_da1 = 0.0;
            var group1 = -SpecialFunction.exp_x_N_y_b(b1, a1, ref dmgroup1_db1, ref dmgroup1_da1, 1.0);
            var a2 = sigma - t;
            var muOneHalfS2 = mu + 0.5 * sigma2;
            var dgroup2_dmuOneHalfS2 = 0.0;
            var dgroup2_da2 = 0.0;
            var group2 = SpecialFunction.exp_x_N_y_b(muOneHalfS2, a2, ref dgroup2_dmuOneHalfS2, ref dgroup2_da2, 1.0);
            var f0 = -1.0 / (2.0 * sqrtPi);
            var dnmt_dmt = 0.0;
            var nmt = SpecialFunction.cdf_normal_b(-t, ref dnmt_dmt, 1.0);
            var res = f0 * sqrt2 * e4 * beta * sigma
                        - betaMinusOneSq * group1
                        + group2
                        - beta * (mu - beta + 2) * (nmt - 1.0);
            var e4_b = res_b * f0 * sqrt2 * beta * sigma;
            beta_b += res_b * f0 * sqrt2 * e4 * sigma;
            beta_b += -res_b * (mu - 2.0 * beta + 2.0) * (nmt - 1.0);
            sigma_b += res_b * f0 * sqrt2 * e4 * beta;
            var betaMinusOneSq_b = -res_b * group1;
            var group1_b = -res_b * betaMinusOneSq;
            var group2_b = res_b;
            mu_b += -res_b * beta * (nmt - 1.0);
            var nmt_b = -res_b * beta * (mu - beta + 2);
            var t_b = -nmt_b * dnmt_dmt;
            var a2_b = group2_b * dgroup2_da2;
            var muOneHalfS2_b = group2_b * dgroup2_dmuOneHalfS2;
            mu_b += muOneHalfS2_b;
            sigma_b += muOneHalfS2_b * sigma;
            sigma_b += a2_b;
            t_b += -a2_b;
            var b1_b = -group1_b * dmgroup1_db1;
            var a1_b = -group1_b * dmgroup1_da1;
            beta_b += -b1_b * mu / betaMinusOneSq;
            mu_b += b1_b * (-beta / betaMinusOneSq + 2.0 * one_over_two_betaMinusOneSq);
            sigma_b += b1_b * 2.0 * sigma * one_over_two_betaMinusOneSq;
            betaMinusOneSq_b += b1_b * beta * mu / betaMinusOneSq / betaMinusOneSq;
            var one_over_two_betaMinusOneSq_b = b1_b * (sigma2 + 2.0 * mu);
            t_b += a1_b;
            var betaMinusOne_b = -a1_b * sigma / betaMinusOne / betaMinusOne;
            sigma_b += a1_b / betaMinusOne;
            var t2_b = -0.5 * e4_b * e4;
            betaMinusOneSq_b += -one_over_two_betaMinusOneSq_b * 0.5 / betaMinusOneSq / betaMinusOneSq;
            t_b += t2_b * 2.0 * t;
            mu_b += -t_b / sigma;
            sigma_b += t_b * mu / sigma / sigma;
            betaMinusOne_b += betaMinusOneSq_b * 2.0 * betaMinusOne;
            beta_b += betaMinusOne_b;
            return res;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double inv_m1_beta_black_convexity(double beta, double f, double sigma)
        {
            var s2 = sigma * sigma;
            var alpha = 0.5 * beta * s2;
            var cvx = alpha == 0.0 ? // rule of thumb for convexity decay
                            1.0 : Math.Min(1.0, Math.Max(0.0, (f + alpha) / alpha));
            return inv_g(beta, f) - 0.5 * s2 * cvx;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double inv_m1_beta_black_cxv(double beta, double f, double sigma)
        {
            var s2 = sigma * sigma;
            var igb = inv_g(beta, f);
            var fac = -0.5 * s2;
            var rf = FS.fun((double cxv) =>
            {
                var mu = igb + fac * cxv;
                var f0 = m1_beta_black(beta, mu, sigma);
                return f - f0;
            });

            var alpha = 0.5 * beta * s2;
            var guess = alpha == 0.0 ? 1.0 : Math.Min(1.0, Math.Max(0.0, (f + alpha) / alpha));
            var lo = guess - 1e-3;
            var hi = guess + 1e-3;

            var tGuess = rf(alpha);
            if (DoubleUtils.ApproximatelyEqual(tGuess, 0.0, 100.0 * DoubleUtils.DoubleEpsilon))
                return igb + fac * alpha;

            NumericalRecipes.zbrac(rf, ref lo, ref hi);
            var cxvRoot = Brent.FindRoot(rf, lo, hi);
            return igb + fac * cxvRoot;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double inv_m1_beta_black(double beta, double f, double sigma)
        {
            return inv_m1_beta_black_newton(beta, f, sigma);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double inv_m1_beta_black_newton(double beta, double f, double sigma)
        {
            //! Quite accurate but does not fix large sigma issue
            //! Maybe second order correction?
            var mu0 = inv_m1_beta_black_convexity(beta, f, sigma);

            var muMax = 37.0 - 0.5 * sigma * sigma; //  inv_g(beta, f); 
                                                    //! Note, can have bespoke function that only does mu derivative.
                                                    //! Forward mode in this case will be faster.
            var beta_b = 0.0;
            var mu_b = 0.0;
            var sigma_b = 0.0;
            var f0 = m1_beta_black_b(beta, mu0, sigma, ref beta_b, ref mu_b, ref sigma_b, 1.0);
            var maxIt = 20;
            while (maxIt > 0 && !DoubleUtils.ApproximatelyEqual(f0, f, 1.0 * Math.Sqrt(DoubleUtils.DoubleEpsilon)))
            {
                --maxIt;
                mu0 -= (f0 - f) / mu_b;
                if (mu0 > muMax)
                {
                    // over shoot, default to brent, probably large vol
                    return inv_m1_beta_black_brent(beta, f, sigma);
                }
                mu_b = 0.0;
                f0 = m1_beta_black_b(beta, mu0, sigma, ref beta_b, ref mu_b, ref sigma_b, 1.0);
            }
            if (maxIt == 0)
            {
                return inv_m1_beta_black_brent(beta, f, sigma);
            }
            mu0 -= (f0 - f) / mu_b; //! perform an extra iteration
            return mu0;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double inv_m1_beta_black_brent(double beta, double f, double sigma)
        {
            var guess = inv_m1_beta_black_convexity(beta, f, sigma);
            var shift = 1.5e-8 * (Math.Max(1.0, Math.Abs(guess)));
            var lo = guess - shift;
            var hi = guess + shift;
            var rf = FS.fun((double mu) =>
            {
                var rr = m1_beta_black(beta, mu, sigma);
                return rr - f;
            });
            return Brent.FindRootExpand(rf, lo, hi);
        }

        // checked
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double inv_m1_beta_black_b(double beta, double f, double sigma, ref double beta_b, ref double f_b, ref double sigma_b, double res_b)
        {
            // very neat
            var res = inv_m1_beta_black(beta, f, sigma);
            var beta_b_loc = 0.0;
            var mu_b_loc = 0.0;
            var sigma_b_loc = 0.0;
            m1_beta_black_b(beta, res, sigma, ref beta_b_loc, ref mu_b_loc, ref sigma_b_loc, 1.0);
            f_b += res_b / mu_b_loc;
            beta_b += res_b * (-beta_b_loc / mu_b_loc);
            sigma_b += res_b * (-sigma_b_loc / mu_b_loc);
            return res;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double beta_black(double beta, OptionType optType, double F, double K, double T, double sigma)
        {
            if (T < 0.0) throw FS.FailWith($"T < 0 not allowed, {T}");
            if (sigma < 0.0) throw FS.FailWith($"sigma < 0 not allowed {sigma}");
            var stdev = sigma * Math.Sqrt(T);
            if (stdev == 0.0) return Math.Max(optType.ToEpsilon() * (F - K), 0.0);
            var mu = inv_m1_beta_black(beta, F, stdev);
            var d = (inv_g(beta, K) - mu) / stdev;
            var t1 = K * SpecialFunction.cdf_normal(d);
            var t2 = beta_black_trunc(beta, mu, stdev, d);
            var put = t1 - t2;

            var intrinsic = Math.Max(optType.ToEpsilon() * (F - K), 0.0);
            if (optType == OptionType.Put) return Math.Max(put, intrinsic);
            else if (optType == OptionType.Call) return Math.Max(put + F - K, intrinsic);
            else throw FS.E_CASE(optType);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double beta_black_normal_vol(double beta, double F, double K, double T, double sigma)
        {
            var optType = OptionType.Put; // F > K ? OptionType.Put : OptionType.Call;
            var bbPrice = beta_black(beta, optType, F, K, T, sigma);

            try
            {

                return ClosedForm.ibachelier(optType, F, K, T, bbPrice);
            }
            catch (Exception e)
            {
                FS.WithInfo(e, $"beta ={beta}, optType = {optType}, F = {F}, K = {K}, T = {T}, sigma = {sigma}, bbPrice = {bbPrice}");
                throw;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double beta_black_b(double beta, OptionType optType, double F, double K, double T, double sigma
            , ref double beta_b, ref double F_b, ref double K_b, ref double T_b, ref double sigma_b, double res_b)
        {
            if (T < 0.0) throw FS.FailWith($"T < 0 not allowed, {T}");
            if (sigma < 0.0) throw FS.FailWith($"sigma < 0 not allowed {sigma}");
            var sqrtT = Math.Sqrt(T);
            var stdev = sigma * sqrtT;
            if (stdev == 0.0) return Math.Max(optType.ToEpsilon() * (F - K), 0.0);
            var dmu_dbeta = 0.0;
            var dmu_dF = 0.0;
            var dmu_dstdev = 0.0;
            var mu = inv_m1_beta_black_b(beta, F, stdev, ref dmu_dbeta, ref dmu_dF, ref dmu_dstdev, 1.0);
            var digb_dbeta = 0.0;
            var digb_dK = 0.0;
            var igb = inv_g_b(beta, K, ref digb_dbeta, ref digb_dK, 1.0);
            var d = (igb - mu) / stdev;
            var dNd_dd = 0.0;
            var Nd = SpecialFunction.cdf_normal_b(d, ref dNd_dd, 1.0);
            var t1 = K * Nd;
            var dt2_dbeta = 0.0;
            var dt2_dmu = 0.0;
            var dt2_dstdev = 0.0;
            var dt2_dd = 0.0;
            var t2 = beta_black_trunc_b(beta, mu, stdev, d, ref dt2_dbeta, ref dt2_dmu, ref dt2_dstdev, ref dt2_dd, 1.0);
            var res = t1 - t2;
            // now rewind
            var t1_b = res_b;
            var t2_b = -res_b;
            beta_b += t2_b * dt2_dbeta;
            var mu_b = t2_b * dt2_dmu;
            var stdev_b = t2_b * dt2_dstdev;
            var d_b = t2_b * dt2_dd;
            K_b += t1_b * Nd;
            var Nd_b = t1_b * K;
            d_b += Nd_b * dNd_dd;
            var igb_b = d_b / stdev;
            mu_b += -d_b / stdev;
            stdev_b += -d_b * (igb - mu) / stdev / stdev;
            beta_b += igb_b * digb_dbeta;
            K_b += igb_b * digb_dK;
            beta_b += mu_b * dmu_dbeta;
            F_b += mu_b * dmu_dF;
            stdev_b += mu_b * dmu_dstdev;
            sigma_b += stdev_b * sqrtT;
            T_b += stdev_b * sigma / (2.0 * sqrtT);
            if (optType == OptionType.Put) return res;
            else if (optType == OptionType.Call)
            {
                F_b += res_b;
                K_b += -res_b;
                return res + F - K;
            }
            else throw FS.E_CASE(optType);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double beta_black_quick_delta_to_strike(double beta, double F, double T, double sigma, double qd)
        {
            if (qd < -1.0 || qd > 1.0) throw FS.FailWith($"delta must be in (-1, 1), got {qd}");
            var qdUsed = qd < 0.0 ? 1.0 + qd : qd;
            return g(beta, inv_g(beta, F) - SpecialFunction.cdfinv_normal(qdUsed) * sigma * Math.Sqrt(T));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double beta_black_strike_to_quick_delta(double beta, double F, double T, double sigma, double strike)
        {
            var qd = SpecialFunction.cdf_normal((inv_g(beta, F) - inv_g(beta, strike)) / (sigma * Math.Sqrt(T)));
            if (qd <= -0.5) return 1.0 + qd;
            if (qd > 0.5) return qd - 1.0;
            return qd;
        }

        public static double beta_black_delta_to_strike(double beta, double F, double T, double sigma, double delta)
        {
            if (delta < -1.0 || delta > 1.0) throw FS.FailWith($"delta must be in (-1, 1), got {delta}");
            var mu = inv_m1_beta_black(beta, F, sigma * Math.Sqrt(T));
            var usedDelta = delta < 0.0 ? -delta : 1.0 - delta;
            return g(beta, SpecialFunction.cdfinv_normal(usedDelta) * sigma * Math.Sqrt(T) + mu);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double ibeta_black(double beta, OptionType optType, double F, double K, double T, double premium)
        {
            var rf = FS.fun((double sigma) => beta_black(beta, optType, F, K, T, sigma) - premium);

            var lo = 1.0;
            var hi = 5.0; // If very deep in-out of the money, will fail. f(lo) = f(hi), zbrac only moves low
            NumericalRecipes.zbrac(rf, ref lo, ref hi);
            return Brent.FindRoot(rf, lo, hi);
        }

        ///// EXPERIMENTAL
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double m1_beta_black_approx_2(double beta, double mu, double sigma)
        {
            return g(beta, mu) + 0.5 * sigma * sigma * d2g(beta, mu);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double inv_m1_beta_black_approx_2(double beta, double f, double sigma)
        {
            //! Numerical 
            var expArg = -1.0 * (SpecialFunction.pow2(beta) - 2.0 * beta + f) / ((-1.0 + beta) * beta);
            if (expArg > 37.0) return inv_g(beta, f);
            var aw = -0.5000000000 * (2.0 * SpecialFunction.pow2(beta) + SpecialFunction.pow2(sigma) - 4.0 * beta + 2.0)
                * Math.Exp(expArg)
                / ((-1.0 + beta) * beta);
            var wb = SpecialFunction.wapr(aw);
            return beta * wb - wb + beta - 2.0 + f / beta;
        }

        //! This should provide a very accurate guess
        //! However, not robust to high volatility..
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double inv_m1_approx_smart(double beta, double f, double sigma)
        {
            // Compute some relevant bounds
            var log_001 = -4.605170186;
            var lo = log_001 * (1 - beta) - 0.5 * sigma * sigma;
            var hi = Math.Min(37 - 0.5 * sigma * sigma, 0.5 * sigma * sigma); //! < 10
                                                                                // Evaluate the function on those points
                                                                                // Using the actual function to too costly 
            var zz = 0.0;
            var vlo = m1_beta_black(beta, lo, sigma);
            var dhi_dmu = 0.0;
            var vhi = m1_beta_black_b(beta, hi, sigma, ref zz, ref dhi_dmu, ref zz, 1.0);
            // Treat bounds
            if (f < vlo) return inv_g(beta, f);  // no convexity
            if (f > vhi) return inv_g(beta, f) - 0.5 * sigma * sigma;    // full convexity
                                                                            //! Middle region, fit approximately
                                                                            // a + b x + c exp(d x)
            var d = dhi_dmu / vhi;
            var c = vhi / Math.Exp(dhi_dmu * hi / vhi);
            var b = -(vlo - c * Math.Exp(d * lo)) / (hi - lo);
            var a = -b * hi;
            //! The invert in closed form
            var aw = Math.Exp(-d * (a - f) / b) * c * d / b;
            return -(b * SpecialFunction.wapr(aw) + a * d - d * f) / (b * d);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double beta_black_break_even(double beta, double f, double t, double sigma, double dt = 1.0 / 365.0)
        {
            try
            {
                //var pack = new OptionBetaBlack(beta, OptionType.Put, f, f, t, sigma);
                //var gamma = pack.d2PV_dS2;
                //if (gamma == 0.0) throw FS.FailWith($"Cannot compute breakeven, gamma is zero.");
                //var res = Math.Sqrt(2.0 * pack.dPV_dT * dt / pack.d2PV_dS2);

                var zz = 0.0;
                var T_b = 0.0;
                var f_b = 0.0;
                beta_black_b(beta, OptionType.Put, f, f, t, sigma, ref zz, ref f_b, ref zz, ref T_b, ref zz, 1.0);

                var df = Math.Sqrt(DoubleUtils.DoubleEpsilon) * (1.0 + Math.Abs(f));
                var fup_b = 0.0;
                beta_black_b(beta, OptionType.Put, f + df, f, t, sigma, ref zz, ref fup_b, ref zz, ref zz, ref zz, 1.0);

                var gamma = (fup_b - f_b) / df;
                if (gamma == 0.0) throw FS.FailWith($"Cannot compute Break even, gamma is zero");
                var res = Math.Sqrt(2.0 * T_b * dt / gamma);
                // if (double.IsNaN(res)) throw FS.FailWith("Nan");
                return res;
            }
            catch (Exception e)
            {
                FS.WithInfo(e, $"beta ={beta}, f = {f}, t = {t}, sigma = {sigma}, dt = {dt}");
                throw;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double beta_black_break_even_to_sigma(double beta, double f, double t, double breakEven, double dt = 1.0 / 365.0)
        {
            var rf = FS.fun((double sigma) => beta_black_break_even(beta, f, t, sigma, dt) - breakEven);
            var x1 = 1.0;
            var x2 = 5.0;
            NumericalRecipes.zbrac(rf, ref x1, ref x2);
            return Brent.FindRoot(rf, x1, x2);
        }
    }
}
