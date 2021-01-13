using System;
using System.Runtime.CompilerServices;

namespace betablack
{
    public static class ClosedForm
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double ibachelier(OptionType optType, double F, double K, double T, double premium, bool safe = false)
        {
            if (T <= 0.0) throw FS.FailWith($"T must be postive, got {T}");
            var epsilon = optType.ToEpsilon();
            var intrinsic = Math.Max(epsilon * (F - K), 0.0);
            if (premium <= intrinsic)
            {
                if (DoubleUtils.ApproximatelyEqual(premium, intrinsic, DoubleUtils.DoubleEpsilon))
                    return 0.0;


                if (safe)
                    return 0.0;
                throw FS.FailWith($"premium must be larger than intrinsic. Got premium = {premium},  vs intrinsic = {intrinsic}"); //! 
            }

            if (F == K) return premium * Math.Sqrt(2.0 * Math.PI / T);

            var sqrtT = Math.Sqrt(T);
            var phiStar = -Math.Abs((premium - intrinsic) / (K - F));

            var xBar = 0.0;
            if (phiStar < -0.001882039271)
            {
                var g = 1.0 / (phiStar - 0.5);
                var g2 = g * g;
                var zeta = (0.032114372355 - g2 * (0.016969777977 - g2 * (2.6207332461e-3 - 9.6066952861e-5 * g2)))
                            / (1.0 - g2 * (0.6635646938 - g2 * (0.14528712196 - 0.010472855461 * g2)));
                xBar = g * (1.0 / Math.Sqrt(2.0 * Math.PI) + zeta * g2);
            }
            else
            {
                var h = Math.Sqrt(-Math.Log(-phiStar));
                xBar = (9.4883409779 - h * (9.6320903635 - h * (0.58556997323 + 2.1464093351 * h)))
                          / (1 - h * (0.65174820867 + h * (1.5120247828 + 6.6437847132e-5 * h)));
            }

            var xBar2 = xBar * xBar;
            var nxBar = SpecialFunction.pdf_normal(xBar);
            var phiTildaXBar = SpecialFunction.cdf_normal(xBar) + nxBar / xBar;
            var q = (phiTildaXBar - phiStar) / nxBar;
            var xStar = xBar + (3.0 * q * xBar2 * (2.0 - q * xBar * (2.0 + xBar2)))
                            / (6.0 + q * xBar * (-12.0 + xBar * (6.0 * q + xBar * (-6.0 + q * xBar * (3.0 + xBar2)))));
            return Math.Abs((K - F) / (xStar * sqrtT));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double black76(OptionType optType, double F, double K, double sdev)
        {
            var num = optType.ToEpsilon();
            if (sdev < 0.0)
                throw FS.FailWith($"sdev must be positive. Got {sdev}");
            if (sdev == 0.0 || K <= 0.0)
                return Math.Max(num * (F - K), 0.0);
            else
            {
                var num2 = (Math.Log(F / K) + 0.5 * sdev * sdev) / sdev;
                var num3 = num2 - sdev;
                return num * (F * SpecialFunction.cdf_normal(num * num2) -
                                K * SpecialFunction.cdf_normal(num * num3));
            }

        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double black76(OptionType optType, double F, double K, double T, double sigma)
        {
            var sdev = 0.0;
            if (T > 0.0)
            {
                sdev = sigma * Math.Sqrt(T);
            }
            return black76(optType, F, K, sdev);
        }
    }
}
