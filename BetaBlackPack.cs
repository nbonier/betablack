using System;
namespace betablack
{

    public class BetaBlackPack : IVanillaPack
    {
        public BetaBlackPack(double beta, OptionType callPut, double f, double k, double t, double sigma)
        {
            if (beta < 0 || beta >= 1.0) throw FS.FailWith($"Beta must be in [0, 1.0), got {beta}");
            if (sigma <= 0.0) throw FS.FailWith($"Sigma must be positive, got {sigma}");
            Beta = beta;
            CallPut = callPut;
            F = f;
            K = k;
            T = t;
            Sigma = sigma;

            //! Nota bene: Use AD for first order derivatives and FD on AD for second order

            // 
            // Mu = BetaBlack.inv_m1_beta_black(beta, f, sigma * Math.Sqrt(t));
            double beta_b = 0.0;
            double f_b = 0.0;
            double k_b = 0.0;
            double t_b = 0.0;
            double sigma_b = 0.0;
            Price = BetaBlack.beta_black_b(beta, callPut, f, k, t, sigma, ref beta_b, ref f_b, ref k_b, ref t_b, ref sigma_b, 1.0);
            dPV_dS = f_b;
            dPV_dK = k_b;
            dPV_dT = t_b;
            dPV_dVol = sigma_b;
            dPV_dBeta = beta_b;

            var eps = Math.Sqrt(DoubleUtils.DoubleEpsilon);
            {
                var dF = eps * (1.0 + Math.Abs(f));
                double betaUp_b = 0.0;
                double fUp_b = 0.0;
                double kUp_b = 0.0;
                double tUp_b = 0.0;
                double sigmaUp_b = 0.0;
                double pUp = BetaBlack.beta_black_b(beta, callPut, f + dF, k, t, sigma, ref betaUp_b, ref fUp_b, ref kUp_b, ref tUp_b, ref sigmaUp_b, 1.0);
                double betaDown_b = 0.0;
                double fDown_b = 0.0;
                double kDown_b = 0.0;
                double tDown_b = 0.0;
                double sigmaDown_b = 0.0;
                double pDown = BetaBlack.beta_black_b(beta, callPut, f - dF, k, t, sigma, ref betaDown_b, ref fDown_b, ref kDown_b, ref tDown_b, ref sigmaDown_b, 1.0);
                d2PV_dFdBeta = (betaUp_b - betaDown_b) / (2.0 * dF);
                d2PV_dS2 = (fUp_b - fDown_b) / (2.0 * dF);
                d2PV_dKdS = (kUp_b - kDown_b) / (2.0 * dF);
                d2PV_dSdT = (tUp_b - tDown_b) / (2.0 * dF);
                d2PV_dSdVol = (sigmaUp_b - sigmaDown_b) / (2.0 * dF);
            }

            {
                var dK = eps * (1.0 + Math.Abs(k));
                double betaUp_b = 0.0;
                double fUp_b = 0.0;
                double kUp_b = 0.0;
                double tUp_b = 0.0;
                double sigmaUp_b = 0.0;
                double pUp = BetaBlack.beta_black_b(beta, callPut, f, k + dK, t, sigma, ref betaUp_b, ref fUp_b, ref kUp_b, ref tUp_b, ref sigmaUp_b, 1.0);
                double betaDown_b = 0.0;
                double fDown_b = 0.0;
                double kDown_b = 0.0;
                double tDown_b = 0.0;
                double sigmaDown_b = 0.0;
                double pDown = BetaBlack.beta_black_b(beta, callPut, f, k - dK, t, sigma, ref betaDown_b, ref fDown_b, ref kDown_b, ref tDown_b, ref sigmaDown_b, 1.0);
                d2PV_dKdBeta = (betaUp_b - betaDown_b) / (2.0 * dK);
                // d2PV_dKdS = (fUp_b - fDown_b) / (2.0 * dK);
                d2PV_dK2 = (kUp_b - kDown_b) / (2.0 * dK);
                d2PV_dKdT = (tUp_b - tDown_b) / (2.0 * dK);
                d2PV_dVoldK = (sigmaUp_b - sigmaDown_b) / (2.0 * dK);
            }

            {
                var dSigma = eps * (1.0 + Math.Abs(sigma));
                double betaUp_b = 0.0;
                double fUp_b = 0.0;
                double kUp_b = 0.0;
                double tUp_b = 0.0;
                double sigmaUp_b = 0.0;
                double pUp = BetaBlack.beta_black_b(beta, callPut, f, k, t, sigma + dSigma, ref betaUp_b, ref fUp_b, ref kUp_b, ref tUp_b, ref sigmaUp_b, 1.0);
                double betaDown_b = 0.0;
                double fDown_b = 0.0;
                double kDown_b = 0.0;
                double tDown_b = 0.0;
                double sigmaDown_b = 0.0;
                double pDown = BetaBlack.beta_black_b(beta, callPut, f, k, t, sigma - dSigma, ref betaDown_b, ref fDown_b, ref kDown_b, ref tDown_b, ref sigmaDown_b, 1.0);
                d2PV_dVoldBeta = (betaUp_b - betaDown_b) / (2.0 * dSigma);
                // d2PV_dSdVol = (fUp_b - fDown_b) / (2.0 * dSigma);
                // d2PV_dVoldK = (kUp_b - kDown_b) / (2.0 * dSigma);
                d2PV_dVoldT = (tUp_b - tDown_b) / (2.0 * dSigma);
                d2PV_dVol2 = (sigmaUp_b - sigmaDown_b) / (2.0 * dSigma);
            }
        }

        // this
        public double Beta { get; }
        public OptionType CallPut { get; }
        public double F { get; }
        public double K { get; }
        public double T { get; }
        // public double Mu { get; }
        public double Sigma { get; }
        public double dPV_dBeta { get; }
        public double d2PV_dFdBeta { get; }
        public double d2PV_dKdBeta { get; }
        public double d2PV_dVoldBeta { get; }

        // IVanillaPack
        public double Intrinsic => Math.Max(CallPut.ToEpsilon() * (F - K), 0.0);
        public double Price { get; }
        public double dPV_dS { get; }
        public double dPV_dK { get; }
        public double dPV_dVol { get; }
        public double dPV_dT { get; }
        public double d2PV_dS2 { get; }
        public double d2PV_dSdVol { get; }
        public double d2PV_dVol2 { get; }
        public double d2PV_dK2 { get; }
        public double d2PV_dSdT { get; }
        public double d2PV_dVoldK { get; }
        public double d2PV_dVoldT { get; }
        public double d2PV_dKdT { get; }
        public double d2PV_dKdS { get; }
        public double d3PV_dK2dS => throw FS.E_NOTIMPL;
        public double d3PV_dKdS2 => throw FS.E_NOTIMPL;
        public double d4PV_dK2dS2 => throw FS.E_NOTIMPL;
    }
}
