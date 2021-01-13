using System;

namespace betablack
{
    class Program
    {
        public static void test_beta_black_perf()
        {
            Console.WriteLine("Computing ... ");
            var N = 100000;
            var tic = DateTime.Now;
            for (var i = 0; i < N; ++i)
            {
                ClosedForm.black76(OptionType.Put, 10.0, 4.0, 0.5, 0.4);
                ClosedForm.black76(OptionType.Put, 10.0, 8.0, 1.5, 0.4);
                ClosedForm.black76(OptionType.Put, 10.0, 4.0, 0.5, 0.4);
                ClosedForm.black76(OptionType.Put, 10.0, 8.0, 1.5, 0.4);
                ClosedForm.black76(OptionType.Put, 10.0, 4.0, 0.5, 0.4);
                ClosedForm.black76(OptionType.Put, 10.0, 8.0, 1.5, 0.4);
            }
            var toc = DateTime.Now;
            var ticBB = DateTime.Now;
            for (var i = 0; i < N; ++i)
            {
                BetaBlack.beta_black(0.1, OptionType.Put, 10.0, 4.0, 0.5, 0.4);
                BetaBlack.beta_black(0.1, OptionType.Put, 8.0, 4.0, 0.5, 0.4);
                BetaBlack.beta_black(0.5, OptionType.Put, 10.0, 4.0, 0.5, 0.4);
                BetaBlack.beta_black(0.5, OptionType.Put, 8.0, 4.0, 0.5, 0.4);
                BetaBlack.beta_black(0.5, OptionType.Put, -10.0, -14.0, 0.5, 0.4);
                BetaBlack.beta_black(0.5, OptionType.Put, -8.0, -10.0, 0.5, 0.4);
            }
            var tocBB = DateTime.Now;

            var ticBBP = DateTime.Now;

            var ignore = 0.0;

            for (var i = 0; i < N; ++i)
            {
                BetaBlack.beta_black_b(0.1, OptionType.Put, 10.0, 4.0, 0.5, 0.4, ref ignore, ref ignore, ref ignore,
                ref ignore, ref ignore, 1.0);
                BetaBlack.beta_black_b(0.1, OptionType.Put, 8.0, 4.0, 0.5, 0.4, ref ignore, ref ignore, ref ignore,
                ref ignore, ref ignore, 1.0);
                BetaBlack.beta_black_b(0.5, OptionType.Put, 10.0, 4.0, 0.5, 0.4, ref ignore, ref ignore, ref ignore,
                ref ignore, ref ignore, 1.0);
                BetaBlack.beta_black_b(0.5, OptionType.Put, 8.0, 4.0, 0.5, 0.4, ref ignore, ref ignore, ref ignore,
                ref ignore, ref ignore, 1.0);
                BetaBlack.beta_black_b(0.5, OptionType.Put, -10.0, -14.0, 0.5, 0.4, ref ignore, ref ignore, ref ignore,
                ref ignore, ref ignore, 1.0);
                BetaBlack.beta_black_b(0.5, OptionType.Put, -8.0, -10.0, 0.5, 0.4, ref ignore, ref ignore, ref ignore,
                ref ignore, ref ignore, 1.0);
            }

            var tocBBP = DateTime.Now;
            Console.WriteLine("Done");
            Console.WriteLine($"Black: {(toc - tic).TotalMilliseconds}");
            Console.WriteLine($"B-Black: {(tocBB - ticBB).TotalMilliseconds}");
            Console.WriteLine($"First Order Risk B-Black : {(tocBBP - ticBBP).TotalMilliseconds}");
            Console.WriteLine($"PV time in # blacks: {(tocBB - ticBB).TotalMilliseconds / (toc - tic).TotalMilliseconds}");
            Console.WriteLine($"Risk first order time in # blacks: {(tocBBP - ticBBP).TotalMilliseconds / (toc - tic).TotalMilliseconds}");

        }

        public static void example()
        {
            var beta = 0.99;
            var f = -2.0;
            var k = -2.0;
            var t = 0.5;
            var sigma = 1.0;
            var callPut = OptionType.Put;
            var bb = new BetaBlackPack(beta, callPut, f, k, t, sigma);
            Console.WriteLine($"{bb.Price}, {bb.dPV_dBeta}");
        }

        static void Main(string[] args)
        {
            example();
            test_beta_black_perf();
        }
    }
}
