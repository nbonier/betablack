namespace betablack
{    
    public enum OptionType
    {
        Call,
        Put
    }

    public static class __H
    {
        public static double ToEpsilon(this OptionType opt)
        {
            return opt == OptionType.Call ? 1.0 : -1.0;
        }
    }
}
