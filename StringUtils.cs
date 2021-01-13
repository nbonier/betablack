using System;
using System.Collections.Generic;

namespace betablack
{
    public static class StringUtils
    {
        public static int[] AllIndexesOf(this string str, string substr, bool ignoreCase = false)
        {
            if (string.IsNullOrWhiteSpace(str) || string.IsNullOrWhiteSpace(substr)) return new int[0];
            var list = new List<int>();
            var startIndex = 0;
            while ((startIndex = str.IndexOf(substr, startIndex, ignoreCase ? StringComparison.OrdinalIgnoreCase : StringComparison.Ordinal)) != -1)
                list.Add(startIndex++);
            return list.ToArray();
        }
    }
}