using System;
using System.Linq;
namespace betablack
{
    public static class FS
    {
        public static Func<R> fun<R>(Func<R> f) { return f; }
        public static Func<A0, R> fun<A0, R>(Func<A0, R> f) { return f; }
        public static Func<A0, A1, R> fun<A0, A1, R>(Func<A0, A1, R> f) { return f; }

        public static Exception FailWith(string s) => new Exception(s);
        public static Exception E_CASE<T>(T p) => new Exception($"Match case error {typeof(T).Name}-{p}");

        public static Exception E_NOTIMPL = new Exception("E_NOTIMPL");

        public static void WithInfo(Exception e, string info)
        {
            var stack = e.StackTrace.Split('\n');
            var myEntry = stack.Last();
            var idx = myEntry.AllIndexesOf(":line").Last();
            var key = myEntry.Substring(0, idx);
            e.Data.Add(key, info); // possible exception: An argument with the same key already exists
        }
    }
}
