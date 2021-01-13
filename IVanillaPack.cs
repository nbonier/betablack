namespace betablack
{
    public interface IVanillaPack
    {
        double Intrinsic { get; }
        double Price { get; }
        double dPV_dS { get; }
        double dPV_dK { get; }
        double dPV_dVol { get; }
        double dPV_dT { get; }
        double d2PV_dS2 { get; }
        double d2PV_dSdVol { get; }
        double d2PV_dVol2 { get; }
        double d2PV_dK2 { get; }
        double d2PV_dSdT { get; }
        double d2PV_dVoldK { get; }
        double d2PV_dVoldT { get; }
        double d2PV_dKdT { get; }
        double d2PV_dKdS { get; }
        double d3PV_dK2dS { get; }
        double d3PV_dKdS2 { get; }
        double d4PV_dK2dS2 { get; }
    }
}
